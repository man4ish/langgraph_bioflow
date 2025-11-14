#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_fastq = 'in/non_redundant_fastq/*.gz'
params.reference = 'in/reference/assembly38.fasta'
params.bed = 'in/reference/Twist_ILMN_Exome_2.0_Plus_Panel.hg38.bed'
params.dbsnp = 'in/reference/Homo_sapiens_assembly38.dbsnp138.vcf'

workflow {

    // -----------------------------
    // FastQC
    // -----------------------------
    fastq_pairs = Channel
        .fromFilePairs(params.input_fastq, size: 2, flat: true)
        .map { sample_id, reads -> tuple(sample_id, reads) }

    fastqc_results = fastqc(fastq_pairs)

    // -----------------------------
    // BWA Alignment
    // -----------------------------
    aligned_bams = bwa_align(fastq_pairs)

    // -----------------------------
    // BAM Sorting
    // -----------------------------
    sorted_bams = sort_bam(aligned_bams)

    // -----------------------------
    // Add Read Groups
    // -----------------------------
    rg_bams = add_read_groups(sorted_bams)

    // -----------------------------
    // Mark Duplicates
    // -----------------------------
    dup_bams = mark_duplicates(rg_bams)

    // -----------------------------
    // Index BAM
    // -----------------------------
    indexed_bams = index_bam(dup_bams)

    // -----------------------------
    // BAM Stats
    // -----------------------------
    bam_stats_files = samtools_stats(indexed_bams)
    stats_table_file = bam_stats_table(bam_stats_files)

    // -----------------------------
    // Base Recalibrator
    // -----------------------------
    recal_csv = base_recalibrator(dup_bams)

    // -----------------------------
    // Apply BQSR
    // -----------------------------
    recal_bams = apply_bqsr(dup_bams, recal_csv)

    // -----------------------------
    // Haplotype Caller
    // -----------------------------
    gvcfs = haplotype_caller(recal_bams)

    // -----------------------------
    // Joint Genotyping
    // -----------------------------
    multi_vcf = joint_genotyping(gvcfs)

    // -----------------------------
    // Variant Recalibrator
    // -----------------------------
    recal_files = variant_recalibrator(multi_vcf)

    // -----------------------------
    // Apply VQSR
    // -----------------------------
    apply_vqsr(multi_vcf, recal_files)
}

// =============================
// Processes
// =============================

process fastqc {
    tag { sample_id }
    container 'docker.io/library/exome:latest'
    cpus 8

    input:
    tuple val(sample_id), path(reads)

    output:
    path "out/fastqc/*"

    script:
    """
    mkdir -p out/fastqc
    fastqc -t ${task.cpus} ${reads.join(' ')} --outdir out/fastqc
    """
}

process bwa_align {
    tag { sample_id }
    container 'biocontainers/bwa:v0.7.17-3-deb_cv1'
    cpus 8

    input:
    tuple val(sample_id), path(reads)

    output:
    path "out/aligned_bam/${sample_id}.bam"

    script:
    """
    bwa mem -t ${task.cpus} ${params.reference} ${reads[0]} ${reads[1]} | samtools view -b -o out/aligned_bam/${sample_id}.bam
    """
}

process sort_bam {
    tag { bam_file.baseName }
    container 'docker.io/library/exome:latest'
    cpus 4

    input:
    path bam_file

    output:
    path "out/bamfiles/new_bam/${bam_file.baseName}.sorted.bam"

    script:
    """
    mkdir -p out/bamfiles/new_bam
    samtools sort -@ ${task.cpus} $bam_file > out/bamfiles/new_bam/${bam_file.baseName}.bam.tmp
    mv out/bamfiles/new_bam/${bam_file.baseName}.bam.tmp out/bamfiles/new_bam/${bam_file.baseName}.sorted.bam
    """
}

process add_read_groups {
    tag { bam_file.baseName }
    container 'docker.io/library/exome:latest'
    cpus 4

    input:
    path bam_file

    output:
    path "out/bamfiles/renamed_bam/${bam_file.baseName}.rg.bam"
    path "out/bamfiles/renamed_bam/${bam_file.baseName}.rg.bai"

    script:
    """
    mkdir -p out/bamfiles/renamed_bam
    java -XX:ParallelGCThreads=${task.cpus} -jar /usr/picard/picard.jar AddOrReplaceReadGroups \
        I=$bam_file \
        O=out/bamfiles/renamed_bam/${bam_file.baseName}.rg.bam \
        RGID=${bam_file.baseName} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=${bam_file.baseName} \
        CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT
    """
}

process mark_duplicates {
    tag { bam_file.baseName }
    container 'docker.io/library/exome:latest'
    cpus 8

    input:
    path bam_file

    output:
    path "out/bamfiles/new_bam/${bam_file.baseName}.dup.bam"
    path "out/bamfiles/dup_log2/${bam_file.baseName}_metrics.txt"

    script:
    """
    mkdir -p out/bamfiles/dup_log2
    java -XX:ParallelGCThreads=${task.cpus} -jar /usr/picard/picard.jar MarkDuplicates \
        INPUT=out/bamfiles/renamed_bam/${bam_file.baseName}.rg.bam \
        OUTPUT=out/bamfiles/new_bam/${bam_file.baseName}.dup.bam \
        METRICS_FILE=out/bamfiles/dup_log2/${bam_file.baseName}_metrics.txt \
        VALIDATION_STRINGENCY=LENIENT
    """
}

process index_bam {
    tag { bam_file.baseName }
    container 'docker.io/library/exome:latest'

    input:
    path bam_file

    output:
    path "${bam_file}.bai"

    script:
    """
    mkdir -p out/bamfiles/csv_log
    samtools index $bam_file
    """
}

process samtools_stats {
    tag { bam_file.baseName }
    container 'docker.io/library/exome:latest'
    cpus 8

    input:
    path bam_file

    output:
    path "out/bamfiles/stats/${bam_file.baseName}.stats.txt"

    script:
    """
    mkdir -p out/bamfiles/stats
    samtools stats --threads ${task.cpus} $bam_file > out/bamfiles/stats/${bam_file.baseName}.stats.txt
    """
}

process bam_stats_table {
    tag "bam_stats_table"
    container 'docker.io/library/exome:latest'

    input:
    path stats_files from samtools_stats.out.collect()

    output:
    path "out/bamfiles/stats/table.txt"

    script:
    """
    mkdir -p out/bamfiles/stats
    output_file="out/bamfiles/stats/table.txt"
    echo -e "Sample\\tTotal_Sequences\\tReads_Mapped\\tReads_Mapped_and_Paired" > \$output_file
    for file in ${stats_files.join(' ')}
    do
        sample_name=\$(basename "\$file" .stats.txt)
        total_sequences=\$(awk '/^SN raw total sequences:/ {print \$5}' "\$file")
        reads_mapped=\$(awk '/^SN    reads mapped:/ {print \$4}' "\$file")
        reads_mapped_and_paired=\$(awk '/^SN reads mapped and paired:/ {print \$6}' "\$file")
        echo -e "\$sample_name\\t\$total_sequences\\t\$reads_mapped\\t\$reads_mapped_and_paired" >> \$output_file
    done
    """
}

process base_recalibrator {
    tag { bam_file.baseName }
    container 'docker.io/library/exome:latest'
    cpus 8
    memory '8 GB'

    input:
    path bam_file

    output:
    path "out/bamfiles/renamed_bam/${bam_file.baseName}_recal_data.csv"

    script:
    """
    gatk --java-options '-Xmx4g' BaseRecalibrator \
        -R ${params.reference} \
        -I $bam_file \
        -L ${params.bed} \
        --known-sites ${params.dbsnp} \
        -O out/bamfiles/renamed_bam/${bam_file.baseName}_recal_data.csv
    """
}

process apply_bqsr {
    tag { bam_file.baseName }
    container 'docker.io/library/exome:latest'
    cpus 8
    memory '8 GB'

    input:
    tuple path(bam_file), path(recal_csv)

    output:
    path "out/bamfiles/renamed_bam/${bam_file.baseName}.recal.bam"

    script:
    """
    gatk --java-options '-Xmx4g' ApplyBQSR \
        -R ${params.reference} \
        -I $bam_file \
        -L ${params.bed} \
        --bqsr-recal-file $recal_csv \
        -O out/bamfiles/renamed_bam/${bam_file.baseName}.recal.bam
    """
}

process haplotype_caller {
    tag { bam_file.baseName }
    container 'docker.io/library/exome:latest'
    cpus 8
    memory '16 GB'

    input:
    path bam_file

    output:
    path "out/bamfiles/renamed_bam/${bam_file.baseName}.g.vcf.gz"

    script:
    """
    gatk HaplotypeCaller \
        -R ${params.reference} \
        --native-pair-hmm-threads ${task.cpus} \
        -I $bam_file \
        -L ${params.bed} \
        -O out/bamfiles/renamed_bam/${bam_file.baseName}.g.vcf.gz \
        --dbsnp ${params.dbsnp} \
        -ERC GVCF
    """
}

process joint_genotyping {
    tag "joint_genotyping"
    container 'docker.io/library/exome:latest'
    cpus 8
    memory '16 GB'

    input:
    path gvcfs from haplotype_caller.out.collect()

    output:
    path "out/output.vcf"

    script:
    """
    VARIANTS="\$(for f in ${gvcfs.join(' ')}; do echo -n "--variant \$f "; done)"
    gatk CombineGVCFs -R ${params.reference} \$VARIANTS --dbsnp ${params.dbsnp} -L ${params.bed} -O out/raw_samples.vcf
    gatk GenotypeGVCFs -R ${params.reference} --variant out/raw_samples.vcf --dbsnp ${params.dbsnp} -L ${params.bed} -O out/output.vcf
    """
}

process variant_recalibrator {
    tag "variant_recalibrator"
    container 'docker.io/library/exome:latest'
    cpus 8
    memory '16 GB'

    input:
    path vcf_file from joint_genotyping.out

    output:
    path "out/cohort_snps.recal"
    path "out/cohort_snps.tranches"
    path "out/cohort_snps.plots.R"

    script:
    """
    gatk VariantRecalibrator \
        -R ${params.reference} \
        -V $vcf_file \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${params.dbsnp} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /scratch/skalki/GatkResourceBundle/v0/hapmap_3.3.hg38.vcf.gz \
        --resource:omni,known=false,training=true,truth=true,prior=12.0 /scratch/skalki/GatkResourceBundle/v0/1000G_omni2.5.hg38.vcf.gz \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 /scratch/skalki/GatkResourceBundle/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
        -mode SNP \
        -O out/cohort_snps.recal \
        --tranches-file out/cohort_snps.tranches \
        --rscript-file out/cohort_snps.plots.R
    """
}

process apply_vqsr {
    tag "apply_vqsr"
    container 'docker.io/library/exome:latest'
    cpus 8
    memory '16 GB'

    input:
    path vcf_file from joint_genotyping.out
    path recal_file from variant_recalibrator.out.filter { it.name.endsWith('.recal') }
    path tranches_file from variant_recalibrator.out.filter { it.name.endsWith('.tranches') }

    output:
    path "out/output.vqsr.vcf"

    script:
    """
    gatk ApplyVQSR \
        -R ${params.reference} \
        -V $vcf_file \
        -O out/output.vqsr.vcf \
        --truth-sensitivity-filter-level 99.0 \
        --tranches-file $tranches_file \
        --recal-file $recal_file \
        -mode SNP
    """
}

process snpeff_annotation {
    tag "snpeff_annotation"
    container 'docker.io/library/exome:latest'
    cpus 4
    memory '8 GB'

    input:
    path vcf_file from apply_vqsr.out

    output:
    path "out/output.vqsr.ann.vcf"

    script:
    """
    mkdir -p out
    java -Xmx4g -jar /usr/snpEff/SnpEff.jar GRCh38.86 $vcf_file > out/output.vqsr.ann.vcf
    """
}



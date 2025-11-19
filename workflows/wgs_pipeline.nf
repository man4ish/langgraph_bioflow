#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Workflow: Whole Exome Sequencing (WES) Analysis
 * Docker image: man4ish/exome:latest
 * Steps: FastQC -> Trimming -> BWA -> BQSR -> HaplotypeCaller -> GenotypeGVCFs -> VariantFiltration
 */

params.fastq1        = ""  // Input FASTQ R1
params.fastq2        = ""  // Input FASTQ R2
params.adapter_file  = ""  // Adapter sequences for trimming
params.reference     = ""  // Reference genome
params.known_sites   = ""  // Known variants for BQSR

params.memory_trim   = 4
params.cpu_trim      = 2
params.memory_align  = 8
params.cpu_align     = 4
params.memory_bqsr   = 16
params.cpu_bqsr      = 4
params.memory_hc     = 16
params.cpu_hc        = 4
params.memory_gg     = 16
params.cpu_gg        = 4
params.memory_vf     = 8
params.cpu_vf        = 2

// -----------------------------
// Process: FastQC
// -----------------------------
process FastQC {
    tag "FastQC"
    publishDir "results/fastqc", mode: 'copy'

    input:
    path fastq1
    path fastq2

    output:
    path "*.zip"

    container 'man4ish/exome:latest'
    memory "${params.memory_trim} GB"
    cpus "${params.cpu_trim}"

    script:
    """
    mkdir -p fastqc_results
    fastqc ${fastq1} ${fastq2} --outdir fastqc_results
    mv fastqc_results/*.zip .
    """
}

// -----------------------------
// Process: Trimmomatic
// -----------------------------
process Trimmomatic {
    tag "Trimmomatic"
    publishDir "results/trimmed", mode: 'copy'

    input:
    path fastq1
    path fastq2
    path adapter_file

    output:
    path "trimmed_paired_1.fastq"
    path "trimmed_paired_2.fastq"

    container 'man4ish/exome:latest'
    memory "${params.memory_trim} GB"
    cpus "${params.cpu_trim}"

    script:
    """
    mkdir -p trimmed_results
    trimmomatic PE -phred33 ${fastq1} ${fastq2} \\
        trimmed_results/trimmed_paired_1.fastq trimmed_results/trimmed_unpaired_1.fastq \\
        trimmed_results/trimmed_paired_2.fastq trimmed_results/trimmed_unpaired_2.fastq \\
        ILLUMINACLIP:${adapter_file}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36
    mv trimmed_results/trimmed_paired_1.fastq .
    mv trimmed_results/trimmed_paired_2.fastq .
    """
}

// -----------------------------
// Process: BWA Alignment
// -----------------------------
process BWAAlign {
    tag "BWAAlign"
    publishDir "results/alignment", mode: 'copy'

    input:
    path trimmed1
    path trimmed2
    path reference

    output:
    path "aligned_sorted.bam"
    path "aligned_sorted.bam.bai"

    container 'man4ish/exome:latest'
    memory "${params.memory_align} GB"
    cpus "${params.cpu_align}"

    script:
    """
    bwa mem -t ${task.cpus} ${reference} ${trimmed1} ${trimmed2} | samtools view -Sb - > aligned.bam
    samtools sort -o aligned_sorted.bam aligned.bam
    samtools index aligned_sorted.bam
    """
}

// -----------------------------
// Process: BQSR
// -----------------------------
process BQSR {
    tag "BQSR"
    publishDir "results/bqsr", mode: 'copy'

    input:
    path bam
    path reference
    path known_sites

    output:
    path "recalibrated.bam"

    container 'man4ish/exome:latest'
    memory "${params.memory_bqsr} GB"
    cpus "${params.cpu_bqsr}"

    script:
    """
    gatk BaseRecalibrator -I ${bam} -R ${reference} --known-sites ${known_sites} -O recalibration_report.grp
    gatk ApplyBQSR -R ${reference} -I ${bam} --bqsr-recal-file recalibration_report.grp -O recalibrated.bam
    """
}

// -----------------------------
// Process: HaplotypeCaller
// -----------------------------
process HaplotypeCaller {
    tag "HaplotypeCaller"
    publishDir "results/hc", mode: 'copy'

    input:
    path recal_bam
    path reference

    output:
    path "raw_variants.vcf.gz"

    container 'man4ish/exome:latest'
    memory "${params.memory_hc} GB"
    cpus "${params.cpu_hc}"

    script:
    """
    gatk HaplotypeCaller -R ${reference} -I ${recal_bam} -O raw_variants.vcf.gz -ERC GVCF
    """
}

// -----------------------------
// Process: GenotypeGVCFs
// -----------------------------
process GenotypeGVCFs {
    tag "GenotypeGVCFs"
    publishDir "results/genotyped", mode: 'copy'

    input:
    path gvcf
    path reference

    output:
    path "genotyped_variants.vcf.gz"

    container 'man4ish/exome:latest'
    memory "${params.memory_gg} GB"
    cpus "${params.cpu_gg}"

    script:
    """
    gatk GenotypeGVCFs -R ${reference} -V ${gvcf} -O genotyped_variants.vcf.gz
    """
}

// -----------------------------
// Process: VariantFiltration
// -----------------------------
process VariantFiltration {
    tag "VariantFiltration"
    publishDir "results/filtered", mode: 'copy'

    input:
    path vcf
    path reference

    output:
    path "filtered_variants.vcf.gz"

    container 'man4ish/exome:latest'
    memory "${params.memory_vf} GB"
    cpus "${params.cpu_vf}"

    script:
    """
    gatk VariantFiltration -R ${reference} -V ${vcf} -O filtered_variants.vcf.gz \\
        --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0" --filter-name "filter"
    """
}

// -----------------------------
// Workflow
// -----------------------------
workflow {
    fastqc_out = FastQC(params.fastq1, params.fastq2)
    trimmed = Trimmomatic(params.fastq1, params.fastq2, params.adapter_file)
    aligned = BWAAlign(trimmed.trimmed_paired_1_fastq, trimmed.trimmed_paired_2_fastq, params.reference)
    recal_bam = BQSR(aligned.aligned_sorted_bam, params.reference, params.known_sites)
    raw_vcf = HaplotypeCaller(recal_bam.recalibrated_bam, params.reference)
    genotyped = GenotypeGVCFs(raw_vcf.raw_variants_vcf_gz, params.reference)
    filtered = VariantFiltration(genotyped.genotyped_variants_vcf_gz, params.reference)

    emit:
    filtered.filtered_variants_vcf_gz
}

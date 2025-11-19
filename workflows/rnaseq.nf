#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * RNA-seq Analysis Pipeline
 * Docker image: docker.io/man4ish/rnaseq:latest
 * Steps: Trim -> Kallisto quantification -> STAR alignment -> BAM sort/index -> Picard RNA metrics
 */

params.fastq1              = ""    // Forward FASTQ
params.fastq2              = ""    // Reverse FASTQ
params.adapters            = ""    // Adapter sequences
params.prefix              = "sample"
params.skewer_threads      = 4
params.minimum_read_length = 36

params.kallisto_threads    = 4
params.bootstrap_samples   = 100
params.idx                  = ""   // Kallisto index
params.gtf                  = ""   // GTF annotation

params.STAR_threads         = 8
params.ref_tar              = ""   // STAR reference genome archive

params.ref_flat             = ""   // Picard reference flat file
params.ribosomal_interval   = ""   // Picard ribosomal interval file
params.ref_seq              = ""   // Picard reference sequence

// -----------------------------
// Process: Trim adapters using Skewer
// -----------------------------
process Trim {
    tag "$prefix-trim"
    container 'docker.io/man4ish/rnaseq:latest'
    cpus params.skewer_threads
    memory '4 GB'

    input:
    path fastq1 from params.fastq1
    path fastq2 from params.fastq2
    path adapters from params.adapters
    val prefix from params.prefix
    val threads from params.skewer_threads
    val minlen from params.minimum_read_length

    output:
    path "${prefix}-trimmed-pair1.fastq", emit: trimmed_r1
    path "${prefix}-trimmed-pair2.fastq", emit: trimmed_r2

    script:
    """
    set -exo pipefail
    skewer -t ${threads} -y ${adapters} -l ${minlen} ${fastq1} ${fastq2} -o ${prefix}
    """
}

// -----------------------------
// Process: Kallisto quantification
// -----------------------------
process KallistoQuant {
    tag "Kallisto-quant"
    container 'docker.io/man4ish/rnaseq:latest'
    cpus params.kallisto_threads
    memory '8 GB'

    input:
    path fastq1 from Trim.out.trimmed_r1
    path fastq2 from Trim.out.trimmed_r2
    path idx from params.idx
    path gtf from params.gtf
    val threads from params.kallisto_threads
    val bootstrap from params.bootstrap_samples
    val prefix from params.prefix

    output:
    path "${prefix}_kallisto_output.tar.gz", emit: kallisto_out

    script:
    """
    set -exo pipefail
    mkdir ${prefix}_output
    /software/utils/kallisto quant -t ${threads} -b ${bootstrap} --rf-stranded --genomebam -i ${idx} -g ${gtf} ${fastq1} ${fastq2} -o ${prefix}_output
    tar -czvf ${prefix}_kallisto_output.tar.gz ${prefix}_output
    """
}

// -----------------------------
// Process: STAR alignment
// -----------------------------
process STARAlign {
    tag "STAR-align"
    container 'docker.io/man4ish/rnaseq:latest'
    cpus params.STAR_threads
    memory '16 GB'

    input:
    path fastq1 from Trim.out.trimmed_r1
    path fastq2 from Trim.out.trimmed_r2
    path ref_tar from params.ref_tar
    val threads from params.STAR_threads
    val prefix from params.prefix

    output:
    path "${prefix}_sample.bam", emit: aligned_bam

    script:
    """
    set -exo pipefail
    mkdir ref_bundle
    tar -xzf ${ref_tar} -C ref_bundle --no-same-owner
    ref_path=\$(realpath ref_bundle)
    mv ref_bundle/*/* ref_bundle
    /software/utils/STAR --genomeDir \${ref_path} --runThreadN ${threads} --outSAMtype BAM Unsorted --readFilesIn ${fastq1} ${fastq2} --outFileNamePrefix ${prefix}_sample
    cp ${prefix}_sampleAligned.out.bam ${prefix}_sample.bam
    """
}

// -----------------------------
// Process: Sort and index BAM
// -----------------------------
process SortIndexBam {
    tag "SortIndex-BAM"
    container 'docker.io/man4ish/rnaseq:latest'
    cpus 4
    memory '8 GB'

    input:
    path bam_file from STARAlign.out.aligned_bam
    val prefix from params.prefix

    output:
    path "sorted_${prefix}_sample.bam", emit: sorted_bam
    path "sorted_${prefix}_sample.bam.bai", emit: sorted_bai

    script:
    """
    set -exo pipefail
    samtools sort ${bam_file} -o sorted_${prefix}_sample.bam
    samtools index sorted_${prefix}_sample.bam
    """
}

// -----------------------------
// Process: Picard RNA-seq metrics
// -----------------------------
process PicardMetrics {
    tag "Picard-RNA"
    container 'docker.io/man4ish/rnaseq:latest'
    cpus 2
    memory '4 GB'

    input:
    path bam_file from SortIndexBam.out.sorted_bam
    path ref_flat from params.ref_flat
    path ribosomal_interval from params.ribosomal_interval
    path ref_seq from params.ref_seq
    val prefix from params.prefix

    output:
    path "${prefix}.rna.summary", emit: rna_summary
    path "${prefix}_position.vs.coverage.plot.pdf", emit: coverage_plot

    script:
    """
    set -exo pipefail
    java -jar /software/utils/picard.jar CollectRnaSeqMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS \\
        REF_FLAT=${ref_flat} RIBOSOMAL_INTERVALS=${ribosomal_interval} \\
        STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND \\
        CHART_OUTPUT=${prefix}_position.vs.coverage.plot.pdf \\
        INPUT=${bam_file} OUTPUT=${prefix}.rna.summary \\
        REFERENCE_SEQUENCE=${ref_seq} VALIDATION_STRINGENCY=LENIENT
    """
}

// -----------------------------
// Workflow
// -----------------------------
workflow {
    trimmed = Trim()
    kallisto = KallistoQuant()
    aligned = STARAlign()
    sorted = SortIndexBam()
    metrics = PicardMetrics()

    emit:
    metrics.rna_summary
    metrics.coverage_plot
    kallisto.kallisto_out
    sorted.sorted_bam
    sorted.sorted_bai
}


// nextflow run rnaseq.nf \
//  --fastq1 sample_R1.fastq.gz \
//  --fastq2 sample_R2.fastq.gz \
//  --adapters adapters.fa \
//  --prefix sample \
//  --idx kallisto_index.idx \
//  --gtf annotation.gtf \
//  --ref_tar star_reference.tar.gz \
//  --ref_flat refFlat.txt \
//  --ribosomal_interval rRNA.interval_list \
//  --ref_seq reference.fa

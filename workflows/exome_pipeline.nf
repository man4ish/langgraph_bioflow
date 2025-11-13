nextflow.enable.dsl=2

// Include modules
include { BWA_ALIGN }              from './modules/bwa/main.nf'
include { SAMTOOLS_SORT }          from './modules/samtools/main.nf'
include { PICARD_MARKDUP }         from './modules/picard/main.nf'
include { GATK_HAPLOTYPECALLER }   from './modules/gatk/main.nf'
include { SNPEFF_ANNOTATE }        from './modules/snpeff/main.nf' // optional if annotation is next

workflow exome_pipeline {

    // Channel for paired FASTQ reads (FastQC already done)
    reads_ch = Channel.fromFilePairs(params.reads, flat: true)
    reference = file(params.reference)

    log.info "📁 Using FASTQ input: ${params.reads}"
    log.info "📁 Using reference genome: ${params.reference}"
    log.info "📁 Skipping FastQC (already in results/fastqc)"

    /*
     * STEP 1: Alignment with BWA
     */
    aligned_sams = reads_ch.map { sample_id, reads ->
        tuple(sample_id, reads, reference)
    } | BWA_ALIGN

    /*
     * STEP 2: Sort BAM files using SAMTOOLS
     */
    sorted_bams = SAMTOOLS_SORT(aligned_sams)

    /*
     * STEP 3: Mark duplicates using PICARD
     */
    dedup_bams = PICARD_MARKDUP(sorted_bams)

    /*
     * STEP 4: Variant calling with GATK HaplotypeCaller
     * ❗ Only pass the dedup_bams channel, not the reference separately
     */
    vcfs = GATK_HAPLOTYPECALLER(dedup_bams.map { sample_id, bam -> tuple(sample_id, bam, reference) })


    /*
     * (Optional) STEP 5: Annotation with SnpEff
     */
    annotated_vcfs = SNPEFF_ANNOTATE(vcfs)

    annotated_vcfs.view { vcf -> "✅ Annotated VCF generated: ${vcf}" }
}

workflow {
    exome_pipeline()
}

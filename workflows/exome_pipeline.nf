nextflow.enable.dsl=2

include { FASTQC } from './modules/fastqc/main.nf'
include { SAMTOOLS_SORT } from './modules/samtools/main.nf'
include { PICARD_MARKDUP } from './modules/picard/main.nf'
include { GATK_HAPLOTYPECALLER } from './modules/gatk/main.nf'

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, flat: true)

    qc_results = reads_ch | map { sample_id, reads -> tuple(sample_id, reads) }
    sorted_bams = SAMTOOLS_SORT(qc_results)
    dedup_bams  = PICARD_MARKDUP(sorted_bams)
    vcfs        = GATK_HAPLOTYPECALLER(dedup_bams)

    vcfs.view { it -> "VCF generated: ${it}" }
}

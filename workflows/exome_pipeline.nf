nextflow.enable.dsl=2

include { SAMTOOLS_SORT } from '../modules/samtools/main.nf'
include { PICARD_MARKDUP } from '../modules/picard/main.nf'
include { GATK_HAPLOTYPECALLER } from '../modules/gatk/main.nf'

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, flat: true)

    sorted_bams = SAMTOOLS_SORT(reads_ch)
    dedup_bams  = PICARD_MARKDUP(sorted_bams)
    vcfs        = GATK_HAPLOTYPECALLER(dedup_bams)

    vcfs.view { vcf -> "VCF generated: ${vcf}" }
}

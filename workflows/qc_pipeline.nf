nextflow.enable.dsl=2

include { FASTQC } from '../modules/fastqc/main.nf'

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, flat: true)
    qc_results = FASTQC(reads_ch)

    qc_results.view { file -> "QC output generated: ${file}" }
}

nextflow.enable.dsl=2

include { FASTQC } from './modules/fastqc/main.nf'

workflow qc_pipeline {
    Channel
        .fromFilePairs("data/fastq/*_{R1,R2}_001.fastq.gz", flat: true)
        .set { reads_ch }

    qc_results = FASTQC(reads_ch)

    qc_results.view { it -> "✅ QC output generated: ${it}" }
}

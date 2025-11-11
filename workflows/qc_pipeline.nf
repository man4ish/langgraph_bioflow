#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Minimal Nextflow Test
 * ----------------------
 * Confirms Nextflow is working correctly in the langgraph_bioflow project.
 */

process HELLO_QC {
    input:
    path sample_file

    output:
    path "qc_summary.txt"

    script:
    """
    echo "Running mock QC for file: ${sample_file}" > qc_summary.txt
    echo "Mean_Q30=0.92" >> qc_summary.txt
    echo "Duplication_Rate=0.10" >> qc_summary.txt
    echo "Alignment_Rate=0.87" >> qc_summary.txt
    """
}

workflow {
    samples_ch = Channel.fromPath("samples/*.fastq")

    qc_results = HELLO_QC(samples_ch)

    qc_results.view { file ->
        println "✅ QC summary generated: ${file}"
    }
}

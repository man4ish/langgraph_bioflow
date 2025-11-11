#!/usr/bin/env nextflow

/*
 * Mock QC pipeline for LangGraph BioFlow testing
 * Works with Nextflow 25.x (DSL2)
 */

nextflow.enable.dsl = 2

params.input = "${baseDir}/samples/sample1.fastq"
params.output_dir = "${baseDir}/results"

workflow {
    main:
        Channel
            .fromPath(params.input)
            .set { fastq_files }

        QC_MOCK(fastq_files)

        QC_MOCK.out.view { result ->
            println "Output file inside sandbox: ${result}"

            def outdir = file(params.output_dir)
            outdir.mkdirs()

            // Copy file to results folder using Groovy NIO
            result.each { f ->
                java.nio.file.Files.copy(
                    f.toPath(),
                    outdir.resolve(f.getName()),
                    java.nio.file.StandardCopyOption.REPLACE_EXISTING
                )
            }

            println "File copied to: ${outdir}/qc_summary.txt"
        }
}

process QC_MOCK {
    tag "$fastq.simpleName"

    input:
        path fastq

    output:
        path "qc_summary.txt"

    script:
    """
    echo "Running mock QC for file: $fastq" > qc_summary.txt
    echo "Mean_Q30=0.92" >> qc_summary.txt
    echo "Duplication_Rate=0.10" >> qc_summary.txt
    echo "Alignment_Rate=0.87" >> qc_summary.txt
    """
}

process FASTQC {
    tag "${sample_id}"
    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.zip"), path("*.html")

    script:
    """
    fastqc --threads 4 --outdir ./ ${reads}
    """
}


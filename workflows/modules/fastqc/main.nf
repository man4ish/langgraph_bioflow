process FASTQC {
    tag "$sample_id"
    publishDir "results/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads1), path(reads2)

    output:
    tuple val(sample_id), path("*.html"), path("*.zip")

    script:
    """
    fastqc --threads 4 --outdir ./ ${reads1} ${reads2}
    """
}

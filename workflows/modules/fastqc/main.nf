process FASTQC {
    tag "${sample_id}"
    container 'biocontainers/fastqc:v0.12.1_cv8'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.html"
    path "*.zip"

    script:
    """
    fastqc ${reads} -o ./
    """
}

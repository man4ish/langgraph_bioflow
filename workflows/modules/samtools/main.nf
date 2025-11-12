process SAMTOOLS_SORT {
    tag "${sample_id}"
    container 'biocontainers/samtools:v1.19.2_cv1'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    path "${sample_id}.sorted.bam"

    script:
    """
    samtools sort -o ${sample_id}.sorted.bam ${bam_file}
    """
}

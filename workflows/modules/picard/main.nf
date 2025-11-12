process PICARD_MARKDUP {
    tag "${sample_id}"
    container 'broadinstitute/picard:2.27.5'

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    path "${sample_id}.dedup.bam"
    path "${sample_id}.metrics.txt"

    script:
    """
    picard MarkDuplicates \
        I=${sorted_bam} \
        O=${sample_id}.dedup.bam \
        M=${sample_id}.metrics.txt \
        REMOVE_DUPLICATES=true
    """
}

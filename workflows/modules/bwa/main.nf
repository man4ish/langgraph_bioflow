process BWA_ALIGN {
    tag "$sample_id"
    publishDir "results/bwa", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    // Remove container line

    script:
    """
    bwa mem -t 4 ${reference} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    """
}



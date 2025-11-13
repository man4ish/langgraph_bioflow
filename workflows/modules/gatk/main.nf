process GATK_HAPLOTYPECALLER {
process GATK_HAPLOTYPECALLER {
    tag "${sample_id}"
    container 'broadinstitute/gatk:4.4.0.0'

    input:
    tuple val(sample_id), path(bam), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller \
      -R ${reference} \
      -I ${bam} \
      -O ${sample_id}.vcf.gz
    """
}

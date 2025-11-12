process GATK_HAPLOTYPECALLER {
    tag "${sample_id}"
    container 'broadinstitute/gatk:4.4.0.0'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz")

    script:
    """
    echo "Simulating GATK HaplotypeCaller for sample: ${sample_id}"
    echo "VCF_HEADER" > ${sample_id}.vcf
    echo "${sample_id}\tchr1\t12345\tA\tG" >> ${sample_id}.vcf
    bgzip ${sample_id}.vcf
    """
}

process SNPEFF_ANNOTATION {

    tag "$sample_id"
    publishDir "${params.outdir}/snpeff", mode: 'copy'

    container 'biocontainers/snpeff:v5.2.1_cv1'

    input:
    tuple val(sample_id), path(vcf_file)
    val genome

    output:
    tuple val(sample_id), path("${sample_id}.snpeff.vcf"), path("${sample_id}.snpeff.summary.txt")

    script:
    """
    echo "Running SnpEff annotation for sample: ${sample_id}"
    
    snpEff -Xmx8g -v ${genome} ${vcf_file} > ${sample_id}.snpeff.vcf

    # Generate simplified summary for downstream AI modules (RAG input)
    grep -v '^#' ${sample_id}.snpeff.vcf | \
        awk -F'\t' '{ split(\$8,a,\"|\"); print \$1\"\\t\"\$2\"\\t\"a[1]\"\\t\"a[2]\"\\t\"a[3] }' \
        > ${sample_id}.snpeff.summary.txt

    echo "SnpEff annotation completed for ${sample_id}"
    """
}

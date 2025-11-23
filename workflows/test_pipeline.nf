#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_fastq = 'data/fastq/*.gz'
params.reference = 'data/reference/assembly38.fasta'
params.bed = 'data/reference/Twist_ILMN_Exome_2.0_Plus_Panel.hg38.bed'
params.dbsnp = 'data/reference/Homo_sapiens_assembly38.dbsnp138.vcf'


process bwa_align {
    tag { sample_id }
    container 'biocontainers/bwa:v0.7.17-3-deb_cv1'
    cpus 8

    input:
    tuple val(sample_id), path(reads)

    output:
    path "out/aligned_bam/${sample_id}.bam"

    script:
    """
    bwa mem -t ${task.cpus} ${params.reference} ${reads[0]} ${reads[1]} | samtools view -b -o out/aligned_bam/${sample_id}.bam
    """
}

workflow {
    aligned_bams = bwa_align(fastq_pairs)
}



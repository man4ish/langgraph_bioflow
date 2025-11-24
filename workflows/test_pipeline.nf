#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reference = "data/reference/assembly38.fasta"

Channel.fromFilePairs("data/fastq/*_{R1,R2}_001.fastq.gz", size: 2)
       .set { fastq_pairs }

process bwa_align {
    tag { sample_id }
    container 'man4ish/bwa-arm64-native:latest'
    cpus 8
    memory '32 GB'
    publishDir 'results/aligned_sam', mode: 'copy', overwrite: true

    input:
    tuple val(sample_id), path(reads)

    output:
    path "${sample_id}.sam"

    script:
    """
    echo "=== Reference and index files in container ==="
    ls -la assembly38.fasta*
    
    echo "=== Starting BWA MEM for $sample_id ==="
    bwa mem -t ${task.cpus} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA' \
        assembly38.fasta ${reads[0]} ${reads[1]} > ${sample_id}.sam
    
    echo "=== DONE === SAM file size and lines:"
    ls -lh ${sample_id}.sam
    grep -v "^@" ${sample_id}.sam | wc -l | awk '{print \$1 " alignment lines"}'
    """
}

workflow {
    fastq_pairs | bwa_align
}

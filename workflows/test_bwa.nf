nextflow.enable.dsl=2

// Include the BWA module
include { BWA_ALIGN } from './modules/bwa/main.nf'

workflow test_bwa {

    // Input channel: paired-end reads
    reads_ch = Channel.fromFilePairs('data/fastq/*_{R1,R2}_001.fastq.gz', flat: true)

    // Reference genome
    reference = file('data/reference/genome.fa')

    log.info "🧬 Testing BWA alignment..."
    log.info "📁 FASTQ input: data/fastq/*_{R1,R2}_001.fastq.gz"
    log.info "📁 Reference genome: data/reference/genome.fa"

    // Prepare input as tuples of (sample_id, R1, R2, reference)
    aligned_input = reads_ch.map { sample_id, reads ->
        tuple(sample_id, reads[0], reads[1], reference)
    }

    // Run BWA
    aligned_bams = BWA_ALIGN(aligned_input)

    // View output
    aligned_bams.view { "✅ BWA alignment complete for: ${it[0]}" }
}

workflow {
    test_bwa()
}

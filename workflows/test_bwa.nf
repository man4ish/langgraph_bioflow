nextflow.enable.dsl=2

include { BWA_ALIGN } from './modules/bwa/main.nf'

workflow {

    log.info "📁 Using FASTQ input: ${params.reads}"
    log.info "📁 Using reference genome: ${params.reference}"

    // 1) Load read pairs
    reads_ch = Channel.fromFilePairs(
                    params.reads,
                    flat: true)

    // 2) Wrap every input in a tuple explicitly
    aligned_inputs = reads_ch.map { pair ->
        def sample_id = pair[0]
        def reads     = pair[1]
        def ref       = file(params.reference)

        return tuple(sample_id, reads, ref)
    }

    // 3) Run BWA
    result = BWA_ALIGN(aligned_inputs)

    result.view()
}

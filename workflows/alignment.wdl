version 1.0

workflow alignment {
  input {
    File reads
  }

  call align_task {
    input: reads = reads
  }

  output {
    File aligned_bam = align_task.output_bam
  }
}

task align_task {
  input {
    File reads
  }

  command <<<
    echo "Simulating alignment for ~{reads}" > aligned.bam
  >>>

  output {
    File output_bam = "aligned.bam"
  }

  runtime {
    docker: "biocontainers/bwa:v0.7.17_cv1"
  }
}

#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ===========================================
// PARAMETERS (can be overridden via CLI)
// ===========================================
params.sid           = ""                       // Sample ID
params.fastq_tar     = ""                       // Tarball of FASTQ files
params.ref_tar       = ""                       // Tarball of reference genome
params.nosecondary   = false                    // Skip secondary analysis
params.localcores    = null                     // Optional local cores
params.localmem      = null                     // Optional local memory
params.localvmem     = null                     // Optional local vmem
params.memory_count  = 32                       // Task memory (GB)
params.cpu_count     = 8                        // Task CPU
params.storage_count = 100                      // Task disk (GB)
params.aggr_id       = ""                       // Aggregation ID
params.data_tar      = ""                       // Tarball of aggregated data
params.aggr_csv_file = ""                       // CSV listing sample paths
params.norm_type     = null                     // Optional normalization
params.memory_aggr   = 32
params.cpu_aggr      = 8
params.storage_aggr  = 100

// ===========================================
// CHANNELS
// ===========================================
Channel.fromPath(params.fastq_tar).set { fastq_tar_ch }
Channel.fromPath(params.ref_tar).set { ref_tar_ch }
Channel.fromPath(params.data_tar).set { data_tar_ch }
Channel.fromPath(params.aggr_csv_file).set { aggr_csv_ch }

// ===========================================
// PROCESS: Count
// ===========================================
process CELLRANGER_COUNT {

    tag "${params.sid}"

    input:
    path fastq_tar from fastq_tar_ch
    path ref_tar   from ref_tar_ch

    output:
    path "${params.sid}_count_result.tar.gz"

    script:
    """
    set -exo pipefail

    mkdir fastq_bundle
    tar -xzf $fastq_tar -C fastq_bundle --no-same-owner
    mv fastq_bundle/*/* fastq_bundle
    fastq_path=\$(realpath fastq_bundle)

    mkdir ref_bundle
    tar -xzf $ref_tar -C ref_bundle --no-same-owner
    mv ref_bundle/*/* ref_bundle
    ref_path=\$(realpath ref_bundle)

    cmd="count --id=${params.sid} --fastqs=\$fastq_path --transcriptome=\$ref_path"

    if [ "${params.nosecondary}" = "true" ]; then
        cmd+=" --nosecondary"
    fi
    if [ ! -z "${params.localcores}" ]; then
        cmd+=" --localcores ${params.localcores}"
    fi
    if [ ! -z "${params.localmem}" ]; then
        cmd+=" --localmem ${params.localmem}"
    fi
    if [ ! -z "${params.localvmem}" ]; then
        cmd+=" --localvmem ${params.localvmem}"
    fi

    /software/reboot-utils/cellranger/bin/cellranger \$cmd

    tar -czvf ${params.sid}_count_result.tar.gz "\$ref_path/${params.sid}"
    """
    cpus params.cpu_count
    memory "${params.memory_count} GB"
    time '24h'
    container 'docker.io/man4ish/cellranger:latest'
}

// ===========================================
// PROCESS: Aggregate
// ===========================================
process CELLRANGER_AGGR {

    tag "${params.aggr_id}"

    input:
    path data_tar from data_tar_ch
    path aggr_csv_file from aggr_csv_ch

    output:
    path "${params.aggr_id}_result.tar.gz"

    script:
    """
    set -exo pipefail

    mkdir data_bundle
    tar -xzf $data_tar -C data_bundle --no-same-owner
    mv data_bundle/*/* data_bundle
    data_path=\$(realpath data_bundle)

    while IFS= read -r line; do
        if [[ "\$line" == "sample_id"* ]]; then
            echo "\$line" >> ${aggr_csv_file}.tmp
        else
            A="\$(cut -d',' -f1 <<<"\$line")"
            echo "\$A,\$data_path/\$A/outs/molecule_info.h5" >> ${aggr_csv_file}.tmp
        fi
    done < ${aggr_csv_file}

    mv ${aggr_csv_file}.tmp ${aggr_csv_file}

    cmd="aggr --id=${params.aggr_id} --csv=${aggr_csv_file}"
    if [ ! -z "${params.norm_type}" ]; then
        cmd+=" --norm_type ${params.norm_type}"
    fi

    /software/reboot-utils/cellranger/bin/cellranger \$cmd

    tar -czvf ${params.aggr_id}_result.tar.gz "\$data_path/${params.aggr_id}"
    """
    cpus params.cpu_aggr
    memory "${params.memory_aggr} GB"
    time '24h'
    container 'docker.io/man4ish/cellranger:latest'
}

// ===========================================
// WORKFLOW
// ===========================================
workflow {

    CELLRANGER_COUNT()
    CELLRANGER_AGGR()
}

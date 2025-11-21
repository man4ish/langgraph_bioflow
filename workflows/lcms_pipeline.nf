#!/usr/bin/env nextflow

/*
 * proteomics_analysis.nf
 * Nextflow pipeline for LC-MS proteomics analysis
 * Steps: QC -> Differential Expression -> Heatmap/Venn Visualization
 */

nextflow.enable.dsl=2

params.input = "$PWD/data/*.tsv"   // Input proteomics TSV files
params.outdir = "$PWD/out"         // Output directory
params.container = "man4ish/lcms:latest"  // Docker image

// -----------------------------
// QC process
// -----------------------------
process QC {
    tag { file.getName() }
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path file from params.input

    output:
    path "preprocessed.tsv" into qc_out

    container params.container

    script:
    """
    python3 /software/lcms_analysis/1_qc.py --input ${file} --output preprocessed.tsv
    """
}

// -----------------------------
// Differential Expression
// -----------------------------
process DiffExp {
    tag { file.getName() }
    publishDir "${params.outdir}/diff_exp", mode: 'copy'

    input:
    path qc_file from qc_out

    output:
    path "diff_exp.tsv" into diff_out

    container params.container

    script:
    """
    python3 /software/lcms_analysis/2_diff_exp.py --input ${qc_file} --output diff_exp.tsv
    """
}

// -----------------------------
// Heatmap and Venn Visualization
// -----------------------------
process HeatmapVenn {
    tag { file.getName() }
    publishDir "${params.outdir}/figures", mode: 'copy'

    input:
    path diff_file from diff_out

    output:
    path "heatmap.png"
    path "venn.png"

    container params.container

    script:
    """
    mkdir -p figures
    python3 /software/lcms_analysis/3_heatmap_venn.py --input ${diff_file} --output figures
    """
}

// -----------------------------
// Workflow definition
// -----------------------------
workflow {
    QC()
    DiffExp()
    HeatmapVenn()
}

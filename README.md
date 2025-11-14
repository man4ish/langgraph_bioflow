# LangGraph BioFlow

**AI-Orchestrated Bioinformatics Workflow Supervisor (Nextflow Edition)**

---

## Overview

**LangGraph BioFlow** is an experimental orchestration framework that integrates **LangGraph** with **Nextflow** to build an **intelligent, adaptive bioinformatics workflow manager**.

It introduces **AI-assisted decision-making** into traditional bioinformatics pipelines — dynamically deciding whether to proceed, rerun, or halt based on QC outcomes.

The goal is to evolve into a **smart pipeline supervisor** capable of:

* Monitoring and controlling multi-stage computational pipelines
* Evaluating data quality in real time
* Making dynamic rerun or stop decisions based on QC metrics
* Generating concise, interpretable summaries for researchers

---

## Current Capabilities

✅ **Sequential execution**: FastQC → Alignment → Deduplication → Variant Calling
✅ **LangGraph-based control flow** with retry limits and state memory
✅ **QC-driven branching** — automatically re-runs QC if metrics fail thresholds
✅ **Nextflow integration** for modular process execution
✅ **Scalable design** — can extend to RNA-seq, GWAS, or WGS workflows

---

## Architecture

```
Raw Reads
   ↓
FASTQC → QC Decision (LangGraph)
   ├── PASS → Exome Pipeline (SAMTOOLS_SORT → PICARD_MARKDUP → GATK_HAPLOTYPECALLER)
   └── FAIL → Re-run QC (max 3 attempts)
```

Each stage is represented as a **LangGraph node**, which calls a corresponding **Nextflow workflow**.
The **QC Decision Node** analyzes the results of `qc_pipeline.nf` and decides the next action.

---

## Example Workflow Execution

```
project_root/
│
├── workflows/                      # Nextflow pipelines
│   ├── qc_pipeline.nf
│   ├── exome_pipeline.nf
│   └── __init__.py
│
├── src/                            # Python scripts
│   ├── __init__.py
│   │
│   ├── pipeline/                   # LangGraph workflow orchestration
│   │   ├── langgraph_workflow.py   # Full workflow: QC → Exome → Pathway → RAG → GNN → Viz
│   │   └── __init__.py
│   │
│   ├── data/                       # Data extraction / graph building / features
│   │   ├── __init__.py
│   │   ├── neo4j_extractor.py      # Extract nodes/edges from Neo4j KG
│   │   ├── graph_builder.py        # Build PyTorch Geometric graphs
│   │   └── feature_engineering.py  # Node/edge features
│   │
│   ├── analysis/                   # Domain-specific analysis
│   │   ├── pathway_enrichment_gene_list.py
│   │   └── __init__.py
│   │
│   ├── rag/                        # RAG + Neo4j KG scripts
│   │   ├── rag_pipeline.py         # LLM + Neo4j structured data generation
│   │   └── structured_drug_kg.py  # Knowledge graph builder
│   │
│   ├── gnn/                        # GNN prediction
│   │   ├── gnn_predictor.py
│   │   └── __init__.py
│   │
│   └── utils/                      # Utility functions (logging, config, etc.)
│       ├── __init__.py
│       └── helpers.py
│
├── out/                             # Output data, reports, visualization
│   ├── bamfiles/
│   ├── figures/
│   └── tables/
│
├── logs/                            # Logs for workflow, GNN, RAG
│
├── requirements.txt                 # Python dependencies
├── environment.yml                  # Conda environment (optional)
└── README.md

```

Example output:

```
=== Running QC (attempt 1) ===
QC Decision: FAIL
QC failed on attempt 1. Retrying with adjusted parameters...

=== Running QC (attempt 2) ===
QC Decision: PASS
Proceeding to exome pipeline.

=== Running Exome Pipeline ===
VCF generated: results/gatk/sample.vcf

Final Workflow State:
decision='COMPLETE' step='exome' attempt=2
```

---

## Project Goals

* Extend to end-to-end workflows (Alignment → QC → Imputation → GWAS → Reporting)
* Replace rule-based QC thresholds with **LLM-driven evaluation**
* Integrate **LangGraph memory and persistence** for auditability
* Add **Nextflow monitoring hooks** for runtime introspection
* Develop a **FastAPI/Streamlit dashboard** for live graph visualization

---

## Directory Layout

```
langgraph_bioflow/
├── main.py                     # LangGraph controller orchestrating Nextflow workflows
├── workflows/
│   ├── qc_pipeline.nf          # Runs FastQC only
│   └── exome_pipeline.nf       # Runs sort → markdup → variant calling
├── modules/
│   ├── fastqc/main.nf
│   ├── samtools/main.nf
│   ├── picard/main.nf
│   └── gatk/main.nf
├── utils/
│   ├── logger.py
│   └── helpers.py
├── config/
│   └── settings.yaml
├── requirements.txt
└── README.md
```

---

## Quick Start

### 1. Clone and Install Dependencies

```bash
git clone https://github.com/yourusername/langgraph_bioflow.git
cd langgraph_bioflow
pip install -r requirements.txt
```

### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

### 3. Run the Orchestrated Workflow

```bash
python main.py
```

---

## Future Enhancements

* [ ] Parse **real FastQC results** instead of random PASS/FAIL simulation
* [ ] Implement **LLM-driven QC reasoning** (LangGraph + GPT)
* [ ] Introduce **Nextflow event hooks** for progress tracking
* [ ] Store LangGraph state in **SQLite** for persistent runs
* [ ] Add **FastAPI dashboard** for interactive graph visualization
* [ ] Extend for **RNA-seq, metagenomics, and multi-omics pipelines**

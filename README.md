# LangGraph BioFlow

**AI-Orchestrated Bioinformatics Workflow Supervisor (Nextflow + RAG + GNN Edition)**

---

## Overview

LangGraph BioFlow is an experimental orchestration framework that integrates **LangGraph**, **Nextflow**, and **AI/ML modules** to build an intelligent, adaptive bioinformatics workflow manager.

It introduces **AI-assisted decision-making** into traditional bioinformatics pipelines — dynamically deciding which pipeline to run, whether to proceed, rerun, or halt based on QC outcomes, pathway enrichment, or predicted drug-target interactions.

The goal is to evolve into a **smart pipeline supervisor** capable of:

* Monitoring and controlling multi-stage computational pipelines
* Evaluating data quality and gene/pathway significance in real time
* Making dynamic rerun or stop decisions based on QC, enrichment, or prediction scores
* Generating interpretable visual summaries and drug predictions

---

## Current Capabilities

✅ Sequential execution: FastQC → Alignment/Dedup → Variant Calling (Exome/WGS/RNA-Seq)

✅ LangGraph-based control flow with retry limits and persistent state memory

✅ QC-driven branching — automatically retries QC (max 3 attempts)

✅ Pathway enrichment from SNPeff gene lists

✅ Retrieval-Augmented Generation (RAG) querying the Drug–Target–Cancer knowledge graph

✅ GNN predictor for drug ranking and scoring

✅ Visualization of pathways, networks, and predicted drugs

✅ Modular and scalable design — easily extendable to multi-omics workflows

---

## Architecture

```
Raw Reads
   ↓
FASTQC → QC Decision (LangGraph)
   ├── PASS → Choose Pipeline: WGS / Exome / RNA-Seq → Variant Calling
   └── FAIL → Re-run QC (max 3 attempts)
         ↓
Pathway Enrichment → RAG → GNN Predictor → Visualization
```

Each stage is represented as a **LangGraph node**, which calls either a Nextflow workflow or a Python module.

---

### Example Workflow Execution

```
=== Running QC (attempt 1) ===
QC Decision: FAIL
QC failed on attempt 1. Retrying with adjusted parameters...

=== Running QC (attempt 2) ===
QC Decision: PASS
Proceeding to Exome pipeline.

=== Running Exome Pipeline ===
VCF generated: results/gatk/sample.vcf

=== Pathway Enrichment ===
Gene list: ["TP53", "APP", "SNCA", "LRRK2"]

=== RAG Pipeline ===
Neo4j KG queried: 100 nodes, 150 edges

=== GNN Predictor ===
Drug Predictions: ["DrugA", "DrugB"], Scores: [0.95, 0.88]

=== Visualization ===
Figures saved to out/figures/

Final Workflow State:
decision='COMPLETE' step='exome' attempt=2
```

---

## Project Goals

* Extend to end-to-end workflows (QC → Alignment → Variant Calling → Pathway → RAG → GNN → Reporting)
* Replace rule-based QC thresholds with **LLM-driven evaluation**
* Integrate **GNN predictor** for drug-target discovery
* Store **LangGraph state** for auditability and reproducibility
* Add **visualization** of pathway enrichment, network, and drug predictions
* Build **FastAPI/Streamlit dashboard** for live workflow graph and prediction results
* Support multi-omics pipelines: RNA-Seq, WGS, metagenomics

---

## Directory Layout

```
project_root/
│
├── workflows/                      # Nextflow pipelines
│   ├── qc_pipeline.nf
│   ├── exome_pipeline.nf
│   ├── wgs_pipeline.nf
│   ├── rnaseq_pipeline.nf
│   └── __init__.py
│
├── src/                            # Python modules
│   ├── __init__.py
│   │
│   ├── pipeline/                   # LangGraph workflow orchestration
│   │   ├── langgraph_workflow.py   # Full workflow: QC → Pipeline → Pathway → RAG → GNN → Viz
│   │   └── __init__.py
│   │
│   ├── data/                       # Data extraction / graph building / features
│   │   ├── __init__.py
│   │   ├── neo4j_extractor.py
│   │   ├── graph_builder.py
│   │   └── feature_engineering.py
│   │
│   ├── analysis/                   # Domain-specific analysis
│   │   ├── pathway_enrichment_gene_list.py
│   │   └── __init__.py
│   │
│   ├── rag/                        # RAG + Neo4j KG scripts
│   │   ├── rag_pipeline.py
│   │   └── structured_drug_kg.py
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

---

## Quick Start

1. **Clone and Install Dependencies**

```bash
git clone https://github.com/yourusername/langgraph_bioflow.git
cd langgraph_bioflow
pip install -r requirements.txt
```

2. **Install Nextflow**

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

3. **Run the Orchestrated Workflow**

```bash
python main.py
```

> Note: Ensure Neo4j is running for RAG + GNN modules.

---

## Future Enhancements

* Parse real FastQC results instead of simulation
* Implement **LLM-driven QC reasoning**
* Introduce Nextflow event hooks for progress tracking
* Store LangGraph state in **SQLite** for persistent runs
* Add **interactive dashboard** (FastAPI/Streamlit) for workflow graph and predictions
* Extend to RNA-Seq, metagenomics, and multi-omics pipelines
* Integrate more advanced GNN architectures for novel drug prediction

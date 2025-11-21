# LangGraph BioFlow

**AI-Orchestrated Bioinformatics Workflow Supervisor (Nextflow + RAG + GNN Edition)**

---

## Overview

LangGraph BioFlow is an experimental orchestration framework that integrates **LangGraph**, **Nextflow**, and **AI/ML modules** to build an intelligent, adaptive bioinformatics workflow manager.

It introduces **AI-assisted decision-making** into bioinformatics pipelines — dynamically deciding which pipeline to run, whether to proceed, rerun, or halt based on QC outcomes, pathway enrichment, or predicted drug-target interactions.

The goal is to evolve into a **smart pipeline supervisor** capable of:

* Monitoring and controlling multi-stage computational pipelines
* Evaluating data quality and gene/pathway significance in real time
* Making dynamic rerun or stop decisions based on QC, enrichment, or prediction scores
* Generating interpretable visual summaries, drug predictions, and network outputs

---

## Current Capabilities

✅ Sequential execution: FastQC → Alignment/Dedup → Variant Calling (Exome/WGS/RNA-Seq)  
✅ LangGraph-based control flow with retry limits and persistent state memory  
✅ QC-driven branching — automatically retries QC (max 3 attempts)  
✅ Pathway enrichment from SNPeff gene lists (`src/analysis/pathway_enrichment_gene_list.py`)  
✅ Gene → Pathway network visualization (`src/analysis/gene_pathway_network_visualization_rag.py`)  
✅ Retrieval-Augmented Generation (RAG) querying Drug–Target–Cancer knowledge graph (`src/rag`)  
✅ GNN predictor for drug ranking and scoring (`src/gnn`)  
✅ Visualization of pathways, networks, and predicted drugs  
✅ Modular and scalable design — easily extendable to multi-omics workflows

---

## Architecture

![Multi-omics Architecture Diagram](https://raw.githubusercontent.com/man4ish/langgraph_bioflow/dev-nextflow-modular/images/multi-omics_architecture_diagram.png)


```

Raw Reads
↓
FASTQC → QC Decision (LangGraph)
├── PASS → Choose Pipeline → Variant Calling → Pathway Enrichment
│                                      ↓
│                               Gene → Pathway Network
│                                      ↓
│                                    RAG → Neo4j KG
│                                      ↓
│                                GNN Predictor
│                                      ↓
│                                Visualization
└── FAIL → Re-run QC (max 3 attempts)

```


Each stage is represented as a **LangGraph node**, which calls either a Nextflow workflow or a Python module.

---

## Directory Layout

```
project_root/
│
├── workflows/                      
│   ├── qc_pipeline.nf
│   ├── exome_pipeline.nf
│   ├── wgs_pipeline.nf
│   ├── rnaseq_pipeline.nf
│   ├── lcms_pipeline.nf
│   ├── scrnaseq_pipeline.nf
│   └── __init__.py
│
├── src/                            
│   ├── __init__.py
│   │
│   ├── pipeline/                   
│   │   ├── langgraph_workflow.py
│   │   └── __init__.py
│   │
│   ├── data/                       
│   │   ├── neo4j_extractor.py
│   │   ├── graph_builder.py
│   │   └── feature_engineering.py
│   │
│   ├── analysis/                   
│   │   ├── pathway_enrichment_gene_list.py
│   │   ├── gene_pathway_network_visualization_rag.py
│   │   └── __init__.py
│   │
│   ├── rag/                        
│   │   ├── rag_pipeline_langchain.py
│   │   ├── structured_drug_kg.py
│   │   └── __init__.py
│   │
│   ├── gnn_predictor/               
│   │   ├── models/
│   │   │   ├── gnn_predictor.py
│   │   │   └── __init__.py
│   │   ├── utils/
│   │   │   ├── loggers.py
│   │   │   └── __init__.py
│   │   ├── visualization/
│   │   │   ├── visualization.py
│   │   │   └── __init__.py
│   │   └── data/
│   │       ├── feature_engineering.py
│   │       ├── graph_builder.py
│   │       ├── neo4j_extractor.py
│   │       └── __init__.py
│   │
│   └── utils/                      
│       ├── helpers.py
│       ├── loggers.py
│       └── __init__.py
│
├── out/                             
│   ├── bamfiles/
│   ├── figures/
│   └── tables/
│
├── logs/                            
├── requirements.txt                 
├── environment.yml                  
└── README.md

````

---

## Example Usage

### 1. Pathway Enrichment
```bash
python src/analysis/pathway_enrichment_gene_list.py
````

* Input: gene list from VCF/annotation
* Output: enrichment JSON & top pathway plots

### 2. Gene → Pathway Network

```bash
python src/analysis/gene_pathway_network_visualization_rag.py
```

* Input: enrichment JSON
* Output: directed graph of genes → pathways

### 3. RAG Pipeline

```bash
python -m src.rag.rag_pipeline_langchain --query "genes linked to Parkinson's disease" --structured
```

* Queries Neo4j knowledge graph
* Output: structured JSON, optional KG population

### 4. GNN Predictor

```python
from src.gnn.gnn_predictor import GNNPredictor
from src.data.graph_builder import build_graph_from_json
from src.data.feature_engineering import generate_node_features

graph = build_graph_from_json("results/enrichment/enrichment_results.json")
features = generate_node_features(graph)

model = GNNPredictor(input_dim=features.shape[1], hidden_dim=64, output_dim=2)
model.train(graph, features)
predictions = model.predict(graph, features)
```

---

## Quick Start

1. **Clone repository and install dependencies**

```bash
git clone https://github.com/man4ish/langgraph_bioflow.git
cd langgraph_bioflow
pip install -r requirements.txt
```

2. **Install Nextflow**

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version
```

3. **Run the orchestrated workflow**

```bash
python main.py
```

> Ensure Neo4j is running for RAG + GNN modules.

---

## Future Enhancements

* Real FastQC parsing and LLM-driven QC evaluation
* Interactive dashboards (FastAPI / Streamlit) for workflow visualization and prediction results
* Multi-omics support: RNA-Seq, WGS, metagenomics
* Advanced GNN architectures for novel drug prediction
* Persistent LangGraph state storage for reproducibility


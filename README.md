# LangGraph BioFlow

### *AI-Orchestrated Bioinformatics Workflow Supervisor (Nextflow + LangGraph + RAG + GNN)*

---

## Overview

**LangGraph BioFlow** is a next-generation **intelligent pipeline orchestration framework** that combines:

* **Nextflow** for scalable bioinformatics workflows
* **LangGraph** for AI-driven control flow and decision-making
* **RAG (Retrieval-Augmented Generation)** for biological knowledge queries
* **Graph Neural Networks (GNNs)** for drug–gene interaction prediction
* **FastAPI + Streamlit** for real-time monitoring and visualization

---

## Key Features

### **1. AI-Supervised Pipeline Execution**

LangGraph manages full workflow logic:

```
QC → Nextflow Pipeline → Pathway Enrichment → RAG → GNN → Visualization
```

Dynamic behaviors include:

* Automatic QC retries (max 3)
* AI-based decisions: continue, stop, or switch pipeline
* Persistent workflow state
* Real-time step-level monitoring

---

### **2. Multi-Omics Pipelines (Nextflow)**

Fully modular workflows:

* WGS
* Exome
* RNA-Seq
* LC-MS
* scRNA-Seq
* QC pipeline

All pipelines run using tested Docker images built for this project.

---

### **3. Intelligent Understanding Layer**

The system integrates:

* **Pathway enrichment**
* **Gene → Pathway network visualization**
* **RAG querying against a biological KG (Neo4j)**
* **GNN predictor for drug ranking**

RAG: Finds relationships like
“genes linked to Parkinson’s disease”, “drug–target interactions”, etc.

GNN: Predicts drug scores using graph embeddings.

---

### **4. Real-Time Dashboard (Streamlit)**

A live dashboard displays:

* Current pipeline status
* Step-by-step logs
* Plots / gene-pathway networks
* RAG structured outputs
* GNN predicted drug candidates
* Final summary & visualizations

The dashboard auto-refreshes during pipeline execution.

---

### **5. REST API (FastAPI)**

A production-ready API provides:

* Run full pipeline: `/run/full`
* Check status: `/status/{run_id}`
* Retrieve outputs: `/results/{run_id}`

Designed for integration with:

* LIMS
* Lab workflow managers
* Hospital EHR systems
* Research analytics dashboards

---

## Project Architecture

![Multi-omics Architecture Diagram](https://raw.githubusercontent.com/man4ish/langgraph_bioflow/dev-nextflow-modular/images/multi-omics_architecture_diagram.png)

### Workflow Logic

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

Each stage is defined as a **LangGraph node**, triggering either:

* a Nextflow pipeline
* a Python analysis module

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
│   ├── pipeline/langgraph_workflow.py
│   ├── analysis/
│   ├── rag/
│   ├── gnn_predictor/
│   └── utils/
│
├── api/                             
│   ├── main.py   ← FastAPI backend
│   └── __init__.py
│
├── dashboard/                       
│   └── app.py   ← Streamlit real-time visualization
│
├── logs/
├── out/
├── requirements.txt
├── environment.yml
└── README.md
```

---

## Running the Full System

### 1. Start FastAPI backend

```bash
uvicorn api.main:app --reload --port 8000
```

### 2. Start dashboard

```bash
streamlit run dashboard/app.py
```

### 3. Trigger full workflow

```bash
curl -X POST http://localhost:8000/run/full
```

Dashboard will auto-update.

---

## Roadmap

* Websocket-based live pipeline streaming
* LIMS-X integration (REST)
* Multi-user workflow launches
* Cloud execution (AWS Batch / Google Life Sciences)
* Multi-omics fusion models
* Real QC scoring using deep learning


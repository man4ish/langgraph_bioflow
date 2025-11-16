# RAG-Powered Gene Discovery Package

## Overview
This package provides a **Retrieval-Augmented Generation (RAG) pipeline** for biomedical knowledge discovery.  
It integrates gene and drug information from pathway enrichment, PubMed literature, and structured databases, producing outputs suitable for visualization, network analysis, or GNN prediction.

---

## Package Structure

```

src/rag/
├── **init**.py
├── rag_pipeline_langchain.py     # Main RAG pipeline
├── data_loader.py                # Utilities to load gene, drug, and literature data
├── embedding_engine.py           # Handles embeddings (Ollama / HuggingFace)
├── structured_drug_kg.py         # Structured Drug-Target-Cancer KG population
├── plot_rag_outputs_network.py   # Generate network plots from RAG output
├── visualize_network.py          # General network visualization utilities

````

---

## Module Descriptions

### 1. `rag_pipeline_langchain.py`
- **Purpose:** Main pipeline for RAG-powered gene discovery.
- **Features:**
  - Retrieves top PubMed abstracts using FAISS via a custom LangChain retriever.
  - Generates answers using Ollama LLM.
  - Supports structured extraction and populates Neo4j knowledge graphs.
- **Usage Example:**

```bash
python -m src.rag.rag_pipeline_langchain \
    --query "genes linked to Parkinson's disease" \
    --structured
````

---

### 2. `data_loader.py`

* **Purpose:** Load and preprocess gene lists, drug-target data, and PubMed abstracts.
* **Functionality:** Provides standardized input for the RAG pipeline.

### 3. `embedding_engine.py`

* **Purpose:** Generate embeddings for queries and documents.
* **Supports:** Ollama embeddings, HuggingFace embeddings.
* **Role:** Essential for retrieval in FAISS-based RAG.

### 4. `structured_drug_kg.py`

* **Purpose:** Populate a structured **Drug-Target-Cancer Knowledge Graph** in Neo4j.
* **Functionality:** Takes structured extraction from RAG and builds nodes and edges in KG.

### 5. `plot_rag_outputs_network.py`

* **Purpose:** Plot the gene-drug-pathway-disease network derived from RAG output.
* **Features:** Directed graph with color-coded node types (gene, pathway, drug, disease).

### 6. `visualize_network.py`

* **Purpose:** General-purpose utilities to visualize networks.
* **Integration:** Can visualize both RAG-derived and pathway enrichment networks.

---

## Workflow

The typical workflow for this package is:

```
           Gene / Drug List
                   │
                   ▼
          Pathway Enrichment
                   │
                   ▼
       ---------------------------
       | RAG Pipeline Module     |
       | (rag_pipeline_langchain)|
       ---------------------------
                   │
       ┌───────────┴───────────┐
       ▼                       ▼
 Structured Output JSON     Neo4j KG Population
 (gene → pathway → disease)    (nodes + edges)
       │
       ▼
   Network Visualization
 (plot_rag_outputs_network.py / visualize_network.py)
       │
       ▼
 GNN Predictions / Literature Summaries
```

---

## Installation

Install required dependencies (Python 3.9+ recommended):

```bash
pip install langchain gseapy networkx matplotlib seaborn faiss-cpu pandas ollama-client
```

---

## Example Usage

```python
from src.rag.rag_pipeline_langchain import RAGAssistant

query = "genes linked to Parkinson's disease"
assistant = RAGAssistant()
summary, pmids, structured_results = assistant.run_pipeline(query, top_k=5, structured=True)

print("Summary:", summary)
print("Top PMIDs:", pmids)
print("Structured Results:", structured_results)
```

**Plot the network:**

```python
from src.rag.plot_rag_outputs_network import plot_network

plot_network(structured_results, save_path="results/rag_network.png")
```

---

## Notes

* FAISS index and Ollama embeddings must be properly configured before running the pipeline.
* Structured outputs are JSON files that can be reused for visualization or GNN analysis.
* Neo4j credentials should be set in `structured_drug_kg.py` or configuration files for KG population.


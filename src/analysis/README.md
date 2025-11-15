# Bioinformatics Analysis Package

## Overview
This package provides tools to perform **pathway enrichment analysis** from a gene list, build a **gene → pathway network**, visualize it, and generate structured queries suitable for a **RAG (Retrieval-Augmented Generation) pipeline**.  

It is designed to be integrated into a bioinformatics workflow where you have gene lists from **VCF annotation**, **RNA-Seq**, or other high-throughput experiments.

---

## Package Structure

```

src/analysis/
├── **init**.py
├── pathway_enrichment_gene_list.py     # Pathway enrichment & visualization
├── gene_pathway_network_visualization_rag.py   # Build network, visualize, generate RAG queries

````

---

## Modules

### 1. `pathway_enrichment_gene_list.py`
- **Purpose:** Perform pathway enrichment analysis on a gene list.
- **Features:**
  - Supports multiple databases (Reactome, KEGG, GO Biological Process, etc.).
  - Returns top enriched pathways.
  - Saves results in JSON format for downstream analysis.
  - Generates bar plots of top enriched pathways.
- **Usage Example:**

```python
from analysis.pathway_enrichment_gene_list import GeneListEnrichment

gene_list = ["TP53","BRCA1","BRCA2","LRRK2","SNCA","APP"]
enrichment = GeneListEnrichment(gene_list)
enrichment.run_enrichment(database='Reactome_2022')
enrichment.plot_top_pathways(10, save_path="results/enrichment/top_pathways.png")
````

---

### 2. `gene_pathway_network_visualization_rag.py`

* **Purpose:** Build a **structured gene → pathway network** from enrichment JSON and generate RAG queries.
* **Features:**

  * Builds a `networkx` directed graph from enrichment JSON.
  * Visualizes the network with genes and pathways as different node types.
  * Generates structured RAG queries like:

    ```
    Genes: TP53, BRCA1, BRCA2 in pathway Cell Cycle
    ```
* **Usage Example:**

```python
from analysis.gene_pathway_network_visualization_rag import (
    build_gene_pathway_network,
    plot_gene_pathway_network,
    generate_rag_queries_from_json
)

json_file = "results/enrichment/enrichment_results.json"
G = build_gene_pathway_network(json_file)
plot_gene_pathway_network(G)

rag_queries = generate_rag_queries_from_json(json_file)
for q in rag_queries[:5]:
    print(q)
```

---

## Workflow

Here is the **high-level workflow** of the package:

```
      Annotated VCF / Gene List
                │
                ▼
  -------------------------------
  | Pathway Enrichment Module   |
  | (pathway_enrichment_gene_list.py) |
  -------------------------------
                │
                ▼
     Enrichment Results JSON
                │
                ▼
  -------------------------------
  | Network & RAG Module        |
  | (gene_pathway_network_visualization_rag.py) |
  -------------------------------
                │
         ┌───────┴────────┐
         ▼                ▼
  Gene → Pathway       Structured RAG Queries
     Network                (for RAG pipeline)
         │
         ▼
 Visualization / GNN / Further Analysis
```

---

## Installation

Install required dependencies (Python 3.9+ recommended):

```bash
pip install gseapy networkx matplotlib seaborn
```

---

## Notes

* **JSON Output:** `results/enrichment/enrichment_results.json` is compatible with the network builder and RAG pipeline.
* **Visualization:** Bar plots and network plots are saved in `results/enrichment/`.
* **Extensibility:** You can easily extend the RAG query generation function to include drugs, diseases, or additional metadata.


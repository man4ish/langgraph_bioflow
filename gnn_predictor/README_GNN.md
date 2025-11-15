
# **GNN Predictor Module**

**Location:** `gnn_predictor/models/gnn_predictor.py`

## **Purpose**

The `gnn_predictor` module implements a Graph Neural Network (GNN) pipeline for predicting drug-target interactions and potential drug repurposing opportunities. It leverages a knowledge graph constructed from structured drug-target-cancer data stored in Neo4j. This module integrates graph-based deep learning with biomedical knowledge to generate predictions and ranking scores for novel associations.

## **Inputs**

* **Nodes:** Drugs, Targets (proteins), and Cancers with attributes such as type, name, and embeddings.
* **Edges:** Relationships like `TARGETS` (Drug→Target) and `USED_FOR` (Drug→Cancer), with additional properties (e.g., pmid, mechanism).
* **Features:** Node and edge features including chemical descriptors, protein embeddings, and graph metrics.

## **Workflow**

1. **Graph Extraction:** Fetch nodes and edges from Neo4j.
2. **Graph Construction:** Build PyTorch Geometric graph objects.
3. **Feature Engineering:** Compute embeddings and graph properties.
4. **GNN Training & Prediction:** Learn node representations and predict potential drug-target/disease links.
5. **Outputs:** Predicted edges, confidence scores, and node embeddings for downstream tasks.

## **Downstream Use**

* Feed predictions into RAG pipelines for literature validation.
* Visualize high-confidence associations in network plots.
* Support drug repurposing and pathway prioritization analyses.

## **Future Extensions**

* Multi-modal embeddings (omics, chemical, clinical).
* Novel GNN architectures for improved predictions.
* Feedback integration from literature-based validation.
* Enhanced visualization for interactive exploration of predictions.


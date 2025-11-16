
# **GNN-Based Drug–Target Node & Edge Significance Predictor**

## **Overview**

This repository contains a Python script that predicts:

1. **Node-level significance** – how important or central a gene, drug, or disease node is within a biomedical knowledge graph.
2. **Edge-level interaction likelihood** – the predicted strength of a drug–target interaction based on graph embeddings learned by a Graph Neural Network (GNN).
```
It is designed to integrate results from **multiple RAG (Retrieval-Augmented Generation) queries**, pathway enrichment outputs, or structured knowledge graphs, merging them into a single unified graph to improve predictive power.

[JSON Files] --> [Graph Construction] --> [Node Feature Engineering] --> [GNN Training] --> [Embeddings]
                                                                                  |
                                                                                  v
            +----------------------+-----------------------+
            |                                              |
   [Node-Level Significance]                     [Edge-Level Interaction Scores]
            |                                              |
   node_significance.csv                        edge_interaction_scores.csv
            |                                              |
            +----------------------+-----------------------+
                                   |
                             [Graph Visualization]
                           gnn_network_edge.png
```
---

## **Scientific Motivation**

* Drug–target interaction (DTI) prediction is central to **drug discovery, repurposing, and personalized medicine**.
* Biomedical data is inherently graph-structured: drugs interact with targets, genes participate in multiple pathways, and diseases share molecular mechanisms.
* Graph Neural Networks (GNNs) allow **learning from both node features and graph topology**, capturing relational dependencies that traditional ML or statistical methods cannot.
* Node-level significance identifies **key biological entities** (genes, drugs, pathways) with high network centrality or influence.
* Edge-level predictions provide **likelihood of interaction between drugs and targets**, enabling prioritization for experimental validation or drug repurposing.

---

## **Technical Description**

The script performs the following workflow:

1. **Input Parsing**

   * Accepts one or more JSON files containing drug–target, drug–disease, or drug–phenotype relationships.
   * Each entry may include metadata such as `pmid` for literature evidence.

2. **Graph Construction**

   * Constructs a **NetworkX graph**:

     * Nodes are labeled as `"drug"` or `"target"`.
     * Edges represent known relationships, annotated with evidence.

3. **Feature Engineering**

   * Node features are computed via `generate_node_features()`.
   * Can be extended to include pathway scores, embeddings, or network centrality measures.

4. **GNN Training**

   * A `GNNPredictor` learns embeddings for each node using the graph structure.
   * Message passing aggregates neighborhood information to compute latent representations.

5. **Node Significance**

   * Node embeddings are used to compute a **significance score**, reflecting:

     * centrality
     * connectivity
     * relevance within the graph
     * potential biological influence

6. **Edge-Level Interaction Prediction**

   * Edge scores are computed from node embeddings (default: dot product).
   * High scores indicate a higher likelihood of **drug-target interaction**.

7. **Outputs**

   * `node_significance.csv`: Node-level importance scores.
   * `edge_interaction_scores.csv`: Predicted likelihood for each drug–target pair.
   * `gnn_network_edge.png`: Graph visualization with optional edge/ node annotations.

---

## **Usage**

### **Command-Line Arguments**

```bash
--json-files   List of one or more JSON files containing drug-target relationships.
--output-dir   Directory to save outputs (default: results).
--node-level   Flag to compute node-level significance.
--edge-level   Flag to compute edge-level interaction scores.
```

### **Example Runs**

1. **Node significance only**

```bash
python test/predict_gnn_significance.py \
    --json-files data/query1.json data/query2.json \
    --node-level
```

2. **Edge interaction scores only**

```bash
python test/predict_gnn_significance.py \
    --json-files data/query1.json data/query2.json \
    --edge-level
```

3. **Both node + edge outputs**

```bash
python test/predict_gnn_significance.py \
    --json-files data/*.json \
    --node-level --edge-level
```

---

## **Scientific Interpretation of Outputs**

### **Node Significance**

* High node scores indicate that a node (drug, gene, or disease) is **highly central and influential** in the graph.
* Example: TCF7L2 may have a high score due to multiple disease associations and connectivity to key pathways.
* Useful for:

  * Identifying key drug targets
  * Prioritizing central genes for experimental follow-up
  * Mapping influential nodes in multi-omics networks

### **Edge Interaction Scores**

* High edge scores indicate a **strong likelihood of drug–target interaction**, accounting for graph structure and node embeddings.
* Useful for:

  * Drug repurposing
  * Prioritizing drug–gene pairs for experimental validation
  * Understanding network propagation effects in biomedical knowledge graphs

### **Combined Interpretation**

* Node scores highlight **important entities**.
* Edge scores highlight **specific actionable interactions**.
* Together they provide **multi-level insight** into biomedical networks.

---

## **File Structure**

```
test/
 ├── predict_gnn_significance.py   # Main script
 ├── README_GNN_SIGNIFICANCE.md    # This README
 └── data/                         # JSON input files for testing
```

---

## **Requirements**

* Python 3.9+
* NetworkX
* NumPy
* PyTorch / PyTorch Geometric (required by GNNPredictor)
* Matplotlib / Seaborn (for visualization)
* Pandas (optional for CSV handling)

```bash
pip install torch torch_geometric networkx matplotlib seaborn pandas
```

---

## **Future Enhancements**

* Multi-relational edge types (e.g., different evidence types, pathway links)
* Edge scoring via MLP on concatenated node embeddings or cosine similarity
* Integration with Neo4j knowledge graphs for large-scale biomedical networks
* GPU acceleration for large graphs
* Interactive visualization dashboards using Plotly or Bokeh



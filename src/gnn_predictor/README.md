# GNN Predictor Package

## Overview
The `gnn_predictor` package provides tools for building, training, and visualizing **graph neural network (GNN) models** applied to biomedical networks such as gene → pathway → disease networks.  
It integrates with structured outputs from **pathway enrichment** or **RAG pipelines**, allowing predictions over biological graphs.

---

## Package Structure

```

gnn_predictor/
├── models/
│   ├── gnn_predictor.py        # Main GNN model class and training routines
│   └── **init**.py
├── utils/
│   ├── loggers.py              # Logging utilities
│   └── **init**.py
├── visualization/
│   ├── visualization.py        # Network and prediction visualization
│   └── **init**.py
├── data/
│   ├── feature_engineering.py  # Feature extraction for graph nodes/edges
│   ├── graph_builder.py        # Build graph from Neo4j or JSON
│   ├── neo4j_extractor.py      # Extract structured graph data from Neo4j
│   └── **init**.py

````

---

## Modules

### 1. `models/gnn_predictor.py`
- **Purpose:** Define and train GNN models on graph data.
- **Features:**
  - Supports node-level and edge-level prediction tasks.
  - Integrates with `torch_geometric` or PyG-style graphs.
- **Usage Example:**
```python
from models.gnn_predictor import GNNPredictor

model = GNNPredictor(input_dim=128, hidden_dim=64, output_dim=2)
model.train(graph_data)  # graph_data from graph_builder or Neo4j extractor
preds = model.predict(graph_data)
````

### 2. `data/graph_builder.py`

* **Purpose:** Build a network graph from enrichment results or RAG outputs.
* **Outputs:** `networkx` or PyG-compatible graph objects.

### 3. `data/feature_engineering.py`

* **Purpose:** Generate features for nodes and edges in the graph.
* **Examples:** gene expression values, pathway scores, node centrality, embedding vectors.

### 4. `data/neo4j_extractor.py`

* **Purpose:** Extract structured gene-drug-pathway networks from a Neo4j KG.

### 5. `visualization/visualization.py`

* **Purpose:** Visualize graphs, predictions, and embeddings.
* **Example Features:** node coloring by type, edge labels, prediction heatmaps.

### 6. `utils/loggers.py`

* **Purpose:** Standardized logging for training, evaluation, and experiments.

---

## Workflow

```
   Pathway / RAG Output JSON or Neo4j KG
                  │
          --------------------
          | Graph Builder     |
          --------------------
                  │
          --------------------
          | Feature Engineering |
          --------------------
                  │
          --------------------
          | GNN Predictor      |
          --------------------
                  │
          --------------------
          | Visualization      |
          --------------------
                  │
         Predictions / Insights
```

* **Input:** Enrichment JSON, RAG output, or Neo4j graph.
* **Processing:** Graph construction → Feature engineering → GNN training/prediction → Visualization.
* **Output:** Node/edge predictions, graphs, plots.

---

## Installation

Install required dependencies (Python 3.9+ recommended):

```bash
pip install torch torch_geometric networkx matplotlib seaborn pandas neo4j
```

---

## Example Usage

```python
from data.graph_builder import build_graph_from_json
from data.feature_engineering import generate_node_features
from models.gnn_predictor import GNNPredictor
from visualization.visualization import plot_graph

# Load and build graph
graph = build_graph_from_json("results/enrichment/enrichment_results.json")

# Generate node/edge features
features = generate_node_features(graph)

# Train GNN
model = GNNPredictor(input_dim=features.shape[1], hidden_dim=64, output_dim=2)
model.train(graph, features)

# Predict
predictions = model.predict(graph, features)

# Visualize predictions
plot_graph(graph, node_predictions=predictions, save_path="results/gnn_network.png")
```

---

## Notes

* The package is designed to integrate **directly with RAG or pathway enrichment outputs**.
* Neo4j KG extraction is optional but recommended for large, structured networks.
* Visualization supports both static plots and interactive network exploration.




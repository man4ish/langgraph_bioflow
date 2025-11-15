import numpy as np
import torch
import pandas as pd
from torch_geometric.data import Data

from gnn_predictor.data.feature_engineering import FeatureEngineer
from gnn_predictor.models.gnn_predictor import GNNEstimator
from gnn_predictor.visualization.visualization import NetworkVisualizer

# -------------------------------------------------------------------

# Synthetic Graph Generator (100 nodes, cluster-based labels)

# -------------------------------------------------------------------

def generate_synthetic_graph(num_nodes=100, num_clusters=3, cluster_sizes=[25, 25, 20]):
rng = np.random.default_rng(42)

```
node_names = [f"Node_{i}" for i in range(num_nodes)]
node_types = rng.choice(["Gene", "Protein", "Drug", "Disease", "Pathway"], size=num_nodes)

labels = np.zeros(num_nodes, dtype=int)
start = 0
for size in cluster_sizes:
    labels[start:start+size] = 1
    start += size

relations = ["activates", "inhibits", "binds", "associated", "interacts"]

edge_list, edge_rel, edge_mech_length, edge_pmid = [], [], [], []

# Dense cluster edges
start = 0
for size in cluster_sizes:
    for i in range(start, start + size):
        for j in range(i + 1, start + size):
            for direction in [(i, j), (j, i)]:
                edge_list.append(direction)
                edge_rel.append(rng.choice(relations))
                edge_mech_length.append(rng.integers(3, 20))
                edge_pmid.append(rng.integers(100000, 999999))
    start += size

# Sparse inter-cluster edges
for _ in range(150):
    i, j = rng.integers(0, num_nodes, size=2)
    if i != j:
        edge_list.append((i, j))
        edge_rel.append(rng.choice(relations))
        edge_mech_length.append(rng.integers(3, 20))
        edge_pmid.append(rng.integers(100000, 999999))

edge_index = torch.tensor(edge_list, dtype=torch.long).t().contiguous()

data = Data(
    x=torch.zeros((num_nodes, 3)),  # placeholder, updated by FeatureEngineer
    edge_index=edge_index,
    node_names=node_names,
    node_types=node_types,
    edge_relations=edge_rel,
    edge_mechanism=edge_mech_length,
    edge_pmid=edge_pmid,
    y=torch.tensor(labels, dtype=torch.float),
)
return data
```

# -------------------------------------------------------------------

# Create synthetic graph

# -------------------------------------------------------------------

print("Generating synthetic data...")
data = generate_synthetic_graph()

# -------------------------------------------------------------------

# Build edges_df for FeatureEngineer

# -------------------------------------------------------------------

edges_df = pd.DataFrame({
"source": data.edge_index[0].cpu().numpy(),
"target": data.edge_index[1].cpu().numpy(),
"relation": data.edge_relations,
"mechanism_length": data.edge_mechanism,
"pmid": data.edge_pmid
})

# -------------------------------------------------------------------

# Feature Engineering

# -------------------------------------------------------------------

print("Building features...")
fe = FeatureEngineer(data, edges_df)
node_features, edge_features = fe.build_features()
data.x = node_features
data.edge_attr = edge_features

# -------------------------------------------------------------------

# Train GNN

# -------------------------------------------------------------------

print("Training model...")
gnn = GNNEstimator(data, hidden_dim=64, lr=0.01, epochs=50)
gnn.train()

# -------------------------------------------------------------------

# Evaluate

# -------------------------------------------------------------------

print("Evaluating...")
auc = gnn.evaluate()
preds = gnn.predict()
print("Predictions:", preds)

# -------------------------------------------------------------------

# Visualize Network

# -------------------------------------------------------------------

print("Saving visualization...")
viz = NetworkVisualizer(data)
viz.build()
viz.show("synthetic_network.html")

print("Done. Output saved as synthetic_network.html")


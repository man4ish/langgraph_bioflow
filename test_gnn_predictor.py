# test_gnn_predictor.py
import pandas as pd
import torch
from torch_geometric.data import Data
from gnn_predictor.data.graph_builder import GraphBuilder
from gnn_predictor.data.feature_engineering import FeatureEngineer
from gnn_predictor.models.gnn_predictor import GNNEstimator
from gnn_predictor.visualization.visualization import NetworkVisualizer

# -------------------------------
# 1. Create dummy nodes and edges
# -------------------------------
nodes_df = pd.DataFrame({
    'node_id': [0, 1, 2, 3],
    'name': ['DrugA', 'Target1', 'CancerX', 'DrugB'],
    'type': ['Drug', 'Target', 'Cancer', 'Drug']
})

edges_df = pd.DataFrame({
    'src': [0, 0, 3],
    'dst': [1, 2, 2],
    'relation': ['TARGETS', 'USED_FOR', 'USED_FOR'],
    'pmid': [12345, 67890, 13579],
    'mechanism': ['mechanism1', 'mechanism2', 'mechanism3']
})

# -------------------------------
# 2. Build PyG graph
# -------------------------------
builder = GraphBuilder(nodes_df, edges_df)
pyg_data = builder.build()

# Add dummy labels for training (binary edges)
pyg_data.y = torch.tensor([1, 0, 1], dtype=torch.float)

print("PyG Data:", pyg_data)
print("Node feature shape:", pyg_data.x.shape)
print("Edge index shape:", pyg_data.edge_index.shape)

# -------------------------------
# 3. Generate features (optional)
# -------------------------------
fe = FeatureEngineer(nodes_df, edges_df)
node_features, edge_features = fe.build_features()

print("Node features:", node_features)
print("Edge features:", edge_features)

# -------------------------------
# 4. Train GNN
# -------------------------------
gnn = GNNEstimator(pyg_data, hidden_dim=16, epochs=20)
gnn.train()
gnn.evaluate()

predictions = gnn.predict()
print("Predictions:", predictions)

# -------------------------------
# 5. Visualization
# -------------------------------
pred_df = pd.DataFrame({
    'drug': ['DrugA', 'DrugB'],
    'target': ['Target1', 'CancerX'],
    'cancer': ['CancerX', 'CancerX'],
    'score': [0.8, 0.9],
    'rag_support': [5, 10]
})

viz = NetworkVisualizer()
viz.build_network(pred_df)
viz.show("test_network.html")

print("Interactive network saved as test_network.html")


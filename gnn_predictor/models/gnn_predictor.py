"""
gnn_predictor.py
----------------
Module to train a GNN on the Drug–Target–Cancer knowledge graph
and predict novel drug-target-cancer associations.

Dependencies:
    - torch
    - torch_geometric
    - pandas
    - sklearn (for evaluation)
"""

import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv, SAGEConv
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
import pandas as pd

class GNNModel(torch.nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(GNNModel, self).__init__()
        self.conv1 = SAGEConv(input_dim, hidden_dim)
        self.conv2 = SAGEConv(hidden_dim, output_dim)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index)
        return x

class GNNEstimator:
    def __init__(self, data: Data, hidden_dim: int = 64, lr: float = 0.01, epochs: int = 100):
        """
        data: PyTorch Geometric Data object containing nodes/features/edges
        """
        self.data = data
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model = GNNModel(data.num_node_features, hidden_dim, 1).to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=lr)
        self.epochs = epochs

    def train(self):
        self.model.train()
        for epoch in range(1, self.epochs + 1):
            self.optimizer.zero_grad()
            out = self.model(self.data.x.to(self.device), self.data.edge_index.to(self.device))
            # Assuming binary labels for edges stored in data.y
            loss = F.binary_cross_entropy_with_logits(out.squeeze(), self.data.y.to(self.device).float())
            loss.backward()
            self.optimizer.step()
            if epoch % 10 == 0:
                print(f"Epoch {epoch}/{self.epochs}, Loss: {loss.item():.4f}")

    @torch.no_grad()
    def predict(self):
        self.model.eval()
        out = self.model(self.data.x.to(self.device), self.data.edge_index.to(self.device))
        probs = torch.sigmoid(out).squeeze().cpu().numpy()
        return probs

    @torch.no_grad()
    def evaluate(self):
        preds = self.predict()
        auc = roc_auc_score(self.data.y.cpu().numpy(), preds)
        print(f"ROC-AUC Score: {auc:.4f}")
        return auc

# ---------------------------
# Example usage
# ---------------------------
if __name__ == "__main__":
    from data.graph_builder import GraphBuilder
    from data.neo4j_extractor import Neo4jExtractor
    from data.feature_engineering import FeatureEngineer

    # 1. Extract KG from Neo4j
    extractor = Neo4jExtractor()
    nodes_df, edges_df = extractor.load_graph()
    extractor.close()

    # 2. Build PyG graph
    gb = GraphBuilder(nodes_df, edges_df)
    data = gb.build_graph()

    # 3. Add features
    fe = FeatureEngineer(data, nodes_df, edges_df)
    data = fe.add_node_features()

    # 4. Train GNN
    gnn = GNNEstimator(data, hidden_dim=64, epochs=50)
    gnn.train()
    gnn.evaluate()
    predictions = gnn.predict()
    print("Predictions:", predictions[:10])

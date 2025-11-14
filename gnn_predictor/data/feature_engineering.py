"""
feature_engineering.py
----------------------
Generate node and edge features for the Drug-Target-Cancer KG.

Inputs:
    - nodes_df: DataFrame with node_id, name, type
    - edges_df: DataFrame with src, dst, relation, mechanism, pmid

Outputs:
    - node_features: torch.Tensor (num_nodes x feature_dim)
    - edge_features: torch.Tensor (num_edges x feature_dim)
"""

import torch
import pandas as pd
from sklearn.preprocessing import OneHotEncoder

class FeatureEngineer:
    def __init__(self, nodes_df: pd.DataFrame, edges_df: pd.DataFrame):
        self.nodes_df = nodes_df
        self.edges_df = edges_df

    # -------------------------------------------------------
    # Node features
    # -------------------------------------------------------
    def build_node_features(self) -> torch.Tensor:
        """
        One-hot encode node types
        """
        types = self.nodes_df['type'].values.reshape(-1, 1)
        encoder = OneHotEncoder(sparse=False)
        node_features = encoder.fit_transform(types)
        return torch.tensor(node_features, dtype=torch.float)

    # -------------------------------------------------------
    # Edge features
    # -------------------------------------------------------
    def build_edge_features(self) -> torch.Tensor:
        """
        Encode relation type (one-hot) and mechanism length as edge features
        """
        # Relation one-hot
        rels = self.edges_df['relation'].values.reshape(-1, 1)
        rel_encoder = OneHotEncoder(sparse=False)
        rel_features = rel_encoder.fit_transform(rels)

        # Mechanism length as numerical feature
        mechanism_length = self.edges_df['mechanism'].fillna("").apply(len).values.reshape(-1, 1)

        # Combine
        edge_features = torch.tensor(
            pd.np.hstack([rel_features, mechanism_length]),
            dtype=torch.float
        )
        return edge_features

    # -------------------------------------------------------
    # Master wrapper
    # -------------------------------------------------------
    def build_features(self) -> (torch.Tensor, torch.Tensor):
        node_features = self.build_node_features()
        edge_features = self.build_edge_features()
        return node_features, edge_features


# -------------------------------------------------------
# Standalone test
# -------------------------------------------------------
if __name__ == "__main__":
    from neo4j_extractor import Neo4jExtractor

    extractor = Neo4jExtractor()
    nodes_df, edges_df = extractor.load_graph()
    extractor.close()

    fe = FeatureEngineer(nodes_df, edges_df)
    node_features, edge_features = fe.build_features()

    print("Node features shape:", node_features.shape)
    print("Edge features shape:", edge_features.shape)

"""
graph_builder.py
----------------
Builds a PyTorch Geometric graph from Neo4j Drug-Target-Cancer KG.

Inputs:
    - nodes_df: DataFrame from Neo4jExtractor (node_id, name, type)
    - edges_df: DataFrame from Neo4jExtractor (src, dst, relation, pmid, mechanism)

Outputs:
    - PyG Data object with node features and edge index
"""

import torch
from torch_geometric.data import Data
import pandas as pd
from typing import Tuple

class GraphBuilder:
    def __init__(self, nodes_df: pd.DataFrame, edges_df: pd.DataFrame):
        self.nodes_df = nodes_df
        self.edges_df = edges_df
        self.node_id_map = {nid: i for i, nid in enumerate(nodes_df['node_id'].tolist())}

    def build_edge_index(self) -> torch.Tensor:
        """
        Convert edges to PyG edge_index tensor (2, num_edges)
        """
        src = [self.node_id_map[s] for s in self.edges_df['src']]
        dst = [self.node_id_map[d] for d in self.edges_df['dst']]
        edge_index = torch.tensor([src, dst], dtype=torch.long)
        return edge_index

    def build_node_features(self) -> torch.Tensor:
        """
        Simple feature: one-hot encoding of node type (Drug, Target, Cancer)
        """
        node_types = sorted(self.nodes_df['type'].unique())
        type_to_idx = {t: i for i, t in enumerate(node_types)}
        features = torch.zeros((len(self.nodes_df), len(node_types)), dtype=torch.float)

        for i, t in enumerate(self.nodes_df['type']):
            features[i, type_to_idx[t]] = 1.0

        return features

    def build(self) -> Data:
        """
        Build PyTorch Geometric Data object
        """
        edge_index = self.build_edge_index()
        x = self.build_node_features()

        # Optionally, you can store node names and edge attributes for reference
        data = Data(x=x, edge_index=edge_index)

        # Store metadata
        data.node_names = self.nodes_df['name'].tolist()
        data.node_types = self.nodes_df['type'].tolist()
        data.edge_relations = self.edges_df['relation'].tolist()
        data.edge_mechanism = self.edges_df.get('mechanism', []).tolist()
        data.edge_pmid = self.edges_df.get('pmid', []).tolist()

        return data

# -------------------------------------------------------
# Standalone test
# -------------------------------------------------------
if __name__ == "__main__":
    from neo4j_extractor import Neo4jExtractor

    extractor = Neo4jExtractor()
    nodes_df, edges_df = extractor.load_graph()
    extractor.close()

    builder = GraphBuilder(nodes_df, edges_df)
    pyg_data = builder.build()

    print("\n=== PyG Data ===")
    print(pyg_data)
    print("Node feature shape:", pyg_data.x.shape)
    print("Edge index shape:", pyg_data.edge_index.shape)

"""
GNN-Based Drug–Target Node & Edge Significance Predictor

This script allows computing:

1. Node-level significance: importance of each node in the graph.
2. Edge-level likelihood: predicted strength of each drug–target pair.

Supports multiple JSON input files for merging into one unified graph.
"""

import argparse
import json
import networkx as nx
import numpy as np
import os

from models.gnn_predictor import GNNPredictor
from data.feature_engineering import generate_node_features
from visualization.visualization import plot_graph


def build_graph_from_json(json_path):
    """Builds a drug–target graph from a single JSON file."""
    with open(json_path, "r") as f:
        data = json.load(f)

    G = nx.Graph()
    for entry in data:
        drug = entry.get("drug") or entry.get("drug_target")
        target = (
            entry.get("target")
            or entry.get("disease")
            or entry.get("condition")
            or entry.get("relevant_phenotype")
        )

        if drug and target:
            G.add_node(drug, type="drug")
            G.add_node(target, type="target")
            G.add_edge(drug, target, pmid=entry.get("pmid"))
    return G


def build_graph_from_json_files(json_paths):
    """Merges multiple JSON-derived graphs into a single unified graph."""
    G = nx.Graph()
    for path in json_paths:
        print(f"Loading JSON: {path}")
        Gi = build_graph_from_json(path)
        G = nx.compose(G, Gi)
    return G


def compute_edge_scores(G, embeddings, node_list):
    """Compute edge-level scores using dot product of node embeddings."""
    node_idx = {n: i for i, n in enumerate(node_list)}
    edge_scores = []
    for u, v in G.edges():
        emb_u = embeddings[node_idx[u]]
        emb_v = embeddings[node_idx[v]]
        score = float(np.dot(emb_u, emb_v))  # simple edge score
        edge_scores.append((u, v, score))
    return edge_scores


def main():
    parser = argparse.ArgumentParser(description="Compute GNN node/edge significance")
    parser.add_argument(
        "--json-files",
        nargs="+",
        required=True,
        help="List of JSON files with drug-target data.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results",
        help="Directory to save predictions and visualization.",
    )
    parser.add_argument(
        "--node-level",
        action="store_true",
        help="Compute node-level significance",
    )
    parser.add_argument(
        "--edge-level",
        action="store_true",
        help="Compute edge-level interaction scores",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Build merged graph
    G = build_graph_from_json_files(args.json_files)
    print(f"Unified graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Generate node features
    X, node_list = generate_node_features(G)

    # Train GNN
    input_dim = X.shape[1]
    model = GNNPredictor(input_dim=input_dim, hidden_dim=64, output_dim=1)
    model.train(G, X)

    # Get node embeddings
    embeddings = model.get_node_embeddings(G, X)  # assume this returns np.array of shape (num_nodes, embedding_dim)

    # Node-level significance
    if args.node_level:
        node_scores = embeddings.mean(axis=1)  # simple: mean embedding value as significance
        node_out = os.path.join(args.output_dir, "node_significance.csv")
        with open(node_out, "w") as f:
            f.write("node,score\n")
            for n, s in zip(node_list, node_scores):
                f.write(f"{n},{s:.4f}\n")
        print(f"Node-level significance saved to {node_out}")

    # Edge-level interaction
    if args.edge_level:
        edge_scores = compute_edge_scores(G, embeddings, node_list)
        edge_out = os.path.join(args.output_dir, "edge_interaction_scores.csv")
        with open(edge_out, "w") as f:
            f.write("drug,target,score\n")
            for u, v, s in edge_scores:
                f.write(f"{u},{v},{s:.4f}\n")
        print(f"Edge-level scores saved to {edge_out}")

        # Optional visualization
        out_png = os.path.join(args.output_dir, "gnn_network_edge.png")
        plot_graph(G, node_predictions=None, edge_scores=edge_scores, save_path=out_png)
        print(f"Graph visualization saved to {out_png}")

    print("Done.")


if __name__ == "__main__":
    main()

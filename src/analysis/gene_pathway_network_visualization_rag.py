"""
Module: Gene → Pathway Network Visualization

Author: Manish Kumar
Date: 2025-11-15
Description: Builds a structured gene → pathway network from pathway 
             enrichment JSON results, visualizes it, and generates RAG queries.
"""

import json
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os

def build_gene_pathway_network(json_file):
    """
    Build a structured Gene → Pathway network from enrichment JSON.

    Args:
        json_file (str): Path to enrichment JSON file.

    Returns:
        G (networkx.DiGraph): Directed graph with genes and pathways as nodes.
    """
    # Load enrichment results
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    G = nx.DiGraph()  # Directed graph

    for entry in data:
        pathway = entry.get("Term")
        genes = entry.get("Genes", "")
        if not genes or not pathway:
            continue

        gene_list = [g.strip() for g in genes.split(";")]

        # Add pathway node
        G.add_node(pathway, type="pathway")

        # Add gene nodes and edges
        for gene in gene_list:
            G.add_node(gene, type="gene")
            G.add_edge(gene, pathway, relation="part_of")
    
    return G

def plot_gene_pathway_network(G, save_path="results/enrichment/gene_pathway_network.png"):
    """
    Plot the Gene → Pathway network.

    Args:
        G (networkx.DiGraph): NetworkX graph.
        save_path (str): Path to save the figure.
    """
    # Assign node colors based on type
    color_map = {"gene": "skyblue", "pathway": "lightgreen"}
    node_colors = [color_map.get(G.nodes[n].get("type", "other"), "grey") for n in G.nodes()]

    # Layout
    pos = nx.spring_layout(G, k=0.5, seed=42)

    # Draw graph
    plt.figure(figsize=(20, 15))
    nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=1500,
            font_size=10, arrows=True)
    
    # Draw edge labels
    edge_labels = nx.get_edge_attributes(G, "relation")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color="red", font_size=8)

    # Add legend
    legend_elements = [
        mpatches.Patch(color='skyblue', label='Gene'),
        mpatches.Patch(color='lightgreen', label='Pathway')
    ]
    plt.legend(handles=legend_elements, loc='best')

    # Save and show
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Network plot saved to: {save_path}")
    plt.show()

def generate_rag_queries_from_json(json_file):
    """
    Generate structured queries from pathway enrichment JSON
    suitable for passing to a RAG pipeline.

    Args:
        json_file (str): Path to enrichment JSON.

    Returns:
        List[str]: List of queries in the format:
                   "Genes: gene1, gene2, ... in pathway pathway_name"
    """
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    queries = []
    for entry in data:
        pathway = entry.get("Term")
        genes = entry.get("Genes", "")
        if not genes or not pathway:
            continue
        gene_list = [g.strip() for g in genes.split(";")]
        query = f"Genes: {', '.join(gene_list)} in pathway {pathway}"
        queries.append(query)
    
    return queries

# -------------------------
# Example usage
# -------------------------
if __name__ == "__main__":
    json_file = "results/enrichment/enrichment_results.json"
    
    # Build network
    G = build_gene_pathway_network(json_file)
    print(f"Network nodes: {len(G.nodes())}, edges: {len(G.edges())}")
    
    # Plot network
    plot_gene_pathway_network(G)

    # Generate RAG queries
    rag_queries = generate_rag_queries_from_json(json_file)
    print("Sample RAG queries:")
    for q in rag_queries[:5]:
        print("-", q)

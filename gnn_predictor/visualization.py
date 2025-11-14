# visualization.py
import networkx as nx
from pyvis.network import Network
import pandas as pd

class NetworkVisualizer:
    """
    Visualize predicted drug-target-cancer network
    """

    def __init__(self, notebook: bool = False):
        self.net = Network(notebook=notebook, height="750px", width="100%")

    def build_network(self, predictions: pd.DataFrame):
        """
        Add nodes and edges from predictions.
        Expects columns: drug, target, cancer, score, rag_support
        """
        for row in predictions.itertuples():
            # Add nodes with type
            self.net.add_node(row.drug, label=row.drug, color="lightblue", title="Drug")
            self.net.add_node(row.target, label=row.target, color="lightgreen", title="Target")
            self.net.add_node(row.cancer, label=row.cancer, color="lightpink", title="Cancer")

            # Add edges
            self.net.add_edge(row.drug, row.target, value=row.score, title=row.rag_support)
            self.net.add_edge(row.drug, row.cancer, value=row.score, title=row.rag_support)

    def show(self, output_file: str = "gnn_network.html"):
        """
        Render interactive network
        """
        self.net.show(output_file)

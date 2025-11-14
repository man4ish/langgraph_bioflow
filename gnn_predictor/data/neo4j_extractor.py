"""
neo4j_extractor.py
-------------------
Extracts nodes and edges from the Neo4j Drug–Target–Cancer knowledge graph.

Output:
    - nodes_df: DataFrame of all nodes with IDs and types
    - edges_df: DataFrame of all edges with source, target, and relationship type

Used in GNN pipeline to build PyTorch Geometric graph.
"""

from neo4j import GraphDatabase
import pandas as pd
from typing import Tuple


class Neo4jExtractor:
    def __init__(self,
                 uri: str = "bolt://localhost:7687",
                 user: str = "neo4j",
                 password: str = "MyNewSecurePassword123"):
        self.driver = GraphDatabase.driver(uri, auth=(user, password))

    # -------------------------------------------------------
    # Internal helpers
    # -------------------------------------------------------
    @staticmethod
    def _records_to_df(records) -> pd.DataFrame:
        """Converts neo4j result records to pandas DataFrame."""
        if not records:
            return pd.DataFrame()

        return pd.DataFrame([r.data() for r in records])

    # -------------------------------------------------------
    # Node extraction
    # -------------------------------------------------------
    def fetch_nodes(self) -> pd.DataFrame:
        """
        Fetch all nodes: Drug, Target, Cancer.
        """
        query = """
        MATCH (n)
        RETURN id(n) AS node_id, labels(n) AS labels, n.name AS name
        """

        with self.driver.session() as session:
            records = session.run(query).data()

        df = pd.DataFrame(records)

        # Simplify node type (Drug / Target / Cancer)
        df["type"] = df["labels"].apply(lambda x: x[0] if isinstance(x, list) else x)
        df = df.drop(columns=["labels"])

        return df

    # -------------------------------------------------------
    # Edge extraction
    # -------------------------------------------------------
    def fetch_edges(self) -> pd.DataFrame:
        """
        Fetch all relationships from the KG:
            - (Drug)-[:TARGETS]->(Target)
            - (Drug)-[:USED_FOR]->(Cancer)
        """
        query = """
        MATCH (a)-[r]->(b)
        RETURN id(a) AS src,
               id(b) AS dst,
               type(r) AS relation,
               r.pmid AS pmid,
               r.mechanism AS mechanism
        """

        with self.driver.session() as session:
            records = session.run(query).data()

        df = pd.DataFrame(records)

        return df

    # -------------------------------------------------------
    # Master extraction wrapper
    # -------------------------------------------------------
    def load_graph(self) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Return nodes_df, edges_df.
        """
        nodes_df = self.fetch_nodes()
        edges_df = self.fetch_edges()

        return nodes_df, edges_df

    # -------------------------------------------------------
    # Cleanup
    # -------------------------------------------------------
    def close(self):
        self.driver.close()


# -------------------------------------------------------
# Standalone test
# -------------------------------------------------------
if __name__ == "__main__":
    extractor = Neo4jExtractor()

    nodes_df, edges_df = extractor.load_graph()

    print("\n=== Nodes ===")
    print(nodes_df.head())

    print("\n=== Edges ===")
    print(edges_df.head())

    extractor.close()

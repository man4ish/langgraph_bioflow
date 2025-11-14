# pathway_enrichment_gene_list.py
"""
Module: Pathway Enrichment and Visualization from Gene List

Author: Manish Kumar
Date: 2025-11-14
Description: Performs pathway enrichment analysis from a gene list
             and visualizes top enriched pathways.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gseapy as gp  # Python library for gene set enrichment
import numpy as np

class GeneListEnrichment:
    def __init__(self, gene_list, organism='Human', outdir='results/enrichment'):
        """
        gene_list: List of gene symbols
        organism: 'Human' or other supported by gseapy
        outdir: directory to save results
        """
        self.gene_list = gene_list
        self.organism = organism
        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        self.enrich_res = None
    
    def run_enrichment(self, database='Reactome_2022', min_size=5, max_size=5000, fdr_threshold=0.05):
        """
        Perform enrichment using gseapy.
        database: pathway database (Reactome, KEGG, GO_Biological_Process)
        """
        print(f"Running enrichment on {len(self.gene_list)} genes using {database}")
        enr = gp.enrichr(
            gene_list=self.gene_list,
            description='gene_list',
            gene_sets=database,
            organism=self.organism,
            outdir=self.outdir,
            cutoff=fdr_threshold
        )
        self.enrich_res = enr.results
        return self.enrich_res
    
    def top_pathways(self, n=10):
        """
        Return top n enriched pathways by adjusted p-value
        """
        if self.enrich_res is None:
            raise ValueError("Run enrichment first")
        return self.enrich_res.sort_values("Adjusted P-value").head(n)
    
    def plot_top_pathways(self, n=10, save_path=None):
        """
        Visualize top n pathways as a horizontal bar plot
        """
        top_df = self.top_pathways(n)
        plt.figure(figsize=(10, max(5, n*0.5)))
        sns.barplot(
            x=-np.log10(top_df['Adjusted P-value']),
            y=top_df['Term'],
            palette='viridis'
        )
        plt.xlabel("-log10(Adjusted P-value)")
        plt.ylabel("Pathway")
        plt.title(f"Top {n} Enriched Pathways")
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300)
            print(f"Saved plot to {save_path}")
        else:
            plt.show()
    
    def get_pathway_genes(self, pathway_name):
        """
        Return genes associated with a specific enriched pathway
        """
        if self.enrich_res is None:
            raise ValueError("Run enrichment first")
        subset = self.enrich_res[self.enrich_res['Term'] == pathway_name]
        if subset.empty:
            return []
        genes = subset.iloc[0]['Genes']
        return [g.strip() for g in genes.split(";")]

# -------------------------
# Example usage
# -------------------------
if __name__ == "__main__":
    # Example gene list (from SnpEff)
    gene_list = ["TP53", "BRCA1", "BRCA2", "LRRK2", "SNCA", "APP"]
    
    enrichment = GeneListEnrichment(gene_list)
    enrichment.run_enrichment(database='Reactome_2022')
    
    # Print top pathways
    print(enrichment.top_pathways(5))
    
    # Plot top pathways
    enrichment.plot_top_pathways(10, save_path="results/enrichment/top_pathways.png")
    
    # Get genes for a specific pathway
    genes = enrichment.get_pathway_genes("Parkinsons disease")
    print("Genes in Parkinsons disease pathway:", genes)

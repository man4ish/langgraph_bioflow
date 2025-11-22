
#!/usr/bin/python
"""
Module: Gene List Pathway Enrichment and Visualization

Author: Manish Kumar
Date: 2025-11-14

Description:
    This module provides a class-based workflow for performing pathway
    enrichment analysis using a list of gene symbols. It leverages the
    gseapy library (Enrichr API) to compute pathway over-representation,
    saves results in JSON format, extracts genes from specific enriched
    pathways, and generates customizable bar plots of the top enriched
    pathways.

Core Features:
    • Run enrichment analysis using Reactome, KEGG, GO Biological Process,
      or any Enrichr-supported gene set collection.
    • Save enrichment output as a JSON file for reproducibility and
      downstream processing.
    • Retrieve top enriched pathways sorted by adjusted p-value.
    • Visualize the most significant pathways using horizontal bar plots
      (−log10 adjusted p-value).
    • Extract gene members associated with any enriched pathway.

Classes:
    GeneListEnrichment:
        Handles all enrichment operations, visualization, and export.

Inputs:
    • List of gene symbols (string identifiers)
    • Optional organism (“Human” by default)
    • User-specified pathway database (e.g., Reactome_2022, KEGG_2021, GO_Biological_Process_2021)

Outputs:
    • JSON file containing full enrichment results
    • Optional PNG/PDF bar plot of top enriched pathways
    • Programmatic access to enriched terms and pathway genes

Example:
    See __main__ section for demonstration with a test gene list.

Dependencies:
    pandas, matplotlib, seaborn, numpy, gseapy, json, os
"""


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import gseapy as gp  # Python library for gene set enrichment
import numpy as np
import json

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
    
    def run_enrichment(self, database='Reactome_2022', fdr_threshold=0.05):
        """
        Perform enrichment using gseapy.enrichr.
        database: pathway database (Reactome, KEGG, GO_Biological_Process)
        fdr_threshold: adjusted p-value cutoff
        """
        print(f"Running enrichment on {len(self.gene_list)} genes using {database}")
        enr = gp.enrichr(
            gene_list=self.gene_list,
            gene_sets=database,
            organism=self.organism,
            outdir=self.outdir,
            cutoff=fdr_threshold
        )
        self.enrich_res = enr.results
        # Save results as JSON for downstream use
        json_file = os.path.join(self.outdir, "enrichment_results.json")
        self.enrich_res.to_json(json_file, orient="records", indent=2)
        print(f"Saved enrichment results to {json_file}")
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
            palette='viridis',
            hue=None,
            dodge=False
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
    # Example gene list (≈100 genes for testing)
    gene_list = [
        "TP53","BRCA1","BRCA2","LRRK2","SNCA","APP","EGFR","VEGFA","PTEN","AKT1",
        "MYC","CDK2","CDK4","RB1","BCL2","MAPK1","MAPK3","JUN","FOS","STAT3"
    ]
    
    enrichment = GeneListEnrichment(gene_list)
    enrichment.run_enrichment(database='Reactome_2022')
    
    # Print top pathways
    print(enrichment.top_pathways(5))
    
    # Plot top pathways
    enrichment.plot_top_pathways(10, save_path="results/enrichment/top_pathways.png")
    
    # Get genes for a specific pathway
    genes = enrichment.get_pathway_genes("Parkinsons disease")
    print("Genes in Parkinsons disease pathway:", genes)

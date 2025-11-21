import pandas as pd
import gseapy as gp

# ----------------------------
# 1. Load differential expression results
# ----------------------------
diff_file = "output/diff_expression.tsv"
gene_col = "GeneID"        # replace with your gene/protein column
fc_col = "logFC"
fdr_col = "FDR"

df = pd.read_csv(diff_file, sep="\t")

# ----------------------------
# 2. Option A: ORA (Over-Representation Analysis)
# ----------------------------
# Select significant genes
sig_genes = df[(df[fdr_col] < 0.05) & (abs(df[fc_col]) >= 1.5)][gene_col].tolist()

# Run enrichment using Reactome (other libraries like KEGG also possible)
ora_results = gp.enrichr(
    gene_list=sig_genes,
    gene_sets=['Reactome_2016'],  # Or KEGG_2019_Human, GO_Biological_Process_2021, etc.
    organism='Human',
    description='ORA_results',
    outdir='output/ora',  # output directory
    cutoff=0.05            # significance threshold
)
print("ORA enrichment finished. Results saved in output/ora")

# ----------------------------
# 3. Option B: GSEA (Ranked List)
# ----------------------------
# Prepare ranked list (all genes) as DataFrame
ranked_df = df[[gene_col, fc_col]].dropna()
ranked_df = ranked_df.sort_values(by=fc_col, ascending=False)

# Save ranked list (optional)
ranked_df.to_csv("output/ranked_genes.rnk", sep="\t", index=False, header=False)

# Run GSEA preranked
pre_res = gp.prerank(
    rnk=ranked_df,                  # ranked dataframe
    gene_sets='Reactome_2016',      # Reactome or other gene sets
    outdir='output/gsea',
    min_size=10,
    max_size=500,
    permutation_num=100,            # reduce for speed
    seed=42
)
print("GSEA preranked finished. Results saved in output/gsea")

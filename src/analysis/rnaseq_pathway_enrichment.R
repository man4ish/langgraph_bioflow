# rnaseq_pathway_enrichment.R
# Author: Manish
# Date: 2025-11-18
# Description: Take Kallisto gene expression matrix and perform pathway enrichment analysis using clusterProfiler and GSEA

# -----------------------------
# Load libraries
# -----------------------------
library(clusterProfiler)
library(org.Hs.eg.db)     # Human gene annotation
library(DOSE)             # For GSEA
library(enrichplot)
library(ggplot2)
library(dplyr)
library(readr)

# -----------------------------
# Input arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("Usage: Rscript pathway_enrichment.R <gene_expression_matrix.tsv>")
}

gene_matrix_file <- args[1]

# -----------------------------
# Load gene expression matrix
# -----------------------------
gene_matrix <- read_tsv(gene_matrix_file)

# Assume first column is gene IDs, remaining columns are samples
gene_ids <- gene_matrix[[1]]
expr_data <- gene_matrix[,-1]  # numeric expression

# -----------------------------
# Preprocessing for enrichment
# -----------------------------
# Example: Compute mean expression across samples
gene_mean <- rowMeans(expr_data)
names(gene_mean) <- gene_ids

# Sort decreasing for GSEA
gene_list <- sort(gene_mean, decreasing = TRUE)

# Map gene symbols to Entrez IDs
gene_entrez <- bitr(names(gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_list <- gene_list[gene_entrez$SYMBOL]
names(gene_list) <- gene_entrez$ENTREZID

# -----------------------------
# GSEA analysis
# -----------------------------
gsea_res <- gseGO(
  geneList = gene_list,
  ont = "BP",              # Biological Process
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Save GSEA table
gsea_table <- as.data.frame(gsea_res)
write_tsv(gsea_table, "GSEA_results.tsv")

# -----------------------------
# Over-Representation Analysis (ORA)
# -----------------------------
# Take top 100 highly expressed genes as example
top_genes <- names(sort(gene_list, decreasing=TRUE)[1:100])

ora_res <- enrichGO(
  gene = top_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Save ORA results
ora_table <- as.data.frame(ora_res)
write_tsv(ora_table, "ORA_results.tsv")

# -----------------------------
# Visualization
# -----------------------------
pdf("Enrichment_Plots.pdf", width=10, height=8)
dotplot(gsea_res, showCategory=20, title="GSEA: Top 20 GO Terms")
dotplot(ora_res, showCategory=20, title="ORA: Top 20 GO Terms")
dev.off()

message("Pathway enrichment analysis complete. Outputs:")
message(" - GSEA_results.tsv")
message(" - ORA_results.tsv")
message(" - Enrichment_Plots.pdf")

#' -----------------------------------------------------------------------------
#' Module: RNA-seq Pathway Enrichment Analysis
#'
#' Author: Manish
#' Date: 2025-11-18
#'
#' Description:
#'   This script performs pathway enrichment analysis from a Kallisto-generated
#'   gene expression matrix. It supports both Gene Set Enrichment Analysis (GSEA)
#'   and Over-Representation Analysis (ORA) using the clusterProfiler package.
#'   The script generates tabular results and visualizations of top enriched pathways.
#'
#' Workflow:
#'   1. Load gene expression matrix (TSV format) with gene symbols as row identifiers.
#'   2. Preprocess expression data (e.g., compute mean expression across samples).
#'   3. Map gene symbols to Entrez IDs for compatibility with enrichment functions.
#'   4. GSEA:
#'        - Uses ranked gene list sorted by mean expression.
#'        - Performs enrichment for GO Biological Process (BP).
#'        - Saves results as "GSEA_results.tsv".
#'   5. ORA:
#'        - Selects top highly expressed genes (default top 100).
#'        - Performs GO BP enrichment.
#'        - Saves results as "ORA_results.tsv".
#'   6. Visualization:
#'        - Dotplots of top enriched pathways for GSEA and ORA.
#'        - Saves plots in "Enrichment_Plots.pdf".
#'
#' Input:
#'   - Command-line argument: path to gene expression TSV file.
#'
#' Output:
#'   - GSEA_results.tsv       : Tabular results from GSEA
#'   - ORA_results.tsv        : Tabular results from ORA
#'   - Enrichment_Plots.pdf   : Dotplots of top pathways
#'
#' Dependencies:
#'   - clusterProfiler
#'   - org.Hs.eg.db
#'   - DOSE
#'   - enrichplot
#'   - ggplot2
#'   - dplyr
#'   - readr
#'
#' Usage:
#'   Rscript rnaseq_pathway_enrichment.R <gene_expression_matrix.tsv>
#' -----------------------------------------------------------------------------


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

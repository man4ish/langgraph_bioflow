#'
#' Single-cell RNA-seq Analysis and Pathway Enrichment Pipeline
#'
#' This script performs an end-to-end analysis of single-cell RNA-seq (scRNA-seq)
#' data processed with CellRanger. It builds a Seurat object, performs quality
#' control, normalization, dimensionality reduction, clustering, marker gene
#' identification, and functional enrichment analysis.
#'
#' Workflow Summary:
#'
#' 1. **Data Loading**
#'    - Reads filtered CellRanger matrices (`filtered_feature_bc_matrix`).
#'    - Creates a Seurat object for downstream analysis.
#'
#' 2. **Quality Control**
#'    - Computes mitochondrial percentage per cell.
#'    - Generates QC violin plots (nFeature_RNA, nCount_RNA, percent.mt).
#'    - Filters low-quality cells based on feature count and mitochondrial content.
#'
#' 3. **Normalization & Feature Selection**
#'    - Log-normalization of counts.
#'    - Identification of highly variable genes (HVGs).
#'    - Scaling of variable genes for PCA.
#'
#' 4. **Dimensionality Reduction**
#'    - Runs PCA on HVGs.
#'    - Saves PCA plot.
#'    - Computes UMAP embeddings for visualization.
#'
#' 5. **Clustering**
#'    - Constructs nearest-neighbor graph.
#'    - Performs clustering at user-specified resolution.
#'    - Saves UMAP plot with cluster labels.
#'
#' 6. **Cluster Marker Identification**
#'    - Identifies marker genes for each cluster using differential expression.
#'    - Saves full marker list and heatmap of top markers per cluster.
#'
#' 7. **Differential Expression Example**
#'    - Includes example contrast (cluster 1 vs cluster 2).
#'
#' 8. **Pathway Enrichment Analysis**
#'    - Converts gene symbols to Entrez IDs.
#'    - Runs:
#'        - GO Biological Process enrichment.
#'        - KEGG pathway enrichment.
#'    - Saves enrichment results and barplots of top enriched pathways.
#'
#' Outputs:
#' - QC plots, PCA plot, UMAP plot
#' - Marker lists, top-marker heatmap
#' - Differential expression results
#' - GO and KEGG enrichment tables and barplots
#' - Final Seurat object (`seurat_obj.rds`)
#'
#' Dependencies:
#' - Seurat, ggplot2, dplyr, patchwork
#' - clusterProfiler, org.Hs.eg.db, enrichplot
#'
#' This script serves as a template for scRNA-seq exploratory analysis and
#' downstream pathway interpretation using standard Seurat and clusterProfiler
#' workflows.
#'


# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)      # For combining plots
library(clusterProfiler)
library(org.Hs.eg.db)   # Human genome annotation
library(enrichplot)     # For enrichment visualization

# =====================
# Define file paths
# =====================
cellranger_directory <- "path/to/CellRanger/output/directory"
output_directory      <- "path/to/output/directory"
dir.create(output_directory, showWarnings = FALSE)

# =====================
# Load data into Seurat
# =====================
seurat_data <- Read10X(data.dir = file.path(cellranger_directory, "filtered_feature_bc_matrix"))
seurat_obj  <- CreateSeuratObject(counts = seurat_data, project = "SingleCellRNASeq")

# =====================
# Quality control
# =====================
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
qc_plots <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file.path(output_directory, "QC_plots.pdf"), plot = qc_plots)

# Filter cells based on QC metrics
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# =====================
# Normalization and variable features
# =====================
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = all_of(VariableFeatures(seurat_obj)))

# =====================
# Dimensional reduction
# =====================
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
pca_plot <- DimPlot(seurat_obj, reduction = "pca") + ggtitle("PCA Plot")
ggsave(file.path(output_directory, "PCA_plot.pdf"), plot = pca_plot)

# =====================
# Clustering
# =====================
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)
umap_plot <- DimPlot(seurat_obj, reduction = "umap") + ggtitle("UMAP Plot")
ggsave(file.path(output_directory, "UMAP_plot.pdf"), plot = umap_plot)

# Save Seurat object
saveRDS(seurat_obj, file = file.path(output_directory, "seurat_obj.rds"))

# =====================
# Cluster markers
# =====================
cluster_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cluster_markers, file = file.path(output_directory, "cluster_markers.csv"))

# Heatmap of top markers
top_markers <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heatmap_plot <- DoHeatmap(seurat_obj, features = top_markers$gene) + ggtitle("Heatmap of Top Markers")
ggsave(file.path(output_directory, "Heatmap_top_markers.pdf"), plot = heatmap_plot)

# =====================
# Differential expression (example: cluster 1 vs 2)
# =====================
diff_exp <- FindMarkers(seurat_obj, ident.1 = 1, ident.2 = 2)
write.csv(diff_exp, file = file.path(output_directory, "diff_exp_cluster_1_vs_2.csv"))

# =====================
# Pathway enrichment analysis
# =====================
# Use cluster_markers for enrichment
genes <- unique(cluster_markers$gene)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO Biological Process enrichment
go_results <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
write.csv(as.data.frame(go_results), file = file.path(output_directory, "GO_enrichment.csv"))

# KEGG pathway enrichment
kegg_results <- enrichKEGG(
  gene         = entrez_ids$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)
write.csv(as.data.frame(kegg_results), file = file.path(output_directory, "KEGG_enrichment.csv"))

# Optional: Visualize top pathways
pdf(file.path(output_directory, "GO_enrichment_barplot.pdf"))
barplot(go_results, showCategory = 20)
dev.off()

pdf(file.path(output_directory, "KEGG_enrichment_barplot.pdf"))
barplot(kegg_results, showCategory = 20)
dev.off()

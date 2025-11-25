# -------------------------------------------------------------------------
#         Test clustering analysis on pbmc_small dataset
# -------------------------------------------------------------------------

# global imports
import::from("magrittr", "%>%", .character_only=TRUE)
import::from(
  "Seurat",
  c("FindAllMarkers", "FindMarkers", "Idents<-", "VariableFeatures"),
  .character_only=TRUE
)
import::from("SeuratDisk", "LoadH5Seurat", .character_only=TRUE)

# local imports
import::from(
  "./scripts/export.R",
  c(
    "export_dotplots", "export_embeddings", "export_heatmap", 
    "export_markers", "export_markers_unique", "FindFullMarkers",
    "visualize_markers"
  ),
  .character_only=TRUE
)
import::from("./scripts/ggtheme.R", "theme_custom", .character_only=TRUE)


# -------------------------------------------------------------------------

# set global ggplot2 theme
ggplot2::theme_set(theme_custom())

# path to the folder where the raw data are stored
datadir <- "./data"

# path to the folder where to export figures and save data
savedir <- "./results/pbmc"

# Load Seurat object
Sobject <- LoadH5Seurat(file = file.path(savedir, "pbmc_small.h5Seurat"))

# -------------------------------------------------------------------------

# set colormap (a valid string from the pals package)
cmap <- "alphabet"

# export UMAP embedding with custom coloring based on clustering
res_vec <- colnames(Sobject@meta.data) %>% 
  strsplit(., split="SCT_snn_res.") %>%
  sapply(., function(x) x[2]) %>% 
  na.omit(.) %>%
  as.numeric(.)
export_embeddings(Sobject, res_vec, cmap, savedir, assay="SCT", device = "png")

# -------------------------------------------------------------------------

# choose clustering resolution (the number of clusters)
ident <- "SCT_snn_res.1"
Idents(Sobject) <- ident

# -------------------------------------------------------------------------

# find all relevant markers with constrained parameters
markers <- FindAllMarkers(
  Sobject, only.pos = TRUE, min.pct = 0.1, min.diff.pct = 0.50,
  logfc.threshold = 0.25
)

# save result table
export_markers(markers, ident, 10, savedir)

# save heatmap of top50 genes
top_genes <- head(VariableFeatures(object = Sobject, assay = "SCT"), 50)
export_heatmap(
  Sobject, top_genes, "SCT", TRUE, cmap, "top50", savedir, device = "png"
)

# -------------------------------------------------------------------------

# export table of markers that are mostly expressed in a single cluster
# min_pct = 5 means expressed in at least 5% of the cells of a cluster
# threshold = 1 means expressed in less than 1% of the remaining cells
export_markers_unique(Sobject, ident, savedir, threshold = 1, min_pct = 5)

# -------------------------------------------------------------------------

# get all markers per cluster irregardless of their expression power
markers <- FindFullMarkers(Sobject, ident, savedir)
export_markers(markers, ident, NULL, savedir, prefix = "full_")

# heatmap of top5 genes per cluster
top5_genes <- markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 5, wt = avg_log2FC)
export_heatmap(
  Sobject, top5_genes$gene, "SCT", TRUE, cmap, "top5_perclust", savedir,
  device = "png"
)

# heatmap of top10 genes per cluster
top10_genes <- markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 10, wt = avg_log2FC)
export_heatmap(
  Sobject, top10_genes$gene, "SCT", TRUE, cmap, "top10_perclust", savedir,
  device = "png"
)

#---------------------------------------------------------------------------
# A collection of utility functions for exporting clustering results
# and visualizations
#--------------------------------------------------------------------------

# global imports
import::here(
  "dplyr", 
  c("filter", "group_by", "left_join", "mutate", "rename", "top_n"),
  .character_only=TRUE
)
import::here(
  "ggplot2",
  c(
    "ggsave",
    "guide_axis",
    "scale_color_gradientn",
    "scale_color_manual",
    "scale_fill_gradientn",
    "scale_x_discrete"
    ),
  .character_only=TRUE)
import::here("magrittr", c("%>%", "%<>%"), .character_only=TRUE)
import::here("matrixStats", "rowVars", .character_only=TRUE)
import::here("pals", "plasma", .character_only=TRUE)
import::here(
  "Seurat",
  c("DoHeatmap", "DimPlot", "DotPlot", "FeaturePlot", "VlnPlot"),
  .character_only=TRUE
)
import::here("SeuratObject", "Misc<-", .character_only=TRUE)
# local imports
import::here(
  "statistics.R",
  c("Cluster_Stats_All_Samples", "Percent_Expressing"),
  .character_only=TRUE
)
import::here("ggtheme.R", "theme_custom", .character_only=TRUE)



#' Export missing markers.
#'
#' Given a target list og markers, export those that are missing from the
#' Seurat dataset.
#'
#' @param object the Seurat class object
#' @param features the features to visualize
#' @param savedir the path to the export directory
#' @param prefix the csv file name prefix
#'
#' @return NULL
export_missing_markers <- function(
    object,
    features,
    savedir,
    prefix,
    overwrite = FALSE
){
  # check whether to process or not
  if (!overwrite | !file.exists(file.path(savedir, prefix))){
    # extract missing markers
    feats_out_RNA <- features[!(features %in% rownames(Sobject@assays$RNA))]
    feats_out_SCT <- features[!(features %in% rownames(Sobject))]
    diff_feats <- feats_out_SCT[!(feats_out_SCT %in% feats_out_RNA)]
    check <- Sobject@assays$RNA[rownames(Sobject@assays$RNA) %in% diff_feats,]
    cts <- rep(0, length(diff_feats))
    cts[as.integer(names(table(check@i)))+1] = table(check@i)
    # write to csv file
    sink(file=file.path(savedir, paste0(prefix, ".csv")))
    cat("Markers missing from dataset:")
    write.csv(data.frame(feats=feats_out_RNA))
    cat('\n\n\n')
    cat("Markers discarded by SCTransform:")
    write.csv(data.frame(feats=diff_feats, ncells=cts))
    sink()
    closeAllConnections()
  }
  
  return(NULL)
}


#' Export top markers.
#'
#' Export the top markers find for each cluster in a csv file.
#'
#' @param object the Seurat class object
#' @param features the features to visualize
#' @param ident the cluster identification (labeling)
#' @param savedir the path to the export directory
#' @param size the number of genes to plot per graph
#' @param device the device to export the graph
#'
#' @return NULL
export_markers <- function(
    markers,
    ident,
    ntop,
    savedir,
    prefix = "",
    wt = "avg_log2FC"
  ){
  # check whether to group per cluster or not
  nclust <- levels(markers$cluster)
  grp_markers <- markers %>% group_by(cluster)
  if (!is.null(ntop)) {
    grp_markers <- grp_markers %>%
      top_n(n = ntop, wt = get(wt))
  }
  # create output folder
  subdir <- file.path(
    "resolutions", format(strsplit(ident, split="res.")[[1]][2], nsmall=2))
  ifelse(!dir.exists(file.path(savedir, subdir)),
         dir.create(file.path(savedir, subdir), recursive=TRUE),
         FALSE)
  # write to csv file
  sink(file=paste0(file.path(savedir, subdir, prefix), "markers.csv"))
  for (idx in nclust) {
    cat(paste0("Markers of Cluster ", idx, "\n"))
    write.csv(grp_markers %>% dplyr::filter(cluster == idx))
    cat('\n\n\n')
  }
  sink()
  closeAllConnections()
  # save markers object
  saveRDS(markers, paste0(file.path(savedir, subdir, prefix), "markers.rds"))

  return(NULL)
}


#' Find markers.
#'
#' Find all the markers for each cluster.
#'
#' @param object the Seurat class object
#' @param ident the cluster identification (labeling)
#' @param savedir the path to the export directory
#' @param only.pos the cluster identification (labeling)
#' @param min_pct the minimum percent of cells (in each of the two groups)
#' where the gene is expressed
#' @param min.diff.pct the minimum difference between each of the two group
#' percent of cells where the gene is expressed
#' @param logfc.threshold restrict testing to genes which show, on average, at
#' least X-fold difference (log-scale) between the two groups of cells
#' @param overwrite whether to overwrite the markers data or not 
#' 
#' @return A DataFrame with a ranked list of putative markers as rows, and
#' associated statistics as columns
FindFullMarkers <- function(
    object,
    ident,
    savedir,
    only.pos = TRUE,
    min.pct = 1e-3,
    min.diff.pct = 1e-3,
    logfc.threshold = 1e-3,
    overwrite = FALSE,
    ...
){
  filepath <- file.path(
    savedir,
    "resolutions",
    format(strsplit(ident, split="res.")[[1]][2], nsmall=2),
    "full_markers.rds"
  )
  if (!overwrite & file.exists(filepath)){
    markers <- readRDS(filepath)
  } else {
    markers <- FindAllMarkers(
      object, only.pos = only.pos, min.pct = min.pct, min.diff.pct = min.diff.pct,
      logfc.threshold = logfc.threshold, ...)
  }

  return(markers)
}

#' Export unique markers.
#'
#' Export the markers that are only expressed in a single cluster.
#'
#' @param object the Seurat class object
#' @param ident the cluster identification (labeling)
#' @param savedir the path to the export directory
#' @param threshold upper bound on the percent of cells not belonging to the
#' cluster where the gene is expressed
#' @param min_pct the minimum percent of cells of the cluster in which the
#' gene is expressed
#'
#' @return NULL
export_markers_unique <- function(
    object,
    ident,
    savedir,
    threshold,
    min_pct
){
  # table of relative marker percent expression
  pct_tab <- Percent_Expressing(object, features = rownames(object))
  # create output folder
  subdir <- file.path(
    "resolutions", format(strsplit(ident, split="res.")[[1]][2], nsmall=2))
  ifelse(!dir.exists(file.path(savedir, subdir)),
         dir.create(file.path(savedir, subdir)),
         FALSE)
  # statistic per cluster
  stats <- Cluster_Stats_All_Samples(object)
  # write to csv file
  sink(file=file.path(savedir, subdir, "unique_markers.csv"))
  cat("Percentages per cluster\n")
  write.csv(stats[,1:3])
  cat('\n\n\n')
  cat(paste0("Threshold for unique markers: ", threshold))
  cat('\n\n')
  # loop over clusters
  nrows <- nrow(stats)
  for(idx in 1:ncol(pct_tab)){
    # get table of absolute marker for other cells
    rel_pct <- stats[-c(idx,nrows),2] / sum(stats[-c(idx,nrows),2])
    abs_tab <- as.data.frame(t(t(pct_tab[,-idx]) * rel_pct))
    cat(paste0("Unique markers of Cluster ", idx, "\n"))
    c_table <- cbind(pct_tab[,idx], abs_tab)
    names(c_table) <- c(names(pct_tab)[idx], names(pct_tab)[-idx])
    cond <- (rowSums(abs_tab) < threshold) & 
      (pct_tab[,idx] >= min_pct)
    write.csv(c_table[cond,])
    cat('\n\n\n')
  }
  sink()
  closeAllConnections()

  return(NULL)
}


#' Export top markers from two-way comparison.
#'
#' Export the top markers find when comparing two clusters. Both the bottom and
#' top markers are exported to respect the two way comparison (i.e. cluster A vs 
#' B and cluster B vs A).
#'
#' @param markers the matrix containing the list of putative markers
#' @param ident the cluster identification (labeling)
#' @param ntop the number of top markers
#' @param savedir the path to the export directory
#' @param wt the variable use for ordering the top markers
#' @param prefix the csv file name prefix
#'
#' @return NULL
export_markers_comp <- function(
    markers,
    ident,
    ntop,
    savedir,
    prefix,
    wt = "avg_log2FC"
){
  # check whether to group per cluster or not
  grp_markers <- rbind(
    markers %>% top_n(n = ntop, wt = get(wt)),
    markers %>% top_n(n = -ntop, wt = get(wt))
  )
  # create output folder
  subdir <- file.path(
    "resolutions", format(strsplit(ident, split="res.")[[1]][2], nsmall=2))
  ifelse(!dir.exists(file.path(savedir, subdir)),
         dir.create(file.path(savedir, subdir), recursive=TRUE),
         FALSE)
  # write to csv file
  sink(file=paste0(file.path(savedir, subdir, prefix), "markers.csv"))
  cat("Top markers")
  cat("\n")
  write.csv(grp_markers)
  sink()
  closeAllConnections()

  return(NULL)
}


#' Export dotplots.
#'
#' Export dotplots graphs for a huge list of genes
#'
#' @param object the Seurat class object
#' @param features the features to visualize
#' @param ident the cluster identification (labeling)
#' @param savedir the path to the export directory
#' @param size the number of genes to plot per graph
#' @param device the device to export the graph
#'
#' @return NULL
export_dotplots <- function(
    object,
    features,
    ident,
    savedir,
    prefix,
    size = 10,
    device = "svg",
    width = 12,
    height = 8,
    ...
){
  # define colormap
  na_cutoff <- 1e-9
  na_color <- "lightgray"
  cmap_vir <- rev(plasma(256))
  # create output folder
  subdir <- file.path(
    "resolutions", format(strsplit(ident, split="res.")[[1]][2], nsmall=2))
  ifelse(!dir.exists(file.path(savedir, subdir, "dotplots", prefix)),
         dir.create(file.path(savedir, subdir, "dotplots", prefix), recursive=TRUE),
         FALSE)
  # keep only the relevant features
  features <- features[features %in% rownames(Sobject)]
  blocks <- seq(1, length(features), by=size)
  if (blocks[length(blocks)] < length(features)) {
    blocks <- c(blocks, length(features))
  }
  # save the dotplots with custom gradient
  filepath <- file.path(savedir, subdir, "dotplots", prefix, "fig")
  for (it in 1:(length(blocks)-1)) {
    gg <- suppressWarnings(suppressMessages(
      DotPlot(object, features = features[blocks[it]:(blocks[it+1]-1)]) & 
        theme_custom() &
        scale_x_discrete(guide = guide_axis(n.dodge = 2)) & 
        scale_color_gradientn(
          colors = cmap_vir, limits = c(na_cutoff, NA), na.value = na_color)
    ))
    ggsave(file=paste0(filepath, it, ".", device), device = device, 
           width = width, height = height, gg, ...)
  }

  return(NULL)
}


#' Color embedding per resolution.
#'
#' Export the embedding result of UMAP for several clustering resolutions.
#'
#' @param object the Seurat class object
#' @param resolutions the vector of clustering resolutions
#' @param colormap a valid colormap name from the "pals" library
#' @param savedir the path to the directory where the plots will be saved
#' @param assay the name of the Assay to use
#' @param device the device to export the graph
#'
#' @return NULL
export_embeddings <- function(
    object,
    resolutions,
    colormap,
    savedir,
    assay = "SCT",
    device = "svg",
    width = 8,
    height = 8,
    ...
){
  # export colormap
  cmap <- get(colormap, envir = loadNamespace("pals"))
  # loop over vector of resolutions
  for (idx in 1:length(resolutions)){
    val <- resolutions[idx]
    colname <- paste0(assay, "_snn_res.", as.character(val))
    cols <- unname(cmap(nlevels(object@meta.data[,colname])))
    filepath <- file.path(savedir, paste0("umap_res", format(val, nsmall=2)))
    gg <- suppressMessages(
      DimPlot(object, reduction = "umap", cols=cols, group.by = colname)
    )
    ggsave(file=paste0(filepath, ".", device), device = device, 
           width = width, height = height, gg, ...)
  }

  return(NULL)
}


#' Features visualization.
#'
#' Export three types of features visualization: violin plot, feature expression
#' in UMAP embedding and dot plot.
#'
#' @param object the Seurat class object
#' @param features the features to visualize
#' @param ident the cluster identification (labeling)
#' @param colormap a valid colormap name from the "pals" library
#' @param prefix the graph name prefix
#' @param savedir the path to the export directory
#' @param device the device to export the graph
#'
#' @return NULL
visualize_markers <- function(
    object,
    features,
    ident,
    colormap,
    prefix,
    savedir,
    device = "svg",
    width = 8,
    height = 8,
    ...
){
  # catch extra arguments
  extra_args <- list(...)
  # define colormap
  cmap <- get(colormap, envir = loadNamespace("pals"))
  cols <- unname(cmap(nlevels(object@active.ident)))
  na_cutoff <- 1e-9
  na_color <- "lightgray"
  cmap_vir <- rev(plasma(256))
  
  # create output folder
  subdir <- file.path(
    "resolutions", format(strsplit(ident, split="res.")[[1]][2], nsmall=2))
  ifelse(!dir.exists(file.path(savedir, subdir)),
         dir.create(file.path(savedir, subdir), recursive=TRUE),
         FALSE)
  # save the violin plot
  for (feat in features){
    filepath <- file.path(savedir, subdir, feat)
    vln_args <- c(
      list("object" = object, "features" = feat, "pt.size" = 0, "cols" = cols),
      extra_args[names(extra_args) %in% names(as.list(args(VlnPlot)))]
    )
    gg <- suppressMessages(do.call(VlnPlot, vln_args))
    ggsave(file=paste0(filepath, "_vlnplot", ".", device), device = device,
           width = width, height = width, gg)
  }
  # save the feature plot with custom gradient
  max_cmap <- SeuratObject::FetchData(
    object = object,
    vars = features,
    slot = "data"
  ) %>% max()
  for (feat in features){
    filepath <- file.path(savedir, subdir, feat)
    feat_args <- c(
      list("object" = object, "features" = feat),
      extra_args[names(extra_args) %in% names(as.list(args(FeaturePlot)))]
    )
    gg <- suppressMessages(
      do.call(FeaturePlot, feat_args) & 
        scale_color_gradientn(
          colors = cmap_vir, limits = c(na_cutoff, max_cmap), na.value = na_color)
    )
    ggsave(file=paste0(filepath, "_featplot", ".", device), device = device,
           width = width, height = height, gg)
  }
  # save the dot plot with custom gradient
  filepath <- file.path(savedir, subdir, prefix)
  dot_args <- c(
    list("object" = object, "features" = features),
    extra_args[names(extra_args) %in% names(as.list(args(DotPlot)))]
  )
  gg <- suppressWarnings(suppressMessages(
    do.call(DotPlot, dot_args) &
      theme_custom() &
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) & 
          scale_color_gradientn(
            colors = cmap_vir, limits = c(na_cutoff, NA), na.value = na_color)
  ))
  ggsave(file=paste0(filepath, "_dotplot", ".", device), device = device,
         width = width, height = height, gg)

  return(NULL)
}


#' Heatmap with cluster.
#'
#' Generates an expression heatmap for given cells and features
#' for each cluster.
#'
#' @param object the Seurat class object
#' @param features the features to visualize
#' @param assay the name of the Assay to use
#' @param legend whether to plot the legend or not
#' @param colormap a valid colormap name from the "pals" library
#' @param prefix the graph name prefix
#' @param savedir the path to the export directory
#' @param device the device to export the graph
#'
#' @return NULL
export_heatmap <- function(
    object,
    features,
    assay,
    legend,
    colormap,
    prefix,
    savedir,
    device = "svg",
    width = 12,
    height = 12,
    ...
){
  # define colormap
  cmap <- get(colormap, envir = loadNamespace("pals"))
  cols <- unname(cmap(nlevels(object@active.ident)))
  # create output folder
  subdir <- file.path(
    "resolutions", format(strsplit(ident, split="res.")[[1]][2], nsmall=2))
  ifelse(!dir.exists(file.path(savedir, subdir)),
         dir.create(file.path(savedir, subdir)),
         FALSE)
  # save the heatmap
  filepath <- file.path(savedir, subdir, prefix)
  gg <- suppressWarnings(suppressMessages(
    DoHeatmap(
      object, assay=assay, features=features, group.colors=cols
      ) &
      scale_color_manual(
        values = cols, name = "Identity", na.translate = FALSE
      )
  ))
  if (!legend) { gg <- gg + NoLengend()}
  ggsave(file=paste0(filepath, "_heatmap", ".", device), device = device,
           width = width, height = height, gg, ...)

  return(NULL)
}

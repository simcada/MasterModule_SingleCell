#------------------------------------------------------------------
# subset of scCustomize statistical functions 
#------------------------------------------------------------------

# global imports
import::here("dplyr", c("left_join", "rename"), .character_only=TRUE)
import::here("janitor", "adorn_totals", .character_only=TRUE)
import::here("magrittr", "%>%", .character_only=TRUE)
import::here("matrixStats", "rowVars", .character_only=TRUE)
import::here("glue", "glue_collapse", .character_only=TRUE)
import::here(
  "Seurat",
  c(
    "Assays", "CellsByIdentities", "DefaultAssay", "Idents", "FetchData",
    "GetAssayData"
  ),
  .character_only=TRUE
)
import::here(
  "tibble",
  c("rownames_to_column", "column_to_rownames"),
  .character_only=TRUE
)
import::here("tidyr", "pivot_wider", .character_only=TRUE)

# local imports
import::here(
  "utilities.R",
  c("%||%", "Gene_Present", "Is_Seurat", "PercentAbove_Seurat"),
  .character_only=TRUE
)

#' Calculate Cluster Stats
#'
#' Calculates both overall and per sample cell number and percentages per cluster based on orig.ident
#'
#' @param seurat_object Seurat object name.
#' @param group_by_var meta data column to classify samples (default = "orig.ident").
#'
#' @importFrom dplyr left_join rename
#' @importFrom janitor adorn_totals
#' @importFrom magrittr "%>%"
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_wider
#'
#' @return A Data Frame
#'
#' @export
#'
#' @concept stats
#'
#' @examples
#' \dontrun{
#' stats <- Cluster_Stats_All_Samples(seurat_object = object, group_by_var = "orig.ident")
#' }
#'
Cluster_Stats_All_Samples <- function(
    seurat_object,
    group_by_var = "orig.ident"
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)
  
  # Check on meta data column
  possible_meta_col <- colnames(seurat_object@meta.data)
  if (!group_by_var %in% possible_meta_col) {
    stop(paste0('"', group_by_var, '"', " was not found in meta.data slot of Seurat Object."))
  }
  
  # Extract total percents
  total_percent <- prop.table(x = table(seurat_object@active.ident)) * 100
  total_percent <- data.frame(total_percent) %>%
    rename(Cluster = Var1)
  
  # Extract total cell number per cluster across all samples
  total_cells <- table(seurat_object@active.ident) %>%
    data.frame() %>%
    rename(Cluster = Var1, Number = Freq)
  
  # Cluster overall stats across all animals
  cluster_stats <- suppressMessages(left_join(total_cells, total_percent))
  
  # Extract cells per metadata column per cluster
  cells_per_cluster_2 <- table(seurat_object@active.ident, seurat_object@meta.data[, group_by_var])
  cells_per_cluster_2 <- data.frame(cells_per_cluster_2) %>%
    rename(Cluster = Var1, group_by_var = Var2, cell_number = Freq)
  
  cells_per_cluster_2 <- cells_per_cluster_2 %>%
    pivot_wider(names_from = group_by_var, values_from = cell_number)
  
  # Merge cells per metadata column per cluster with cluster stats
  cluster_stats_2 <- suppressMessages(left_join(cluster_stats, cells_per_cluster_2))
  
  # Calculate and extract percents of cells per cluster per
  percent_per_cluster_2 <- prop.table(x = table(seurat_object@active.ident, seurat_object@meta.data[, group_by_var]), margin = 2) * 100
  percent_per_cluster_2 <- data.frame(percent_per_cluster_2) %>%
    rename(cluster = Var1, group_by_var = Var2, percent = Freq)
  percent_per_cluster_2 <- percent_per_cluster_2 %>%
    pivot_wider(names_from = group_by_var, values_from = percent) %>%
    column_to_rownames("cluster")
  colnames(percent_per_cluster_2) <- paste(colnames(percent_per_cluster_2), "%", sep = "_")
  
  percent_per_cluster_2 <- percent_per_cluster_2 %>%
    rownames_to_column(var = "Cluster")
  
  # Merge percent cells per metadata column per cluster with cluster stats and add Totals column
  cluster_stats <- suppressMessages(left_join(cluster_stats_2, percent_per_cluster_2)) %>%
    adorn_totals("row")
  return(cluster_stats)
}


#' Calculate percent of expressing cells
#'
#' Calculates the percent of cells that express a given set of features by various grouping factors
#'
#' @param seurat_object Seurat object name.
#' @param features Feature(s) to plot.
#' @param threshold Expression threshold to use for calculation of percent expressing (default is 0).
#' @param group_by Factor to group the cells by.
#' @param split_by Factor to split the groups by.
#' @param entire_object logical (default = FALSE).  Whether to calculate percent of expressing cells
#' across the entire object as opposed to by cluster or by `group_by` variable.
#' @param assay Assay to pull feature data from.  Default is active assay.
#' @param slot Slot to pull feature data for.  Default is "data".
#'
#' @return A Data Frame
#'
#' @references Part of code is modified from Seurat package as used by \code{\link[Seurat]{DotPlot}}
#' to generate values to use for plotting.  Source code can be found here:
#' (https://github.com/satijalab/seurat/blob/4e868fcde49dc0a3df47f94f5fb54a421bfdf7bc/R/visualization.R#L3391) (Licence: GPL-3).
#'
#' @export
#'
#' @concept stats
#'
#' @examples
#' \dontrun{
#' percent_stats <- Percent_Expressing(seurat_object = object, features = "Cx3cr1", threshold = 0)
#' }
#'
Percent_Expressing <- function(
    seurat_object,
    features,
    threshold = 0,
    group_by = NULL,
    split_by = NULL,
    entire_object = FALSE,
    slot = "data",
    assay = NULL
) {
  # Check Seurat
  Is_Seurat(seurat_object = seurat_object)
  
  # set assay (if null set to active assay)
  assay <- assay %||% DefaultAssay(object = seurat_object)
  
  # Check features exist in object
  features_list <- Gene_Present(data = seurat_object, gene_list = features, print_msg = FALSE, case_check = TRUE, seurat_assay = assay)[[1]]
  
  # Check group_by is in object
  if (!is.null(x = group_by)) {
    possible_groups <- colnames(seurat_object@meta.data)
    if (!group_by %in% possible_groups) {
      stop(paste0("Grouping variable: ", '"', group_by, '"', " was not found in Seurat Object"))
    }
  }
  
  # Check split_by is in object
  if (!is.null(x = split_by)) {
    possible_groups <- colnames(seurat_object@meta.data)
    if (!split_by %in% possible_groups) {
      stop(paste0("Splitting variable: ", '"', split_by, '"', " was not found in Seurat Object"))
    }
  }
  
  # Pull Expression Info
  cells <- unlist(x = CellsByIdentities(object = seurat_object, idents = NULL))
  expression_info <- FetchData(object = seurat_object, vars = features_list, cells = cells, slot = slot)
  
  # Add grouping variable
  if (entire_object) {
    expression_info$id <- "All_Cells"
  } else {
    expression_info$id <- if (is.null(x = group_by)) {
      paste("clust", Idents(object = seurat_object)[cells, drop = TRUE])
    } else {
      seurat_object[[group_by, drop = TRUE]][cells, drop = TRUE]
    }
  }
  if (!is.factor(x = expression_info$id)) {
    expression_info$id <- factor(x = expression_info$id)
  }
  id.levels <- levels(x = expression_info$id)
  expression_info$id <- as.vector(x = expression_info$id)
  
  # Split data if split.by is true
  if (!is.null(x = split_by)) {
    splits <- seurat_object[[split_by, drop = TRUE]][cells, drop = TRUE]
    expression_info$id <- paste(expression_info$id, splits, sep = '_')
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  
  # Calculate percent expressing
  percent_expressing <- lapply(
    X = unique(x = expression_info$id),
    FUN = function(ident) {
      data.use <- expression_info[expression_info$id == ident, 1:(ncol(x = expression_info) - 1), drop = FALSE]
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove_Seurat, threshold = threshold)
      return(list(pct.exp = pct.exp))
    }
  )
  names(x = percent_expressing) <- unique(x = expression_info$id)
  
  # Convert & return data.frame
  row_dim_names <- features_list
  col_dim_names <- names(percent_expressing)
  mat_dims <- list(row_dim_names, col_dim_names)
  final_df <- data.frame(matrix(unlist(percent_expressing), nrow = length(features_list), byrow = FALSE, dimnames = mat_dims), stringsAsFactors = FALSE)
  return(final_df)
}


#' PCA variance.
#'
#' Calculate the fraction of total variance explained by each PCA component.
#'
#' @param object the Seurat class object
#' @param assay the name of the Assay to use
#'
#' @return return a vector containing the fractions of the total variance.
varPCA <- function(
    object,
    assay="SCT"
) {
  # get the total variance
  tot_var <- sum(rowVars(object[[assay]]@scale.data))
  # compute the explained variance
  eig_vals <- object[["pca"]]@stdev^2
  res <- eig_vals / tot_var
  
  return(res)
}


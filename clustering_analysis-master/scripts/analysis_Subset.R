# ------------------------------------------------------------------------#
#                 Analysis on a Seurat Object "Name.h5Seurat"             #
# ------------------------------------------------------------------------#
###########################################################################

# -------------------------------------------------------------------------
# Import packages
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

library(Seurat)
library(dplyr)
library(tidyverse)

# set global ggplot2 theme
#maybe need to re isntall ggplot
#renv::install("ggplot2@3.4.1")
#renv::record("ggplot2@3.4.1")
ggplot2::theme_set(theme_custom())


# -------------------------------------------------------------------------
# Change Path accordingly to yours
# -------------------------------------------------------------------------


# path to the folder where the raw data are stored
datadir <- "/Users/groot/Downloads/ForMaster/share/data"

# path to the folder where to export figures and save data
savedir <- "/Users/groot/Downloads/ForMaster/clustering_analysis-master/results"

# -------------------------------------------------------------------------
# Load Seurat Object
# -------------------------------------------------------------------------

# Load Seurat object
Sobject <- LoadH5Seurat(file = file.path(savedir, "Kurmangaliyev_Subset.h5Seurat"))

#check metadata of your Seurat object
colnames(Sobject@meta.data)

### if other dataset are load, remove them to save space
#rm(Sobject_W1118)

# -------------------------------------------------------------------------
# Choose clustering resolution 
# -------------------------------------------------------------------------

# choose clustering resolution (the number of clusters), for example here 3.5
ident <- "SCT_snn_res.3.5"
Idents(Sobject) <- ident

# -------------------------------------------------------------------------
# Look at UMAP graph
# -------------------------------------------------------------------------

#Look at umap, group by cell type
gg <- DimPlot(Sobject, group.by = c("type"))
ggsave(file=paste0(savedir, "Subset_type_UMAP.jpeg"), device = "png",  width = 8, height = 8, gg)

#Look at umap, group by subtypes
gg <- DimPlot(Sobject, group.by = c("subtype"))
ggsave(file=paste0(savedir, "Subset_subtype_UMAP.jpeg"), device = "png",  width = 8, height = 8, gg)

#Look at umap, group by time point
gg <- DimPlot(Sobject, group.by = c("time"))
ggsave(file=paste0(savedir, "Subset_time_UMAP.jpeg"), device = "png",  width = 8, height = 8, gg)

# -------------------------------------------------------------------------
# Look at gene expression
# -------------------------------------------------------------------------

#heatmap top3 genes highly expressed
top3_genes <- markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(n = 3, wt = avg_log2FC)
cmap <- get(colormap, envir = loadNamespace("pals"))
cols <- unname(cmap(nlevels(as.factor(Sobject@meta.data$type))))
gg <- Seurat::DoHeatmap(Sobject, assay="SCT", features=top3_genes$gene, group.by="type", group.colors=cols)
ggsave(file=paste0(savedir, "Top3_heatmap.jpeg"), device = "png",  width = 8, height = 8, gg)

# visualize specific genes expression
features <- c("dac","bi","grn") #change with genes of interest
visualize_markers(
  Sobject, group.by = c("type"), features, ident, cmap, "dac_bi_grn_expression", savedir, device = "png"
) #can change group.by by subtype, time...ect... 

# -------------------------------------------------------------------------
# Look at gene through time
# -------------------------------------------------------------------------

# Choose gene
gene <- "shakB"

# Expression + metadata
data <- FetchData(Sobject, vars = c("type", "time", gene))
data$expressed <- ifelse(data[[gene]] > 0, 1, 0)

# Group and calculate dot size (percentage) and dot color (avg expr)
summary_df <- data %>%
  group_by(type, time) %>%
  summarise(
    pct_expr = mean(expressed) * 100,
    avg_expr = mean(get(gene))
  )

# Create custom DotPlot
p <- ggplot(summary_df, aes(x = type, y = time)) +
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_size(range = c(0, 10)) +  # Adjust size scale
  scale_color_gradient(low = "grey", high = "blue") +
  theme_minimal() +
  labs(
    title = paste(gene, "Expression Across Time and Types"),
    x = "Type",
    y = "Time",
    size = "% of cells expressing",
    color = "Average expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = paste0(savedir, "shakB_in_times.png"), plot = p, width = 10, height = 5, dpi = 300,limitsize = FALSE)


##visualize on umap
p <- FeaturePlot(
  Sobject,
  features = gene,
  split.by = "time",  # ou "time"
  cols = c("lightgrey", "blue"),
  pt.size = 0.5
)
ggsave(file = paste0(savedir, gene, "_umap_time.png"),plot = p, width = 20, height = 5, dpi = 300,limitsize = FALSE)


# -------------------------------------------------------------------------
# Look at multiple genes through time
# -------------------------------------------------------------------------

# Define the genes you want to plot
genes_of_interest <- c("ogre","Inx2","shakB")  # Replace with your gene list

# metada
data <- FetchData(Sobject, vars = c("type", "time", genes_of_interest), assay = "RNA")
long_data <- data %>%
  pivot_longer(cols = all_of(genes_of_interest), names_to = "gene", values_to = "expression")
long_data$expressed <- ifelse(long_data$expression > 0, 1, 0)

# Group and summarize
summary_df <- long_data %>%
  group_by(type, time, gene) %>%
  summarise(
    pct_expr = mean(expressed) * 100,
    avg_expr = mean(expression)
  )

# Plot
p <- ggplot(summary_df, aes(x = type, y = gene)) + 
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_size(range = c(0, 8)) +
  scale_color_gradient(low = "grey", high = "blue") +
  facet_wrap(~ time, ncol = 1) +  # One facet per time point (optional)
  theme_minimal() +
  labs(
    title = "Genes Expression through Time",
    x = "type",
    y = "Gene",
    size = "% of cells expressing",
    color = "Average expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = paste0(savedir, "genes_through_times.png"), plot = p, width = 10, height = 15, dpi = 300,limitsize = FALSE)

# -------------------------------------------------------------------------
# Look at multiple genes through time on a specific cell type
# -------------------------------------------------------------------------

# Define the genes you want to plot
genes_of_interest <- c("ogre","Inx2","shakB")  # Replace with your gene list

# metada
data <- FetchData(Sobject, vars = c("type", "time", genes_of_interest))
genes_present <- genes_of_interest[genes_of_interest %in% colnames(data)]

types_to_plot <- c("T4.T5") #define by your cell type
data_filtered <- data %>% 
  filter(type %in% types_to_plot)

long_data <- data_filtered %>%
  pivot_longer(cols = all_of(genes_present), names_to = "gene", values_to = "expression")
long_data$expressed <- ifelse(long_data$expression > 0, 1, 0)

summary_df <- long_data %>%
  group_by(type, time, gene) %>%
  summarise(
    pct_expr = mean(expressed) * 100,
    avg_expr = mean(expression)
  )

# Plot
p <- ggplot(summary_df, aes(x = type, y = gene)) + 
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  scale_size(range = c(0, 8)) +
  scale_color_gradient(low = "grey", high = "blue") +
  facet_wrap(~ time, ncol = 1) +  # One facet per time point (optional)
  theme_minimal() +
  labs(
    title = "Gene Expression in T4.T5",
    x = "type",
    y = "Gene",
    size = "% of cells expressing",
    color = "Average expression"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = paste0(savedir, "genes_T4T5_times.png"), plot = p, width = 10, height = 15, dpi = 300,limitsize = FALSE)



###################################
# -------------------------------------------------------------------------
# Work on only one time point
# -------------------------------------------------------------------------
###################################

#subset
Sobject_24h <- subset(Sobject, subset = time == "24h")

#visualize
DimPlot(Sobject_24h, group.by = c("time"))
DimPlot(Sobject_24h, group.by = c("type"))

#then you can repeat all the previous steps! 





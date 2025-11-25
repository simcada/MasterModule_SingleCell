# ------------------------------------------------------------------------#
#                             Create a Seurat object                      #
# ------------------------------------------------------------------------#
###########################################################################

# -------------------------------------------------------------------------
# Import packages
# -------------------------------------------------------------------------

#charge all packages DO ONLY ONCE
#renv::restore()

#if you want to install a new package and then save new environment:
#renv::install("mojaveazure/seurat-disk")
#renv::snapshot()


# global imports
import::from("magrittr", "%>%", .character_only=TRUE)
import::from(
  "ggplot2", 
  c(
    "ggsave",
    "scale_colour_gradientn",
    "scale_y_continuous", "theme",
    "ylab"
  ),
  .character_only=TRUE
)
import::from("ggraph", "guide_edge_colourbar", .character_only=TRUE)
import::from(
  "Seurat",
  c(
    "CreateSeuratObject", "ElbowPlot", "FindClusters", "FindNeighbors",
    "LabelPoints", "PercentageFeatureSet", "Read10X", "RunPCA",
    "RunUMAP", "SCTransform", "VariableFeaturePlot", "VariableFeatures",
    "VlnPlot"
  ),
  .character_only=TRUE
)
import::from("SeuratObject", "Misc<-", .character_only=TRUE)
import::from("SeuratDisk", "SaveH5Seurat", .character_only=TRUE)

# relative import
import::from("./scripts/ggtheme.R", "theme_custom", .character_only=TRUE)
import::from("./scripts/statistics.R", "varPCA", .character_only=TRUE)

# set global ggplot2 theme
ggplot2::theme_set(theme_custom())

library(Seurat)
library(dplyr)
library(tidyr)


# -------------------------------------------------------------------------
# Change Path accordingly to yours
# -------------------------------------------------------------------------


# path to read raw data
datadir <- "/Users/groot/Downloads/ForMaster/share/data"

# path to the folder where to export figures and save data
savedir <- "/Users/groot/Downloads/ForMaster/clustering_analysis-master/results/"


# ---------------------------------------------------------------
# DATA LOADING AND PREPROCESSING
# ---------------------------------------------------------------

# load the dataset (and "named it"):
data <- Read10X(
  data.dir = file.path(datadir, "Kurmangaliyev")
)

# initialize the Seurat object with the raw non-normalized data
Sobject <- CreateSeuratObject(
  counts = data, project = "Kurmangaliyev", min.cells = 3, min.features = 200
)
#delete the first raw data to keep only the Seurat object
rm(data)

#add metadata file (that contain annotation of the dataset in our case)
metadata <- read.csv("/Users/groot/Downloads/ForMaster/share/data/Kurmangaliyev/metadata.csv", sep=',')

#add tge row names to our dataset
rownames(metadata) <- colnames(x = Sobject)

#create the full seurat object
Sobject <- SeuratObject::AddMetaData(
  object = Sobject,
  metadata = metadata,
  col.name = colnames(metadata)
)

#####################################################################################
#    to work with a smaller dataset we will subset the big dataset to keep only     #
#    in the W1118 dataset the celltype ("type'): Tm3, Mi9, and T4.T5                #
#####################################################################################

#create the subset Sobject
Sobject_W1118 <- subset(Sobject, subset = set == "W1118")
rm(Sobject)
Sobject <- subset(Sobject_W1118, subset = type %in% c("Tm3", "Mi9", "T4.T5"))

#delete the first Sobject to keep only the subset Seurat object (to save space)
rm(Sobject_W1118)


# ---------------------------------------------------------------
# QUALITY CONTROL AND FILTERING
# ---------------------------------------------------------------

# compute QC metrics
Sobject <- PercentageFeatureSet(
  Sobject, pattern = "(?i)^mt", col.name = "percent.mt"
)

features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

# Extract relevant metadata
df <- Sobject@meta.data %>%
  select(all_of(features), type)

# Convert to long format for ggplot
df_long <- pivot_longer(df, cols = all_of(features),
                        names_to = "feature", values_to = "value")

# Create violin plot grouped by type
gg <- ggplot(df_long, aes(x = type, y = value, fill = type)) +
  geom_violin() +
  facet_wrap(~feature, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none")

# Save plot
ggsave(file.path(savedir, "vlnplot_before_QC.svg"),
       plot = gg, device = "svg", width = 10, height = 4)

#filtering based on QC metrics
Sobject <- subset(Sobject, subset = nFeature_RNA > 500 & nFeature_RNA < 1800)
Sobject <- subset(Sobject, subset = nCount_RNA < 8000 & nCount_RNA > 0)
Sobject <- subset(Sobject, subset = percent.mt < 5)

# other possible filtering based on expert knowledge
#As example: can also select specific cell expression (GFP for example)
#Sobject <- subset(Sobject, subset = GFP != 0)

# export QC metrics after filtering
# Extract relevant metadata
df <- Sobject@meta.data %>%
  select(all_of(features), type)

# Convert to long format for ggplot
df_long <- pivot_longer(df, cols = all_of(features),
                        names_to = "feature", values_to = "value")

# Create violin plot grouped by type
gg <- ggplot(df_long, aes(x = type, y = value, fill = type)) +
  geom_violin() +
  facet_wrap(~feature, scales = "free_y") +
  theme_classic() +
  theme(legend.position = "none")

# Save plot
ggsave(file.path(savedir, "vlnplot_after_QC.svg"),
       plot = gg, device = "svg", width = 10, height = 4)


# ---------------------------------------------------------------
# DATA NORMALIZATION THROUGH SCTRANSFORM
# -> robust alternative to the following common three steps:
#     * NormalizeData,
#     * FindVariableFeatures,
#     * ScaleData
# ---------------------------------------------------------------

# fit negative binomial model
Sobject <- SCTransform(
  Sobject, method = "glmGamPoi", vst.flavor = "v2",
  vars.to.regress = "percent.mt", seed.use = 404
)

# ---------------------------------------------------------------
# HIGHLY VARIABLE GENE SELECTION
# ---------------------------------------------------------------

# get the top 50 genes
top_genes <- head(VariableFeatures(object = Sobject, assay = "SCT"), 50)

# plot the Pearson residual variance against the geometric mean expression
gg_vfp <- VariableFeaturePlot(Sobject, selection.method="sctransform")

# label the top genes
vfp_labels <- LabelPoints(
  plot = gg_vfp, points = top_genes, repel = TRUE, xnudge = 0,
  ynudge = 0, max.overlaps = Inf
)

# export the plot with log10 scaling of y-axis
gg <- vfp_labels & 
  scale_y_continuous(trans='log10') & 
  theme_custom() & 
  theme(legend.position=c(0.5, 0.2))
ggsave(file=file.path(savedir, "Sobject_top_variable_features.svg"), device = "svg",
       width = 10, height = 8, gg)

# ---------------------------------------------------------------
# DIMENSIONALITY REDUCTION : PCA
# ---------------------------------------------------------------

# perform linear dimensional reduction
Sobject <- RunPCA(
  Sobject, features = VariableFeatures(object = Sobject, assay = "SCT"), 
  npcs = 200, seed.use = 42
)

# compute the total variance explained by each PC
Sobject[["pca"]]@stdev <- varPCA(Sobject, assay="SCT") * 100
ncomp <- ncol(Sobject@reductions$pca)

# view the PCA threshold
gg <- ElbowPlot(Sobject, ndims=ncomp) & ylab("Variance explained (%)") &
  theme_custom()
ggsave(file=file.path(savedir, "Sobject_pca.svg"), device="svg", gg)
#choose limit
ncomp <- 25 #based on the PCA graph


# ---------------------------------------------------------------
# DIMENSIONALITY REDUCTION : UMAP
# ---------------------------------------------------------------

# construct the shared nearest neighbors graph (SNN)
n_neighbors <- 20 #how many neighbors to consider for each cell
Sobject <- FindNeighbors(
  Sobject, dims = 1:ncomp, k.param = n_neighbors, annoy.metric = "cosine",
  n.trees = 200, prune.SNN = 1/15
)

# apply Leiden clustering algorithm for several resolutions
res_vec <- seq(2.5, 4, 0.5) #(lower resolution, higher resolution, step)
Sobject <- FindClusters(
  Sobject, algorithm = 4, resolution = res_vec, n.iter = 100,
  random.seed = 666
)

# construct UMAP embedding with custom parameters
Sobject <- RunUMAP(
  Sobject , dims = 1:ncomp, n.neighbors = n_neighbors, metric="cosine",
  reduction = "pca", umap.method = "umap-learn", n.epochs = 1000,
  learning.rate = 1, negative.sample.rate = 9, local.connectivity = 1,
  set.op.mix.ratio = 0.5, min.dist = 0.05, seed.use = 1672071369
)

cmap <- "tol.rainbow"

# export UMAP embedding with custom coloring based on clustering
res_vec <- colnames(Sobject@meta.data) %>% 
  strsplit(., split="SCT_snn_res.") %>%
  sapply(., function(x) x[2]) %>% 
  na.omit(.) %>%
  as.numeric(.) 
export_embeddings(Sobject, res_vec, cmap, savedir, assay="SCT", device = "svg")


# ---------------------------------------------------------------
# SAVE Seurat object
# ---------------------------------------------------------------

SaveH5Seurat(Sobject, filename = file.path(savedir, "Kurmangaliyev_Subset.h5Seurat"))

## End ##


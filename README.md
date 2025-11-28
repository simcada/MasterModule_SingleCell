# Analysis of single-cell dataset

Data from: Kurmangaliyev, et al. (2020), "Transcriptional Programs Of Circuit Assembly In The Drosophila Visual System", Neuron 108, 1045-1057.E6. doi:10.1016/J.Neuron.2020.10.006

Link to download data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156455

---

## Installation

### Prerequisites/Requirements

This project requires to install :
- **R** and **RStudio** (**>4.2.1**),
- a "recent" version of Conda (**> 23.1.0**).

Conda is used to build the self-contained Python environment that hosts the following packages :
- **umap-learn** version **0.5.3**: used for building the umap embedding (instead of the default R package uwot),
- **leidenalg** version **0.9.0 (or later)**: used for applying Leiden clustering algorithm when calling the `Seurat::FindClusters` function.


### Steps

1. Download the repository: https://github.com/simcada/MasterModule_SingleCell.git

2. Open the Rproj file: "ClusteringAnalysis.Rproj"

**IMPORTANT**: (if you work on **Windows**: change "/" into "\ \") 

3. Install the R environment with the R command: 
  * `install.packages("renv")`
  * and then: `renv::init(".")`

Error can happened at this stage, if so => install manually the different packages.
For example with: `renv::install("reticulate")` 
And make a new snapshot of the environment with: `renv::snapshot()`

  * And install the python environment with the command: `renv:::renv_python_conda_restore(".")`

4. Exit RStudio

5. In the .Rprofile file, add on the first row: `Sys.setenv('RETICULATE_PYTHON_ENV'='YourPATH/clustering_analysis-master/renv/python/condaenvs/renv-python')`
(if you work on Windows: change / into \\) 
6. Open Rstudio

7. `reticulate::py_config()`

Error can happened at this stage, if so => install manually the different packages.
For example with: 
  `install.packages("Matrix", repos = "http://R-Forge.R-project.org")` 
  `install.packages('Seurat')`
  `install.packages('remotes')`
  `remotes::install_github("mojaveazure/seurat-disk")`
  `install.packages('janitor')` 
And make a new snapshot of the environment with: `renv::snapshot()`
And restart your R session

Note: if installation with R doesn't work, try installing via conda => example: in your terminal window: `conda install -c conda-forge r-seuratdisk`


---

## Usage

To perform the analysis:

1. In the "share" folder: create a folder named "data",

2. In the new "data"folder: create a folder named "Kurmangaliyev",

3. Copy the dataset from the paper in the "Kurmangaliyev",

4. In the "clustering_analysis-master" folder: create a folder named "results" (where the results figures and files will be stored),

5. Run the "CreateSeuratObject.R" script with the following specification at the beginning:
   - the variable `datadir` should correctly point to the subfolder "data" created at step 2,
   - the variable `savedir` should correctly point to the subfolder "results" created at step 1,
   - the variable `metadata` should correctly point to the subfolder "Kurmangaliyev" created at step 2,

6. Run the "analysis_Subset.R" script (we will work with the Seurat object create with step 3) with the same specification that step 3




## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)


## ----eval=FALSE---------------------------------------------------------------
#  install.packages('harmony')

## -----------------------------------------------------------------------------
## Source required data
data("pbmc_stim")
pbmc <- CreateSeuratObject(counts = cbind(pbmc.stim, pbmc.ctrl), project = "PBMC", min.cells = 5)

## Separate conditions

pbmc@meta.data$stim <- c(rep("STIM", ncol(pbmc.stim)), rep("CTRL", ncol(pbmc.ctrl)))

## ----eval = FALSE, class.source='fold-hide'-----------------------------------
#  library(Matrix)
#  ## Download and extract files from GEO
#  ##setwd("/path/to/downloaded/files")
#  genes =  read.table("GSE96583_batch2.genes.tsv.gz", header = FALSE, sep = "\t")
#  
#  pbmc.ctrl.full = as.readMM("GSM2560248_2.1.mtx.gz")
#  colnames(pbmc.ctrl.full) = paste0(read.table("GSM2560248_barcodes.tsv.gz", header = FALSE, sep = "\t")[,1], "-1")
#  rownames(pbmc.ctrl.full) = genes$V1
#  
#  pbmc.stim.full = readMM("GSM2560249_2.2.mtx.gz")
#  colnames(pbmc.stim.full) = paste0(read.table("GSM2560249_barcodes.tsv.gz", header = FALSE, sep = "\t")[,1], "-2")
#  rownames(pbmc.stim.full) = genes$V1
#  
#  library(Seurat)
#  
#  pbmc <- CreateSeuratObject(counts = cbind(pbmc.stim.full, pbmc.ctrl.full), project = "PBMC", min.cells = 5)
#  pbmc@meta.data$stim <- c(rep("STIM", ncol(pbmc.stim.full)), rep("CTRL", ncol(pbmc.ctrl.full)))
#  
#  
#  
#  
#  # Running Harmony
#  
#  Harmony works on an existing matrix with cell embeddings and outputs its transformed version with the datasets aligned according to some user-defined experimental conditions. By default, harmony will look up the `pca` cell embeddings and use these to run harmony. Therefore, it assumes that the Seurat object has these embeddings already precomputed.
#  
#  ## Calculate PCA cell embeddings
#  
#  Here, using `Seurat::NormalizeData()`, we will be generating a union of highly variable genes using each condition (the control and stimulated cells). These features are going to be subsequently used to generate the 20 PCs with `Seurat::RunPCA()`.
#  

## -----------------------------------------------------------------------------
pbmc <- pbmc %>%
    NormalizeData(verbose = FALSE)

VariableFeatures(pbmc) <- split(row.names(pbmc@meta.data), pbmc@meta.data$stim) %>% lapply(function(cells_use) {
    pbmc[,cells_use] %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
        VariableFeatures()
}) %>% unlist %>% unique

pbmc <- pbmc %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(features = VariableFeatures(pbmc), npcs = 20, verbose = FALSE)

## ---- eval=FALSE--------------------------------------------------------------
#  ## run harmony with default parameters
#  pbmc <- pbmc %>% RunHarmony("stim")
#  ## is equivalent to:
#  pbmc <- RunHarmony(pbmc, "stim")

## ---- fig.width = 4, fig.height = 3, fig.align = "center", out.width="50%", fig.cap="By setting `plot_converge=TRUE`, harmony will generate a plot with its objective showing the flow of the integration. Each point represents the cost measured after a clustering round. Different colors represent different Harmony iterations which is controlled by `max_iter` (assuming that early_stop=FALSE). Here `max_iter=10` and up to 10 correction steps are expected. However, `early_stop=TRUE` so harmony will stop after the cost plateaus."----

pbmc <- pbmc %>% 
    RunHarmony("stim", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

## -----------------------------------------------------------------------------
harmony.embeddings <- Embeddings(pbmc, reduction = "harmony")

## ---- fig.width=7, fig.height=3, out.width="100%", fig.align="center", fig.cap="Evaluate harmonization of stim parameter in the harmony generated cell embeddings"----

p1 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "stim")
p2 <- VlnPlot(object = pbmc, features = "harmony_1", group.by = "stim",  pt.size = .1)
plot_grid(p1,p2)

## ---- fig.width = 6, fig.height=3, out.width="100%"---------------------------

DimHeatmap(object = pbmc, reduction = "harmony", cells = 500, dims = 1:3)

## -----------------------------------------------------------------------------
pbmc <- pbmc %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = 0.5) 

## ---- fig.width=5, fig.height=2.5, fig.align="center", fig.cap="t-SNE Visualization of harmony embeddings"----
pbmc <- pbmc %>%
    RunTSNE(reduction = "harmony")


p1 <- DimPlot(pbmc, reduction = "tsne", group.by = "stim", pt.size = .1)
p2 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = .1)
plot_grid(p1, p2)


## ---- fig.width = 7, fig.height = 7, out.width="100%", fig.cap="Expression of gene panel heatmap in the harmonized PBMC dataset"----
FeaturePlot(object = pbmc, features= c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), 
            min.cutoff = "q9", cols = c("lightgrey", "blue"), pt.size = 0.5)


## ---- fig.width=5, fig.height=2.5, fig.align="center", fig.cap="UMAP Visualization of harmony embeddings"----
pbmc <- pbmc %>%
    RunUMAP(reduction = "harmony",  dims = 1:20)

p1 <- DimPlot(pbmc, reduction = "umap", group.by = "stim", pt.size = .1)
p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE,  pt.size = .1)
plot_grid(p1, p2)


## -----------------------------------------------------------------------------
sessionInfo()


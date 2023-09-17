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


## -----------------------------------------------------------------------------
## Install latest branch of harmony
## devtools::install_github('immunogenomics/harmony', force = TRUE)

## -----------------------------------------------------------------------------
## Source required data
data("pbmc_stim")
pbmc <- CreateSeuratObject(counts = cbind(pbmc.stim, pbmc.ctrl), project = "PBMC", min.cells = 5)

## Separate conditions

pbmc@meta.data$stim <- c(rep("STIM", ncol(pbmc.stim)), rep("CTRL", ncol(pbmc.ctrl)))

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

## ---- echo=FALSE--------------------------------------------------------------
## run harmony with default parameters
pbmc <- pbmc %>% RunHarmony("stim")
## is equivalent to:
pbmc <- RunHarmony(pbmc, "stim")

## ---- fig.width = 8, fig.height = 6, out.width="100%", fig.cap="By setting `plot_converge=TRUE`, harmony will generate a plot with its objective showing the flow of the integration. Each point represents the cost measured after a clustering round. Different colors represent different Harmony iterations which is controlled by `max_iter` (assuming that early_stop=FALSE). Here `max_iter=10` and up to 10 correction steps are expected. However, `early_stop=TRUE` so harmony will stop after the cost plateaus."----

pbmc <- pbmc %>% 
    RunHarmony("stim", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

## -----------------------------------------------------------------------------
harmony.embeddings <- Embeddings(pbmc, reduction = "harmony")

## ---- fig.width = 8, fig.height = 4, out.width="100%", fig.cap="Evaluate harmonization of stim parameter in the harmony generated cell embeddings"----

p1 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "stim")
p2 <- VlnPlot(object = pbmc, features = "harmony_1", group.by = "stim",  pt.size = .1)
plot_grid(p1,p2)

## ---- fig.width = 6, fig.height=6, out.width="100%"---------------------------

DimHeatmap(object = pbmc, reduction = "harmony", cells = 500, dims = 1:6)

## -----------------------------------------------------------------------------
pbmc <- pbmc %>%
    FindNeighbors(reduction = "harmony") %>%
    FindClusters(resolution = 0.5) 

## ---- fig.width = 8, fig.height = 4, out.width = "100%", fig.cap="t-SNE Visualization of harmony embeddings"----
pbmc <- pbmc %>%
    RunTSNE(reduction = "harmony")


p1 <- DimPlot(pbmc, reduction = "tsne", group.by = "stim", pt.size = .1)
p2 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = .1)
plot_grid(p1, p2)


## ---- fig.width = 8, fig.height = 10, out.width="100%", fig.cap="Expression of gene panel heatmap in the harmonized PBMC dataset"----
FeaturePlot(object = pbmc, features= c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), 
            min.cutoff = "q9", cols = c("lightgrey", "blue"), pt.size = 0.5)


## ---- fig.width = 8, fig.height=4, out.width="100%", fig.cap="UMAP Visualization of harmony embeddings"----
pbmc <- pbmc %>%
    RunUMAP(reduction = "harmony",  dims = 1:20)

p1 <- DimPlot(pbmc, reduction = "umap", group.by = "stim", pt.size = .1)
p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE,  pt.size = .1)
plot_grid(p1, p2)


## -----------------------------------------------------------------------------
sessionInfo()


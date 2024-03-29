Harmony <img src="man/figures/logo.png" width="181px" align="right" />
===========

[![Travis-CI Build Status](https://travis-ci.org/immunogenomics/harmony.svg?branch=master)](https://travis-ci.org/immunogenomics/harmony)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/immunogenomics/harmony?branch=master&svg=true)](https://ci.appveyor.com/project/immunogenomics/harmony)
[![DOI](https://zenodo.org/badge/doi/10.1038/s41592-019-0619-0.svg)](https://doi.org/10.1038/s41592-019-0619-0)

*Fast, sensitive and accurate integration of single-cell data with Harmony*

Check out the manuscript in Nature Methods: 
- [nature website](https://www.nature.com/articles/s41592-019-0619-0)
- [read link](https://www.nature.com/articles/s41592-019-0619-0.epdf?shared_access_token=rDg_Rd07lrFXExt_ySj7V9RgN0jAjWel9jnR3ZoTv0NfDJkKCfDV_X9Mq3lweQmKiXEXxhrebQRjJEZdc-xNv6-7ZN1XotlD_mo5TSS4Z4eWn-kUo6mBwA5dEAKlTfR8OT6E10MZY_E-906ajbzvgg%3D%3D)

For Python users, check out the [harmonypy package](https://github.com/slowkow/harmonypy) by Kamil Slowikowski. 

# System requirements 

Harmony has been tested on R versions >= 3.4. Please consult the DESCRIPTION file for more details on required R packages. Harmony has been tested on Linux, OS X, and Windows platforms.

# Installation

To run Harmony, open R and install harmony from CRAN: 

```r
install.packages("harmony")
```

If you'd like the latest development version, install from this github directly: 

```r
devtools::install_github("immunogenomics/harmony", build_vignettes=TRUE)
```


# Usage

 Harmony is designed to be user-friendly and supports some SingleCellExperiment and Seurat R analysis pipelines. Alternatively, it can be used in standalone mode.

## Quick Start 

### Standalone Mode
Check out this [vignette](http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/quickstart.html) for a quick start tutorial which demonstrates the usage of the tool in standalone mode.

At minimum the following parameters need to be specified to achieve an integration. 

```r
library(harmony)
my_harmony_embeddings <- RunHarmony(my_pca_embeddings, meta_data, "dataset")
```


## Seurat Objects

By default, the harmony API works on Seurats PCA cell embeddings and corrects them. You can run Harmony within your Seurat workflow with `RunHarmony()`. Prior `RunHarmony()` the PCA cell embeddings need to be precomputed through Seurat's API. For downstream analyses, use the `harmony` embeddings instead of `pca`.

For example, the following snippet run Harmony and then calculates UMAP of the corrected input embeddings:

```r
seuratObj <- RunHarmony(seuratObj, "dataset")
seuratObj <- RunUMAP(seuratObj, reduction = "harmony")
```

For a more detailed overview of the `RunHarmony()` Seurat interface check, the [Seurat vignette](http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/Seurat.html)

## Harmony with two or more covariates

Harmony can integrate over multiple covariates. To do this, specify a vector covariates to integrate. 

```r
my_harmony_embeddings <- RunHarmony(
  my_pca_embeddings, meta_data, c("dataset", "donor", "batch_id")
)
```

Do the same with your Seurat object: 

```r
seuratObject <- RunHarmony(seuratObject, c("dataset", "donor", "batch_id"))
```

## Advanced tutorial 

The examples above all return integrated PCA embeddings. We created a [detailed walkthrough](http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/detailedWalkthrough.html) that explores the internal data structures and mechanics of the Harmony algorithm.


# Performance Notes

## BLAS vs. OPENBLAS

R distributions can be bundled with different scientific computing libraries. This can drastically impact harmony's performance. Rstudio comes by default with BLAS. In contrast, conda distributions of R are bundled with OPENBLAS. Overall, our benchmarks show that **harmony+OPENBLAS is substantially faster compared harmony+BLAS**. Therefore users with large datasets will benefit using OPENBLAS.

## Multithreading in OPENBLAS

One caveat is that OPENBLAS uses OPENMP to parallelize operations. By default, OPENBLAS will utilize all cores for these operations. While in theory this accelerates runtimes, in practice harmony is not optimized for multi-threaded performance and the unoptimized parallelization granularity may result in significantly slower run times and inefficient resource utilization (wasted CPU cycles). Therefore, by default harmony turns off multi-threading. However, very large datasets >1M may benefit from parallelization. This behavior can be controlled by the `ncores` parameter which expects a number threads which harmony will use for its math operation. Users are advised to increase gradually `ncores` and assess potential performance benefits.


# Reproducing results from manuscript

Code to reproduce Harmony results from the Korsunsky et al 2019 manuscript will be made available on github.com/immunogenomics/harmony2019. 



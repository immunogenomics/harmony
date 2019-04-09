#' List of expression matrix, metadata table, and PCA embeddings for 500 stimulated and 500 unstimulated PBMCs. 
#'
#' @format List of 3 objects: 
#'   1) exprs_norm_vargenes: dgCMatrix of 1636 genes by 1000 cells. 
#'   2) meta_data: data.frame
#'   3) pca_embeddings: matrix of 1000 cells by 20 PCs
#'
#' @source \url{https://www.nature.com/articles/nbt.4042}
"pbmc.small"

#' Seurat object with the same data as pbmc.small. 
#'
#' @format Seurat v2.3.4 object. 
#'
#' @source \url{https://www.nature.com/articles/nbt.4042}
"pbmc.small.seurat"

#' List of metadata table and scaled PCs matrix
#' 
#' @format:
#'   meta_data: data.table of 9478 rows with information about dataset and cell_type
#'   scaled_pcs: data.table of 9478 rows (cells) and 20 columns (PCs)
#' 
#' 
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"celllines"
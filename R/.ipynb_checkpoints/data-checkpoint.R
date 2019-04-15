#' List of metadata table and scaled PCs matrix
#' 
#' @format:
#'   meta_data: data.table of 9478 rows with information about dataset and cell_type
#'   scaled_pcs: data.table of 9478 rows (cells) and 20 columns (PCs)
#' 
#' 
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"cell_lines"

#' Same as cell_lines but smaller (300 cells).
#' 
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"cell_lines_small"

#' Same as cell_lines_small but as Seurat object. Expression matrices filled in with dummy values. 
#' 
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets}
"cell_lines_small_seurat"


#' List of metadata table and scaled PCs matrix
#' 
#' @format:
#'   meta_data: data.table of 9478 rows with defining dataset and cell_type
#'   scaled_pcs: data.table of 9478 rows (cells) and 20 columns (PCs)
#' 
#' @source \url{https://www.10xgenomics.com}
"cell_lines"

#' Same as cell_lines but smaller (300 cells).
#' 
#' @source \url{https://www.10xgenomics.com}
"cell_lines_small"


#' Gene expression data of control PBMC from Kang et al. 2017. This
#' contains a sample of 1000 cells from that condition and is used for
#' the Seurat Vignette.
#' 
#' @source \url{https://doi.org/10.1038/nbt.4042}
"pbmc.ctrl"


#' Gene expression data of stimulated PBMC from Kang et al. 2017. This
#' contains a sample of 1000 cells from that condition and is used for
#' the Seurat Vignette.
#' 
#' @source \url{https://doi.org/10.1038/nbt.4042}
"pbmc.stim"



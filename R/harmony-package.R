

#' @useDynLib harmony
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp loadModule
#' @importFrom methods new
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom cowplot plot_grid
#' @importFrom rlang .data
#' @importFrom rlang `%||%`
#' @importFrom stats model.matrix
loadModule("harmony_module", TRUE)
NULL

#' Harmony: fast, accurate, and robust single cell integration.
#'
#' Algorithm for single cell integration.
#'
#' @section Usage:
#'
#' 
#' ?RunHarmony to run Harmony on cell embeddings matrix, Seurat or
#' SingleCellExperiment objects.
#' 
#' @section Useful links:
#'
#' \enumerate{
#' \item Report bugs at \url{https://github.com/immunogenomics/harmony/issues}
#' \item Read the manuscript
#' \doi{10.64898/2026.03.16.711825}
#' }
#' 
#' @name harmony
#' @docType package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

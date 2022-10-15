#' Harmony: fast, accurate, and robust single cell integration.
#'
#' Algorithm for single cell integration.
#'
#' @section Usage:
#'
#' \enumerate{
#' \item ?HarmonyMatrix to run Harmony on gene expression or PCA
#' embeddings matrix.
#' \item ?RunHarmony to run Harmony on Seurat or SingleCellExperiment objects.
#' }
#' @section Useful links:
#'
#' \enumerate{
#' \item Report bugs at \url{https://github.com/immunogenomics/harmony/issues}
#' \item Read the manuscript
#' \href{https://www.biorxiv.org/content/10.1101/461954v2}{online}.
#' }
#'
#'
#' @name harmony
#' @docType package
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

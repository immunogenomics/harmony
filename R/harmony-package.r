#' harmony.
#'
#' @name harmony
#' @docType package
#' @useDynLib harmony
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp loadModule
#' @importFrom methods new
#' @importFrom methods as
#' @importFrom rlang .data
loadModule("harmony_module", TRUE)
NULL

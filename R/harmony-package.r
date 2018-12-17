#' harmony.
#'
#' @name harmony
#' @docType package
#' @useDynLib harmony
#' @importFrom Rcpp sourceCpp loadModule
loadModule("harmony_module", TRUE)
NULL

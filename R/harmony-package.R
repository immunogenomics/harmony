#' harmony.
#'
#' @name harmony
#' @docType package
#' @useDynLib harmony, .registration = TRUE
#' @importFrom Rcpp sourceCpp loadModule evalCpp
loadModule("harmony_module", TRUE)
NULL

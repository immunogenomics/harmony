#' Set advanced options for RunHarmony
#' @param alpha When setting lambda = NULL and use lambda estimation mode, 
#'     lambda would be determined by the expected number of cells assuming 
#'     idependece between batches and clusters. i.e., lambda = alpha * expected
#'     number of cells, default 0.2 and alpha should be 0 < alpha < 1
#' @param tau Protection against overclustering small datasets with 
#'     large ones. `tau` is the expected number of cells per cluster.
#' @param block.size What proportion of cells to update during clustering. 
#'     Between 0 to 1, default 0.05. Larger values may be faster but less 
#'     accurate.
#' @param max.iter.cluster Maximum number of rounds to run clustering 
#'     at each round of Harmony.
#' @param epsilon.cluster Convergence tolerance for clustering round 
#'     of Harmony. Set to -Inf to never stop early.
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to
#'     never stop early. When `epsilon.harmony` is set to not NULL, then
#'     user-supplied values of `early_stop` is ignored.
#' @returns Return a list for `.options` argument of `RunHarmony`
#' @export
#' @examples
#' ## If want to set lambda to be fixed to 1, do
#' \dontrun{
#' RunHarmony(data_meta, meta_data, vars_use,
#'               .options = harmony_options(lambda = c(1, 1)))
#' }
#' 
harmony_options <- function(
  alpha = 0.2,
  tau = 0,
  block.size = 0.05,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-3,
  epsilon.harmony = 1e-2) {
    
    block.size <- validate_block.size(block.size)
    
    out <- list(
        alpha = alpha,
        tau = tau,
        block.size = block.size,
        max.iter.cluster = max.iter.cluster,
        epsilon.cluster = epsilon.cluster,
        epsilon.harmony = epsilon.harmony
    )
    out <- structure(out, class = "harmony_options")
    return(out)
}

## Validate functions -----------------------------------------------------------
validate_block.size <- function(block.size) {
    if(block.size <= 0 | block.size > 1){
        stop('Error: block.size should be set between 0 and 1 (0 < block.size <= 1)')
    }
    return(block.size)
}


#' @importFrom methods hasArg
check_legacy_args <- function(...) {
    if (hasArg("do_pca") || hasArg("npcs")) legacy_warning("do_pca_npcs")
    if (hasArg("tau")) legacy_warning("tau")
    if (hasArg("block.size")) legacy_warning("block.size")
    if (hasArg("max.iter.harmony")) legacy_warning("max.iter.harmony")
    if (hasArg("max.iter.cluster")) legacy_warning("max.iter.cluster")
    if (hasArg("epsilon.cluster")) legacy_warning("epsilon.cluster")
    if (hasArg("epsilon.harmony")) legacy_warning("epsilon.harmony")
    
}




legacy_warning <- function(param) {
    common_warn <- paste0(
        "Warning: The parameter ", param, " is deprecated. ",
        "It will be ignored for this function call ",
        "and please remove parameter ", param, " in future function calls. ",
        "Advanced users can set value of parameter ", param,
        " by using parameter .options and function harmony_options()."
    )
    do_pca_npcs_warn <- paste0(
        "Warning: The parameters ", "do_pca and npcs", " are deprecated. ",
        "They will be ignored for this function call ",
        "and please remove parameters ", "do_pca and npcs",
        " and pass to harmony cell_embeddings directly."
    )
    max.iter.harmony_warn <- paste0(
        "Warning: The parameter ", "max.iter.harmony ",
        "is replaced with parameter ", "max_iter. ",
        "It will be ignored for this function call ",
        "and please use parameter ", "max_iter ", "in future function calls."
    )
    epsilon.harmony_warn <- paste0(
        "Warning: The parameter ", "epsilon.harmony", " is deprecated. ",
        "It will be ignored for this function call ",
        "and please remove parameter ", "epsilon.harmony",
        " in future function calls. ",
        "If users want to control if harmony would stop early or not, ",
        "use parameter ", "early_stop. ",
        "Advanced users can set value of parameter ", "epsilon.harmony",
        " by using parameter .options and function harmony_options()."
    )


    if (param %in% c("tau", "block.size", "max.iter.cluster",
                     "epsilon.cluster")) {
        warn_str <- common_warn
    }
    if (param == "do_pca_npcs") {
        warn_str <- do_pca_npcs_warn
    }
    if (param == "max.iter.harmony") {
        warn_str <- max.iter.harmony_warn
    }
    if (param == "epsilon.harmony") {
        warn_str <- epsilon.harmony_warn
    }

    rlang::warn(warn_str, .frequency = "once", .frequency_id = param)
}

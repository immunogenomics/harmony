#' Set advanced parameters for RunHarmony
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
#' @param batch.prop.cutoff During the integration step, if a batch
#'     has less of the specified proportion in a harmony cluster it
#'     will be excluded from the integration step. For example,
#'     batch.prop.cutoff=0.01 and a batch has less than 1/100 of its
#'     cells soft-assigned to a cluster this batch won't participating
#'     in the correction step for the particular batch.
#' @returns Return a list for `.options` argument of `RunHarmony`
#' @export
#' @examples
#' ## If want to set max.iter.cluster to be 100, do
#' \dontrun{
#' RunHarmony(data_meta, meta_data, vars_use,
#'               .options = harmony_options(max.iter.cluster = 100))
#' }
#' 
harmony_options <- function(
  alpha = 0.2,
  tau = 0,
  block.size = 0.05,
  max.iter.cluster = 4,
  epsilon.cluster = 1e-3,
  epsilon.harmony = 1e-2,
  batch.prop.cutoff = 1e-5) {
    
    block.size <- validate_block.size(block.size)
    
    out <- list(
        alpha = alpha,
        tau = tau,
        block.size = block.size,
        max.iter.cluster = max.iter.cluster,
        epsilon.cluster = epsilon.cluster,
        epsilon.harmony = epsilon.harmony,
        batch.prop.cutoff = batch.prop.cutoff
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
    all.args = list(...)    
    legarg <- c("do_pca", "npcs", "tau", "block.size",
                "max.iter.harmony", "max.iter.cluster",
                "epsilon.cluster", "epsilon.harmony")
    for(arg in names(all.args)) {
        if (arg %in% legarg) {
            legacy_error(arg)
            all.args[[arg]] = NULL
        }
    }
    if(length(all.args) > 0){
        stop(paste("Argument", names(all.args),"is unhandled. Please refer to the documentation for the valid harmony options!\n"))
    }
}




legacy_error <- function(param) {
    common_warn <- paste0(
        "Error: The parameter ", param, " has been dropped from the RunHarmony API. ",
        "It will be ignored for this function call ",
        "and please remove parameter ", param, " in future function calls. ",
        "Advanced users can set value of parameter ", param,
        " by using parameter .options and function harmony_options()."
    )
    do_pca_npcs_warn <- paste0(
        "Error: The parameters ", "do_pca and npcs", " have been dropped from the RunHarmony API. ",
        "They will be ignored for this function call ",
        "and please remove parameters ", "do_pca and npcs",
        " and pass to harmony cell_embeddings directly."
    )
    max.iter.harmony_warn <- paste0(
        "Error: The parameter ", "max.iter.harmony ",
        "is replaced with parameter ", "max_iter. ",
        "It will be ignored for this function call ",
        "and please use parameter ", "max_iter ", "in future function calls."
    )
    epsilon.harmony_warn <- paste0(
        "Error: The parameter ", "epsilon.harmony", " has been dropped from the RunHarmony API. ",
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
    if (param %in% c("do_pca", "npcs")) {
        warn_str <- do_pca_npcs_warn
    }
    if (param == "max.iter.harmony") {
        warn_str <- max.iter.harmony_warn
    }
    if (param == "epsilon.harmony") {
        warn_str <- epsilon.harmony_warn
    }
    stop(warn_str)
}

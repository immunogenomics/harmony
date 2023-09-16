#' Set advanced options for RunHarmony
#' @param lambda_range Default lambda_range = c(0.1, 10). Lambda is ridge 
#'     regression penalty parameter and smaller values result in more 
#'     aggressive correction. During harmony iterations, the appropriate value 
#'     of lambda is dynamically estimated. And parameter `lambda_range` set the 
#'     allowed range for lambda estimation. e.g. `lambda_range` = c(0.1, 10) 
#'     means that lambda can only vary between 0.1 and 10 when being 
#'     dynamically estimated. Note that when setting the upper and lower bound 
#'     of lambda_range to the same value would result in using a fixed lambda 
#'     throughout harmony iterations. e.g. `lambda_range` = c(1,1) would make 
#'     harmony using a fixed lambda = 1.
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
  lambda_range = c(0.1, 10),
  tau = 0,
  block.size = 0.05,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4) {
    
    lambda_range <- validate_lambda_range(lambda_range)
    block.size <- validate_block.size(block.size)
    
    out <- list(
        lambda_range = lambda_range,
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
validate_lambda_range <- function(lambda_range) {
    if (length(lambda_range) != 2) {
        stop('Error: lambda_range should have length == 2')
    }
    if (lambda_range[2] < lambda_range[1]) {
        stop('Error: lambda_range[2] cannot be smaller than lambda_range[1]')
    }
    if (lambda_range[1] <= 0) {
        stop('Error: lambda_range cannot be smaller or equal to 0')
    }
    
    if (lambda_range[2] == lambda_range[1]) {
        message("Automatic lambda estimation for ridge is disabled. Harmony will have a fixed lambda for all batches")
    }
    return(lambda_range)
}


validate_block.size <- function(block.size) {
    if(block.size <= 0 | block.size > 1){
        stop('Error: block.size should be set between 0 and 1 (0 < block.size <= 1)')
    }
    return(block.size)
}


#' @importFrom methods hasArg
check_legacy_args <- function(...) {
    if (hasArg(do_pca) || hasArg(npcs)) legacy_warning("do_pca_npcs")
    if (hasArg(tau)) legacy_warning("tau")
    if (hasArg(block.size)) legacy_warning("block.size")
    if (hasArg(max.iter.harmony)) legacy_warning("max.iter.harmony")
    if (hasArg(max.iter.cluster)) legacy_warning("max.iter.cluster")
    if (hasArg(epsilon.cluster)) legacy_warning("epsilon.cluster")
    if (hasArg(epsilon.harmony)) legacy_warning("epsilon.harmony")
    
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

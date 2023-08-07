#' Set advanced options for HarmonyMatrix
#' @param lambda Default lambda = c(0.1, 10). In the new release of Harmony,
#'     the appropriate value of lambda is dynamically estimated during harmony
#'     iterations. And parameter lambda set the allowed range for lambda
#'     estimation. e.g. lambda = c(0.1, 10) means that lambda can only vary
#'     between 0.1 and 10.
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
#' @returns Return a list for `.options` argument of `HarmonyMatrix`
#' @export
#' @examples
#' ## If want to set lambda to be fixed to 1, do
#' HarmonyMatrix(data_meta, meta_data, vars_use,
#'               .options = harmony_options(lambda = c(1, 1)))


harmony_options <- function(lambda = c(0.1, 10),
                            tau = 0,
                            block.size = 0.05,
                            max.iter.cluster = 20,
                            epsilon.cluster = 1e-5,
                            epsilon.harmony = NULL){
    lambda <- validate_lambda(lambda)
    tau <- validate_tau(tau)
    block.size <- validate_block.size(block.size)
    max.iter.cluster <- validate_max.iter.cluster(max.iter.cluster)
    epsilon.cluster <- validate_epsilon.cluster(epsilon.cluster)
    epsilon.harmony <- validate_epsilon.harmony(epsilon.harmony)
    
    out <- list(
        lambda = lambda,
        tau = tau,
        block.size = block.size,
        max.iter.cluster = max.iter.cluster,
        epsilon.cluster = epsilon.cluster,
        epsilon.harmony = epsilon.harmony
    )
    out <- structure(out, class = "harmony_options")
    return(out)
}

# Validate functions -----------------------------------------------------------

validate_lambda <- function(lambda){
    if(length(lambda) != 2){
        stop('Error: lambda should have length == 2')
    }
    if (lambda[2] < lambda[1]){
        stop('Error: lambda[2] cannot be smaller than lambda[1]')
    }
    if (lambda[1] <= 0){
        stop('Error: lambda cannot be smaller or equal to 0')
    }
    
    if (lambda[2] == lambda[1]) {
        message("Automatic lambda estimation for ridge is disabled. Harmony will have a fixed lambda for all batches")
    }
    return(lambda)
}

validate_tau <- function(tau){
    return(tau)
}


validate_block.size <- function(block.size){
    if(block.size <= 0 | block.size > 1){
        stop('Error: block.size should be set between 0 and 1 (0 < block.size <= 1)')
    }
    return(block.size)
}

validate_max.iter.cluster <- function(max.iter.cluster){
    if(max.iter.cluster != round(max.iter.cluster)){

    }
    # integer, [1 inf]
    return(max.iter.cluster)
}


validate_epsilon.cluster <- function(epsilon.cluster){
    return(epsilon.cluster)
}

validate_epsilon.harmony <- function(epsilon.harmony){
    return(epsilon.harmony)
}
#' Set advanced options for HarmonyMatrix
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
#' @returns Return a list for `.options` argument of `HarmonyMatrix`
#' @export
#' @examples
#' ## If want to set lambda to be fixed to 1, do
#' HarmonyMatrix(data_meta, meta_data, vars_use,
#'               .options = harmony_options(lambda = c(1, 1)))


harmony_options <- function(lambda_range = c(0.1, 10),
                            tau = 0,
                            block.size = 0.05,
                            max.iter.cluster = 20,
                            epsilon.cluster = 1e-5,
                            epsilon.harmony = NULL) {
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

# Validate functions -----------------------------------------------------------

validate_lambda_range <- function(lambda_range) {
    if (length(lambda_range) != 2){
        stop('Error: lambda_range should have length == 2')
    }
    if (lambda_range[2] < lambda_range[1]){
        stop('Error: lambda_range[2] cannot be smaller than lambda_range[1]')
    }
    if (lambda_range[1] <= 0){
        stop('Error: lambda_range cannot be smaller or equal to 0')
    }
    
    if (lambda_range[2] == lambda_range[1]) {
        message("Automatic lambda estimation for ridge is disabled. Harmony will have a fixed lambda for all batches")
    }
    return(lambda_range)
}


validate_block.size <- function(block.size){
    if(block.size <= 0 | block.size > 1){
        stop('Error: block.size should be set between 0 and 1 (0 < block.size <= 1)')
    }
    return(block.size)
}

#' Main Harmony interface
#' 
#' Use this  to run the Harmony algorithm directly on cell embedding
#' matrix. 
#' 
#' @param data_mat Matrix of cell embeddings. Cells can be rows or
#'     columns and will be inferred by the rows of meta_data.
#' @param meta_data Either (1) Dataframe with variables to integrate
#'     or (2) vector with labels.
#' @param vars_use If meta_data is dataframe, this defined which
#'     variable(s) to remove (character vector).
#' @param theta Diversity clustering penalty parameter. Specify for
#'     each variable in vars_use Default theta=2. theta=0 does not
#'     encourage any diversity. Larger values of theta result in more
#'     diverse clusters.
#' @param sigma Width of soft kmeans clusters. Default
#'     sigma=0.1. Sigma scales the distance from a cell to cluster
#'     centroids. Larger values of sigma result in cells assigned to
#'     more clusters. Smaller values of sigma make soft kmeans cluster
#'     approach hard clustering.
#' @param lambda Ridge regression penalty. Default lambda=1. Bigger
#'     values protect against over correction. If several covariates
#'     are specified, then lambda can also be a vector which needs to
#'     be equal length with the number of variables to be
#'     corrected. In this scenario, each covariate level group will be
#'     assigned the scalars specified by the user. If set to NULL,
#'     harmony will determine lambdas automatically and try to
#'     minimize overcorrection (beta).
#' @param nclust Number of clusters in model. nclust=1 equivalent to
#'     simple linear regression.
#' @param max_iter Maximum number of rounds to run Harmony. One round
#'     of Harmony involves one clustering and one correction step.
#' @param early_stop Enable early stopping for harmony. The
#'     harmonization process will stop when the change of objective
#'     function between corrections drops below 1e-4
#' @param ncores Number of processors to be used for math operations
#'     when optimized BLAS is available. If BLAS is not supporting
#'     multithreaded then this option has no effect. By default,
#'     ncore=1 which runs as a single-threaded process. Although
#'     Harmony supports multiple cores, it is not optimized for
#'     multithreading. Increase this number for large datasets iff
#'     single-core performance is not adequate.
#' @param plot_convergence Whether to print the convergence plot of
#'     the clustering objective function. TRUE to plot, FALSE to
#'     suppress. This can be useful for debugging.
#' @param return_object (Advanced Usage) Whether to return the Harmony
#'     object or only the corrected PCA embeddings.
#' @param verbose Whether to print progress messages. TRUE to print,
#'     FALSE to suppress.
#' @param .options Advanced parameters of RunHarmony. This must be the
#'     result from a call to `harmony_options`. See ?`harmony_options`
#'     for more details.
#' @param ... other parameters that are not part of the API
#' @return By default, matrix with corrected PCA embeddings. If
#'     return_object is TRUE, returns the full Harmony object (R6
#'     reference class type).
#'
#' @export 
#' 
#' @examples
#' 
#' 
#' ## By default, Harmony inputs a cell embedding matrix
#' \dontrun{
#' harmony_embeddings <- RunHarmony(cell_embeddings, meta_data, 'dataset')
#' }
#' 
#' ## If PCA is the input, the PCs need to be scaled
#' data(cell_lines_small)
#' pca_matrix <- cell_lines_small$scaled_pcs
#' meta_data <- cell_lines_small$meta_data
#' harmony_embeddings <- RunHarmony(pca_matrix, meta_data, 'dataset')
#' 
#' ## Output is a matrix of corrected PC embeddings
#' dim(harmony_embeddings)
#' harmony_embeddings[seq_len(5), seq_len(5)]
#' 
#' ## Finally, we can return an object with all the underlying data structures
#' harmony_object <- RunHarmony(pca_matrix, meta_data, 'dataset', return_object=TRUE)
#' dim(harmony_object$Y) ## cluster centroids
#' dim(harmony_object$R) ## soft cluster assignment
#' dim(harmony_object$Z_corr) ## corrected PCA embeddings
#' head(harmony_object$O) ## batch by cluster co-occurence matrix
#' 
RunHarmony.default <- function(
  data_mat,
  meta_data,
  vars_use,
  theta = NULL,
  sigma = 0.1,
  lambda = 1,
  nclust = NULL,
  max_iter = 10,
  early_stop = TRUE,
  ncores = 1,
  plot_convergence = FALSE,
  return_object = FALSE,
  verbose = TRUE,
  .options = harmony_options(),
  ...
  ) {
    
    
    tryCatch({
        ## The following block may fail in some build environments.
        ## We control the flow using the single.thread.mode
        single.thread.mode <- FALSE        
        ## Sanity check for number of cores
        max.cores <- RhpcBLASctl::omp_get_max_threads() 
        if ((ncores != as.integer(ncores)) || (ncores < 1) || (ncores > max.cores)) {
            stop(paste0(
                "Invalid number of ncores provided (", ncores, "). ",
                "Acceptable range of ncores: 1 -", max.cores))
        }
        prev.ncores.blas <- RhpcBLASctl::blas_get_num_procs()
        prev.ncores.omp <- RhpcBLASctl::omp_get_num_procs()
    },
    error = function() {
        warning(paste(
            "RhpcBLASctl was unable to set number of cores.",
            "Running in single-thread mode"
        ))
        single.thread.mode <- TRUE
    })
    
    
    
    tryCatch({
        ## Check legacy arguments
        check_legacy_args(...)

        ## Set threads if BLAS threas are set/detected properly
        if (!single.thread.mode) {
            RhpcBLASctl::blas_set_num_threads(ncores)
            RhpcBLASctl::omp_set_num_threads(ncores)
        }
        
        
        ## Parameter setting --------------------------------------------------------

        if (!inherits(.options, "harmony_options")) {
            stop("Error: .options must be created from harmony_options()!")
        }

        if (early_stop) {
            epsilon.harmony <- .options$epsilon.harmony
        } else {
            epsilon.harmony = -Inf
        }
        max.iter.harmony <- max_iter
        lambda_range <- .options$lambda_range
        tau <- .options$tau
        block.size <- .options$block.size
        max.iter.cluster <- .options$max.iter.cluster
        epsilon.cluster <- .options$epsilon.cluster   
        
        

        ## TODO: check for 
        ##    partially observed batch variables (WARNING)
        ##    batch variables with only 1 level (WARNING)
        ##    if lambda given, check correct length
        ##    if theta given, check correct length
        ##    very small batch size and tau=0: suggest tau>0
        ##    is PCA correct? 
        if (!(is(meta_data, 'data.frame') | is(meta_data, 'DataFrame'))) {
            if (length(meta_data) %in% dim(data_mat)) {
                meta_data <- data.frame(batch_variable = meta_data)
                vars_use <- 'batch_variable'
            } else {
                stop('meta_data must be either a data.frame or a vector with batch 
                values for each cell')
            }
        }
        
        if (is.null(vars_use) | any(!vars_use %in% colnames(meta_data))) {
            msg <- gettextf('must provide variables names (e.g. vars_use=%s)', 
                            sQuote('stim'))
            stop(msg)
        }

        ## Number of cells
        N <- nrow(meta_data)
        
        ## Check if we need to transpose our data
        if (nrow(data_mat) == N) {
            message("Transposing data matrix")
            data_mat <- Matrix::t(data_mat)
        }

        if (ncol(data_mat) != N) {
            stop("number of labels do not correspond to number of 
                samples in data matrix")
        }
        

                                        # determine K if null
        if (is.null(nclust)) {
            nclust <- min(round(N / 30), 100)
        }
        
                                        # determine theta if null
        if (is.null(theta)) {
            theta <- rep(2, length(vars_use))
        } else if (length(theta) != length(vars_use)) {
            stop('Please specify theta for each variable')
        }
        
                                        # determine sigma if it is a scalar
        if (length(sigma) == 1 & nclust > 1) {
            sigma <- rep(sigma, nclust)
        }
        

        ## Pre-compute some useful statistics
        phi <- Reduce(rbind, lapply(vars_use, function(var_use) {
            res <- Matrix::sparse.model.matrix(~0 + as.factor(meta_data[[var_use]]))
            Matrix::t(res)
        }))

        ## ## number of cells per batch
        N_b <- Matrix::rowSums(phi)

        ## Number of factors per covariate
        B_vec <- Reduce(c, lapply(vars_use, function(var_use) {
            nlevels(as.factor(meta_data[[var_use]]))
        }))

        if (is.null(lambda)) {
            lambda_vec <- -1 ## This enables automatic estimation in the backend
        } else if (!is.vector(lambda)) {
            lambda_vec <- c(0, rep(lambda, sum(B_vec)))
        } else if (length(lambda) != length(vars_use)) {
            stop(paste0("You specified a lambda value for each ",
                        "covariate but the number of lambdas specified (",
                        length(lambda), ") and the number of covariates (",
                        length(vars_use),") mismatch."))
            if(any(!(lambda > 0))) {
                stop("Provided lambdas must be positive")
            }
        } else {
            lambda_vec <- Reduce(c, lapply(seq_len(length(B_vec)), function(b)
                rep(lambda[b], B_vec[b])))
            lambda_vec <- c(0, lambda_vec)
        }
        
        ## Calculate theta (#covariates) x (#levels)
        theta <- Reduce(c, lapply(seq_len(length(B_vec)), function(b)
            rep(theta[b], B_vec[b])))

        ## Theta scaling
        theta <- theta * (1 - exp(-(N_b / (nclust * tau))^2))
        
        ## RUN HARMONY
        harmonyObj <- new(harmony)
        
        harmonyObj$setup(
                       data_mat, phi, sigma, theta, lambda_vec,
                       max.iter.cluster, epsilon.cluster,
                       epsilon.harmony, nclust, block.size,
                       lambda_range, B_vec, verbose
                   )
        
        harmonyObj$init_cluster_cpp(0)
        
        harmonize(harmonyObj, max.iter.harmony, verbose)
        
        if (plot_convergence) graphics::plot(HarmonyConvergencePlot(harmonyObj))

        
        ## Return either the R6 Harmony object or the corrected PCA matrix
        if (return_object) {
            return(harmonyObj)
        } else {
            res <- as.matrix(harmonyObj$Z_corr)
            row.names(res) <- row.names(data_mat)
            colnames(res) <- colnames(data_mat)
            return(t(res))
        }

    }, ## main tryCatch block ends here
    
    finally={
        if(!single.thread.mode) {
            RhpcBLASctl::blas_set_num_threads(prev.ncores.blas)
            RhpcBLASctl::omp_set_num_threads(prev.ncores.omp)
        }        
    })
    
    
    
}

#' @rdname RunHarmony
#' @export
HarmonyMatrix <- function(...) {
    .Deprecated("RunHarmony", msg="HarmonyMatrix is deprecated and will be removed in the future from the API in the future")
    RunHarmony(...)
}

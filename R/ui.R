#' Main Harmony interface
#' 
#' Use this to run the Harmony algorithm on gene expression or PCA matrix. 
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
#' @param nclust Number of clusters in model. nclust=1 equivalent to
#'     simple linear regression.
#' @param max_iter Maximum number of rounds to run Harmony. One round of
#'     Harmony involves one clustering and one correction step.
#' @param early_stop When to stop Harmony iteration before reaching max_iter 
#'     when the change in objectie function is small enough (< 1e-4)
#' @param plot_convergence Whether to print the convergence plot of
#'     the clustering objective function. TRUE to plot, FALSE to
#'     suppress. This can be useful for debugging.
#' @param return_object (Advanced Usage) Whether to return the Harmony
#'     object or only the corrected PCA embeddings.
#' @param verbose Whether to print progress messages. TRUE to print,
#'     FALSE to suppress.
#' @param reference_values (Advanced Usage) Defines reference
#'     dataset(s).  Cells that have batch variables values matching
#'     reference_values will not be moved.
#' @param .options Advanced parameters of HarmonyMatrix. This must be the 
#'     result from a call to `harmony_options`. See ?`harmony_options` for more
#'     details.
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
#' harmony_embeddings <- HarmonyMatrix(cell_embeddings, meta_data, 'dataset')
#' }
#' 
#' ## If PCA is the input, the PCs need to be scaled
#' data(cell_lines_small)
#' pca_matrix <- cell_lines_small$scaled_pcs
#' meta_data <- cell_lines_small$meta_data
#' harmony_embeddings <- HarmonyMatrix(pca_matrix, meta_data, 'dataset')
#' 
#' ## Output is a matrix of corrected PC embeddings
#' dim(harmony_embeddings)
#' harmony_embeddings[seq_len(5), seq_len(5)]
#' 
#' ## Finally, we can return an object with all the underlying data structures
#' harmony_object <- HarmonyMatrix(pca_matrix, meta_data, 'dataset', return_object=TRUE)
#' dim(harmony_object$Y) ## cluster centroids
#' dim(harmony_object$R) ## soft cluster assignment
#' dim(harmony_object$Z_corr) ## corrected PCA embeddings
#' head(harmony_object$O) ## batch by cluster co-occurence matrix
#' 
HarmonyMatrix <- function(
    data_mat, meta_data, vars_use, theta = NULL, sigma = 0.1, nclust = NULL,
    max_iter = 10, early_stop = TRUE, plot_convergence = FALSE,
    return_object = FALSE, verbose = TRUE, reference_values = NULL,
    .options = harmony_options(), ...
    ) {

    # Parameter checking -------------------------------------------------------
    if (hasArg(do_pca) || hasArg(npcs)) {
        rlang::warn(
            paste0("Warning: Function parameters do_pca and npcs is deprecated ",
                   "in this new release of harmony. They will be ignored for ",
                   "this function call and please remove any of the do_pca or ",
                   "npcs parameters and pass to harmony cell_embeddings ",
                   "directly"),
            .frequency = "once", .frequency_id = "do_pca_npcs_warning"
        )
    }
    if (hasArg(lambda)) {
        rlang::warn(
            paste0("Warning: Function parameter lambda is deprecated in this new ",
                   "release of harmony. It will be ignored for this function ",
                   "call and please remove parameter lambda in future ",
                   "function calls. For advanced users, value of lambda ",
                   "should be set by using parameter .options."),
            .frequency = "once", .frequency_id = "lambda_warning")
    }
    if (hasArg(tau)) {
        rlang::warn(
            paste0("Warning: Function parameter tau is deprecated in this new ",
                   "release of harmony. It will be ignored for this function ",
                   "call and please remove parameter tau in future ",
                   "function calls. For advanced users, value of tau ",
                   "should be set by using parameter .options."),
            .frequency = "once", .frequency_id = "tau_warning")
    }
    if (hasArg(block.size)) {
        rlang::warn(
            paste0("Warning: Function parameter block.size is deprecated in this new ",
                   "release of harmony. It will be ignored for this function ",
                   "call and please remove parameter block.size in future ",
                   "function calls. For advanced users, value of block.size ",
                   "should be set by using parameter .options."),
            .frequency = "once", .frequency_id = "block.size_warning")
    }
    if (hasArg(max.iter.harmony)) {
        rlang::warn(
            paste0("Warning: Function parameter max.iter.harmony is replaced ",
                   "with parameter max_iter in this new release of harmony. ",
                   "It will be ignored for this function call and please use ",
                   "parameter max_iter in future function calls."),
            .frequency = "once", .frequency_id = "max.iter.harmony_warning")
    }
    if (hasArg(max.iter.cluster)) {
        rlang::warn(
            paste0("Warning: Function parameter max.iter.cluster is deprecated in this new ",
                   "release of harmony. It will be ignored for this function ",
                   "call and please remove parameter max.iter.cluster in future ",
                   "function calls. For advanced users, value of max.iter.cluster ",
                   "should be set by using parameter .options."),
            .frequency = "once", .frequency_id = "max.iter.cluster_warning")
    }
    if (hasArg(epsilon.cluster)) {
        rlang::warn(
            paste0("Warning: Function parameter epsilon.cluster is deprecated in this new ",
                   "release of harmony. It will be ignored for this function ",
                   "call and please remove parameter epsilon.cluster in future ",
                   "function calls. For advanced users, value of epsilon.cluster ",
                   "should be adjusted by using parameter .options."),
            .frequency = "once", .frequency_id = "epsilon.cluster_warning")
    }
    if (hasArg(epsilon.harmony)) {
        rlang::warn(
            paste0("Warning: Function parameter epsilon.harmony is deprecated in this new ",
                   "release of harmony. It will be ignored for this function ",
                   "call and please remove parameter epsilon.harmony in future ",
                   "function calls. For advanced users, value of epsilon.harmony ",
                   "should be adjusted by using parameter .options."),
            .frequency = "once", .frequency_id = "epsilon.harmony_warning")
    }
    
    # Parameter setting --------------------------------------------------------
    if(early_stop == TRUE){
        epsilon.harmony = 1e-4
    }else{
        epsilon.harmony = -Inf
    }
    max.iter.harmony <- max_iter
    lambda <- .options$lambda
    tau <- .options$tau
    block.size <- .options$block.size
    max.iter.cluster <- .options$max.iter.cluster
    epsilon.cluster <- .options$epsilon.cluster
    if(!(is.null(.options$epsilon.harmony))){
        epsilon.harmony <- .options$epsilon.harmony
    }

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
    
    ## Calculate theta (#covariates) x (#levels)
    theta <- Reduce(c, lapply(seq_len(length(B_vec)), function(b)
        rep(theta[b], B_vec[b])))

    ## Theta scaling

    theta <- theta * (1 - exp(-(N_b / (nclust * tau))^2))
    
    ## RUN HARMONY
    harmonyObj <- new(harmony)
    
    harmonyObj$setup(
        data_mat, phi,
        sigma, theta, max.iter.cluster, epsilon.cluster,
        epsilon.harmony, nclust, block.size, lambda, B_vec, verbose
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
}

#' Main Harmony interface
#' 
#' Use this to run the Harmony algorithm on gene expression or PCA matrix. 
#' 
#' @param data_mat Matrix of normalized gene expession (default) or PCA 
#' embeddings (see do_pca). 
#' Cells can be rows or columns. 
#' @param meta_data Either (1) Dataframe with variables to integrate or (2) 
#' vector with labels. 
#' @param vars_use If meta_data is dataframe, this defined which variable(s) 
#' to remove (character vector).
#' @param do_pca Whether to perform PCA on input matrix. 
#' @param npcs If doing PCA on input matrix, number of PCs to compute. 
#' @param theta Diversity clustering penalty parameter. Specify for each
#'  variable in vars_use Default theta=2. theta=0 does not encourage any 
#'  diversity. Larger values of theta result in more diverse clusters. 
#' @param lambda Ridge regression penalty parameter. Specify for each variable
#'  in vars_use. 
#' Default lambda=1. Lambda must be strictly positive. Smaller values result 
#' in more aggressive correction. 
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
#'  the distance from a cell to cluster centroids. Larger values of sigma 
#'  result in cells assigned to more clusters. Smaller values of sigma make 
#'  soft kmeans cluster approach hard clustering. 
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple 
#' linear regression. 
#' @param tau Protection against overclustering small datasets with large ones.
#'  tau is the expected number of cells per cluster. 
#' @param block.size What proportion of cells to update during clustering.
#'  Between 0 to 1, default 0.05. Larger values may be faster but less accurate
#' @param max.iter.cluster Maximum number of rounds to run clustering at each 
#' round of Harmony. 
#' @param epsilon.cluster Convergence tolerance for clustering round of 
#' Harmony. Set to -Inf to never stop early. 
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round
#'  of Harmony involves one clustering and one correction step. 
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to
#'  never stop early. 
#' @param plot_convergence Whether to print the convergence plot of the 
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging. 
#' @param return_object (Advanced Usage) Whether to return the Harmony object 
#' or only the corrected PCA embeddings. 
#' @param verbose Whether to print progress messages. TRUE to print, 
#' FALSE to suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s). 
#' Cells that have batch variables values matching reference_values will not 
#' be moved.  
#' @param cluster_prior (Advanced Usage) Provides user defined clusters for 
#' cluster initialization. If the number of provided clusters C is less than K, 
#' Harmony will initialize K-C clusters with kmeans. C cannot exceed K.  
#' 
#' @return By default, matrix with corrected PCA embeddings. If return_object 
#' is TRUE, returns the full Harmony object (R6 reference class type). 
#'
#' @export 
#' 
#' @examples
#' 
#' 
#' ## By default, Harmony inputs a normalized gene expression matrix
#' \dontrun{
#' harmony_embeddings <- HarmonyMatrix(exprs_matrix, meta_data, 'dataset')
#' }
#' 
#' ## Harmony can also take a PCA embeddings matrix
#' data(cell_lines_small)
#' pca_matrix <- cell_lines_small$scaled_pcs
#' meta_data <- cell_lines_small$meta_data
#' harmony_embeddings <- HarmonyMatrix(pca_matrix, meta_data, 'dataset', 
#'                                     do_pca=FALSE)
#' 
#' ## Output is a matrix of corrected PC embeddings
#' dim(harmony_embeddings)
#' harmony_embeddings[seq_len(5), seq_len(5)]
#' 
#' ## Finally, we can return an object with all the underlying data structures
#' harmony_object <- HarmonyMatrix(pca_matrix, meta_data, 'dataset', 
#'                                     do_pca=FALSE, return_object=TRUE)
#' dim(harmony_object$Y) ## cluster centroids
#' dim(harmony_object$R) ## soft cluster assignment
#' dim(harmony_object$Z_corr) ## corrected PCA embeddings
#' head(harmony_object$O) ## batch by cluster co-occurence matrix
#' 
HarmonyMatrix <- function(
    data_mat, meta_data, vars_use, do_pca = TRUE,
    npcs = 20, theta = NULL, lambda = NULL, sigma = 0.1, 
    nclust = NULL, tau = 0, block.size = 0.05, 
    max.iter.harmony = 10, max.iter.cluster = 200, 
    epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
    plot_convergence = FALSE, return_object = FALSE, 
    verbose = TRUE, reference_values = NULL, cluster_prior = NULL, Y = NULL
) {
    
    
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
    
    if (do_pca) {
        pca_res <- data_mat %>%
            scaleData() %>% 
            irlba::prcomp_irlba(n = npcs, retx = TRUE, center = FALSE, 
                                scale. = FALSE)
        data_mat <- pca_res$rotation %*% diag(pca_res$sdev)
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
    # determine lamda if null
    if (is.null(lambda)) {
        lambda <- rep(1, length(vars_use))
    } else if (length(lambda) != length(vars_use)) {
        stop('Please specify lambda for each variable')
    }
    
    # determine sigma if it is a scalar
    if (length(sigma) == 1 & nclust > 1) {
        sigma <- rep(sigma, nclust)
    }
    
    ## Pre-compute some useful statistics
    phi <- Reduce(rbind, lapply(vars_use, function(var_use) {
        res <- model.matrix(~0 + as.factor(meta_data[[var_use]]))
        t(res)
    }))

    ## ## number of cells per batch
    N_b <- rowSums(phi)
    ## ## Probability of batches Not used anywhere in R!
    ## Pr_b <- N_b / N

    ## Number of factors per covariate
    B_vec <- Reduce(c, lapply(vars_use, function(var_use) {
        nlevels(as.factor(meta_data[[var_use]]))
    }))
    
    ## Calculate theta (#covariates) x (#levels)
    theta <- Reduce(c, lapply(seq_len(length(B_vec)), function(b)
        rep(theta[b], B_vec[b])))

    ## Theta scaling
    theta <- theta * (1 - exp(-(N_b / (nclust * tau))^2))
    
    ## Calculate lambda (#covariates) x (#levels)
    lambda <- Reduce(c, lapply(seq_len(length(B_vec)), function(b) 
        rep(lambda[b], B_vec[b])))
    lambda_mat <- diag(c(0, lambda))
    
    ## ## TODO: check that each ref val matches exactly one covariate
    ## ## TODO: check that you haven't marked all cells as reference! 
    ## if (!is.null(reference_values)) {
    ##     idx <- which(row.names(phi) %in% reference_values)
    ##     cells_ref <- which(colSums(phi[idx, , drop = FALSE] == 1) >= 1)
    ##     b_keep <- which(!row.names(phi) %in% reference_values)
    ##     phi_moe <- phi[b_keep, , drop = FALSE]
    ##     phi_moe[, cells_ref] <- 0
        
    ##     phi_moe <- rbind(rep(1, N), phi_moe)
    ##     lambda_mat <- lambda_mat[c(1, b_keep + 1), c(1, b_keep + 1)]
    ## }
    
    ## RUN HARMONY
    harmonyObj <- new(harmony)
    
    harmonyObj$setup(
        data_mat, phi,
        sigma, theta, max.iter.cluster, epsilon.cluster,
        epsilon.harmony, nclust, tau, block.size, lambda_mat, verbose
        )
    
    if (!is.null(Y)) {
        harmonyObj$setY(Y)
    }
    
    harmonyObj$init_cluster_cpp(0)
    
    ## if (plot_convergence) graphics::plot(HarmonyConvergencePlot(harmonyObj))

    

    harmonize(harmonyObj, max.iter.harmony, verbose)


    
    ## Return either the R6 Harmony object or the corrected PCA matrix
    if (return_object) {
        return(harmonyObj)
    } else {
        res <- as.matrix(harmonyObj$Z_corr)
        row.names(res) <- row.names(data_mat)
        colnames(res) <- colnames(data_mat)
        if (!cells_as_cols) 
            res <- t(res)
        return(res)      
    }
}


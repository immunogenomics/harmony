#' Main Harmony interface
#' 
#' Use this to run the Harmony algorithm on a PCA matrix. See also RunHarmony to run on Seurat object. 
#' 
#' @param data_mat Matrix of normalized gene expession (default) or PCA embeddings (see do_pca). Cells can be rows or columns. 
#' @param meta_data Dataframe with variables to integrate. Cells must be rows. 
#' @param vars_use Which variable(s) to remove (character vector).
#' @param do_pca Whether to perform PCA on input matrix. 
#' @param npcs If doing PCA on input matrix, number of PCs to compute. 
#' @param theta Diversity clustering penalty parameter. Specify for each variable in vars_use Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters. 
#' @param lambda Ridge regression penalty parameter. Specify for each variable in vars_use Default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction. 
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales the distance from a cell to cluster centroids. Larger values of sigma result in cells assigned to more clusters. Smaller values of sigma make soft kmeans cluster approach hard clustering. 
#' @param nclust Number of clusters in Harmony model. nclust=1 equivalent to simple linear regression. 
#' @param tau Protection against overclustering small datasets with large ones. tau is the expected number of cells per cluster. 
#' @param block.size What proportion of cells to update during clustering. Between 0 to 1, default 0.05. Larger values may be faster but less accurate. 
#' @param max.iter.cluster Maximum number of rounds to run clustering at each round of Harmony. 
#' @param epsilon.cluster Convergence tolerance for clustering round of Harmony. Set to -Inf to never stop early. 
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step. 
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to never stop early. 
#' @param plot_convergence Whether to print the convergence plot of the clustering objective function. TRUE to plot, FALSE to suppress. This can be useful for debugging. 
#' @param return_object (Advanced Usage) Whether to return the Harmony object or only the corrected PCA embeddings. 
#' @param verbose Whether to print progress messages. TRUE to print, FALSE to suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s). Cells that have batch variables values matching reference_values will not be moved.  
#' 
#' @return By default, matrix with corrected PCA embeddings. If return_object is TRUE, returns the full Harmony object (R6 reference class type). 
#' 
#' @export 
#' 
#' examples
#' 
#' pca_embeddings <- harmony::pbmc.small$pca_embeddings
#' meta_data <- harmony::pbmc.small$meta_data
#' 
#' ## first, get the corrected PCA matrix
#' harmony_embeddings <- HarmonyMatrix(pca_embeddings, meta_data, 'stim', theta=2, plot_convergence=TRUE, nclust=50, max.iter.cluster=10, max.iter.harmony=4)
#' head(meta_data)
#' dim(pca_embeddings)
#' dim(harmony_embeddings)
#' 
#' 
#' ## now, let's get the full Harmony object
#' harmony_object <- HarmonyMatrix(pca_embeddings, meta_data, 'stim', theta=2, plot_convergence=TRUE, nclust=50, max.iter.cluster=10, max.iter.harmony=4, return_object = TRUE)
#' dim(harmony_object$Y) ## cluster centroids
#' dim(harmony_object$R) ## soft cluster assignment
#' dim(harmony_object$Z_corr) ## corrected PCA embeddings
#' head(harmony_object$O) ## batch by cluster co-occurence matrix
HarmonyMatrix <- function(data_mat, meta_data, vars_use, do_pca = TRUE, npcs=20, 
                          theta = NULL, lambda = NULL, sigma = 0.1, nclust = 100, 
                          tau = 0, block.size = 0.05, max.iter.harmony = 10, 
                          max.iter.cluster = 200, epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                          plot_convergence = FALSE, 
                          return_object = FALSE, verbose = TRUE, reference_values = NULL) {

    
  ## TODO: check for 
  ##    partially observed batch variables (WARNING)
  ##    batch variables with only 1 level (WARNING)
  ##    if lambda given, check correct length
  ##    if theta given, check correct length
  ##    very small batch size and tau=0: suggest tau>0
  ##    is PCA correct? 
  if (!'data.frame' %in% class(meta_data)) {
    if (length(meta_data) %in% dim(data_mat)) {
      meta_data <- data.frame(batch_variable = meta_data)
      vars_use <- 'batch_variable'
    } else {
      stop('meta_data must be either a data.frame or a vector with batch values for each cell.')
    }
  }
  
  if (is.null(vars_use) | any(!vars_use %in% colnames(meta_data))) {
    stop('Must provides variables to integrate over (e.g. vars_use="stim")')
  }

  if (do_pca) {
    if (ncol(data_mat) != nrow(meta_data)) {
        data_mat <- Matrix::t(data_mat)
    }
      
    pca_res <- data_mat %>% 
      scaleData() %>% 
      irlba::prcomp_irlba(n = npcs, retx = TRUE, center = FALSE, scale. = FALSE)
    data_mat <- pca_res$rotation %*% diag(pca_res$sdev)
  } 

  N <- nrow(meta_data)
  cells_as_cols <- TRUE
  if (ncol(data_mat) != N) {
    if (nrow(data_mat) == N) {
      data_mat <- t(data_mat)
      cells_as_cols <- FALSE
    } else {
      stop("ERROR: Number of labels do not correspond to number of samples in PC matrix.")
    }
  }
  
  if (is.null(theta)) {
    theta <- rep(2, length(vars_use))
  }
  if (is.null(lambda)) {
    lambda <- rep(1, length(vars_use))
  }    
  if (length(sigma) == 1 & nclust > 1) {
    sigma <- rep(sigma, nclust)
  }
  ## TODO: if theta or lambda doesn't match number of variables, fix this
  
    
  ## Pre-compute some useful statistics
  phi <- Reduce(rbind, lapply(vars_use, function(var_use) {t(onehot(meta_data[[var_use]]))}))
  N_b <- rowSums(phi)
  Pr_b <- N_b / N
  B_vec <- Reduce(c, lapply(vars_use, function(var_use) {length(unique(meta_data[[var_use]]))}))
  theta <- Reduce(c, lapply(1:length(B_vec), function(b) rep(theta[b], B_vec[b])))
  theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))
  
  lambda <- Reduce(c, lapply(1:length(B_vec), function(b) rep(lambda[b], B_vec[b])))
  lambda_mat <- diag(c(0, lambda))

  ## TODO: check that each ref val matches exactly one covariate
  ## TODO: check that you haven't marked all cells as reference! 
  if (!is.null(reference_values)) {
    cells_ref <- which(colSums(phi[which(row.names(phi) %in% reference_values), , drop = FALSE] == 1) >= 1)
    b_keep <- which(!row.names(phi) %in% reference_values)
    phi_moe <- phi[b_keep, , drop = FALSE]
    phi_moe[, cells_ref] <- 0

    phi_moe <- rbind(rep(1, N), phi_moe)
    lambda_mat <- lambda_mat[c(1, b_keep + 1), c(1, b_keep + 1)]    
  } else {
    phi_moe <- rbind(rep(1, N), phi)
  }
                             
  ## RUN HARMONY
  harmonyObj <- new(harmony, 0) ## 0 is a dummy variable - will change later
  harmonyObj$setup(
    data_mat, phi, phi_moe, 
    Pr_b, sigma, theta, max.iter.cluster,epsilon.cluster,
    epsilon.harmony, nclust, tau, block.size, lambda_mat, verbose
  )
  init_cluster(harmonyObj)
  harmonize(harmonyObj, max.iter.harmony, verbose)
  if (plot_convergence) plot(HarmonyConvergencePlot(harmonyObj))
  
  ## Return either the R6 Harmony object or the corrected PCA matrix (default)
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

#' Harmony wrapper for Seurat
#' 
#' Use this to run the Harmony algorithm on a Seurat object. See also HarmonyMatrix to run directly on PCA matrix. 
#' 
#' @param object Seurat object. Must have PCA computed. 
#' @param group.by.vars Which variable(s) to remove (character vector).
#' @param theta Diversity clustering penalty parameter. Specify for each variable in group.by.vars. Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters. 
#' @param lambda Ridge regression penalty parameter. Specify for each variable in group.by.vars. Default lambda=1. Lambda must be strictly positive. Smaller values result in more aggressive correction. 
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales the distance from a cell to cluster centroids. Larger values of sigma result in cells assigned to more clusters. Smaller values of sigma make soft kmeans cluster approach hard clustering. 
#' @param nclust Number of clusters in Harmony model. nclust=1 equivalent to simple linear regression. 
#' @param tau Protection against overclustering small datasets with large ones. tau is the expected number of cells per cluster. 
#' @param block.size What proportion of cells to update during clustering. Between 0 to 1, default 0.05. Larger values may be faster but less accurate. 
#' @param max.iter.cluster Maximum number of rounds to run clustering at each round of Harmony. 
#' @param epsilon.cluster Convergence tolerance for clustering round of Harmony. Set to -Inf to never stop early. 
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round of Harmony involves one clustering and one correction step. 
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to never stop early. 
#' @param plot_convergence Whether to print the convergence plot of the clustering objective function. TRUE to plot, FALSE to suppress. This can be useful for debugging. 
#' @param verbose Whether to print progress messages. TRUE to print, FALSE to suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s). Cells that have batch variables values matching reference_values will not be moved.  
#' 
#' @return Seurat object. Harmony dimensions placed into object@dr$harmony. For downstream Seurat analyses, use reduction.use='harmony' and reduction.type='harmony'.
#' 
#' @export 
#' 
#' example
#' 
#' pbmc <- RunHarmony(harmony::pbmc.small.seurat, 'stim', theta=2, plot_convergence=TRUE, nclust=50, max.iter.cluster=10, max.iter.harmony=4)
#' pbmc@dr$harmony@cell.embeddings[1:5, 1:10]
#' pbmc@dr$harmony@gene.loadings[1:5, 1:10]
#' p1 <- DimPlot(object = pbmc, reduction.use = 'harmony', pt.size = .1, group.by = 'stim', do.return = T)
#' p2 <- VlnPlot(object = pbmc, features.plot = 'Harmony1', group.by = 'stim', do.return = TRUE)
#' plot_grid(p1,p2)
#' 
RunHarmony <- function(object, group.by.vars, dims.use, theta = NULL, lambda = NULL, sigma = 0.1, 
                       nclust = 100, tau = 0, block.size = 0.05, max.iter.harmony = 10, 
                       max.iter.cluster = 20, epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                       plot_convergence = FALSE, verbose = TRUE, reference_values = NULL) {
  ## CHECK: PCs should be scaled. Unscaled PCs yield misleading results. 
  ##      sqrt(sum((apply(object@dr$pca@cell.embeddings, 2, sd) - object@dr$pca@sdev) ^ 2))  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Failed to load Seurat library.")
  }
  
  if (!"seurat" %in% class(object)) {
    stop("Must pass a Seurat object to RunHarmony function. Did you mean to use the Harmony function?")
  }
  if (!"pca" %in% names(object@dr)) {
    stop("PCA must be computed before running Harmony.")
  }
  if (missing(dims.use)) {
    dims.use <- 1:ncol(object@dr$pca@cell.embeddings)        
  } else if (!all(dims.use %in% 1:ncol(object@dr$pca@cell.embeddings))) {
    stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")        
  }
  
  missing.vars <- setdiff(group.by.vars, colnames(object@meta.data))
  if (length(missing.vars) > 0) {
    stop(sprintf("ERROR: Primary variable(s) [%s] not in meta.data", paste(missing.vars, collapse = ", ")))
  }
  
  if (verbose) {
    message("Run Harmony")
    message(sprintf("Using top %d PCs", length(dims.use)))    
  }
    
  
  harmonyEmbed <- HarmonyMatrix(object@dr$pca@cell.embeddings, object@meta.data, group.by.vars, 
                                FALSE, 0, 
                                theta, lambda, sigma, nclust, tau, block.size, max.iter.harmony, 
                                max.iter.cluster, epsilon.cluster, epsilon.harmony,
                                plot_convergence, FALSE, verbose, reference_values)
  
  rownames(harmonyEmbed) <- row.names(object@meta.data)
  colnames(harmonyEmbed) <- paste0("harmony_", 1:ncol(harmonyEmbed))
  
  object <- object %>%
    Seurat::SetDimReduction(reduction.type = "harmony", slot = "cell.embeddings", new.data = harmonyEmbed) %>%
    Seurat::SetDimReduction(reduction.type = "harmony", slot = "key", new.data = "Harmony") %>%
    Seurat::ProjectDim(reduction.type = "harmony", replace.dim = T, do.print = F)
  
  return(object)
  
}


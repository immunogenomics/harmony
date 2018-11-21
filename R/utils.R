onehot <- function(x) {
  res <- data.frame(x)
  res$id <- row.names(res)
  res <- res %>%
    dplyr::mutate(dummy = 1) %>% 
    tidyr::spread(x, dummy, fill = 0) %>% 
    dplyr::select(-id) %>%
    as.matrix
  return(res)
}

HarmonyConvergencePlot <- function(harmonyObj) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      message("WARNING: Failed to load ggplot2 package. Cannot plot convergence without this package.")
      return(NA)
    }

    ## ignore initial value (random initialization)
    ## break down kmeans objective into rounds
    obj_fxn <- data.frame(
        kmeans_idx = Reduce(c, lapply(harmonyObj$kmeans_rounds, function(rounds) {1:rounds})),
        harmony_idx = Reduce(c, lapply(1:length(harmonyObj$kmeans_rounds), function(i) {rep(i, harmonyObj$kmeans_rounds[i])})),
        val = tail(harmonyObj$objective_kmeans, -1)
    ) %>%
        tibble::rowid_to_column("idx")
##    data.table(obj_fxn)[, (tail(.SD$val, 1) - head(.SD$val, 1)) / head(.SD$val, 1), by = harmony_idx]
    obj_fxn %>% ggplot(ggplot2::aes(idx, val, col = harmony_idx)) + 
      geom_point(shape = 21) + 
      labs(y = "Objective Function", x = "Iteration Number")
}

#' Encode vector of values from categorical variable as one-hot matrix 
#' 
#' @param
#' @param
#' @param
#' @param
#' @param
#' 
#' @return
#' 
#' @export
#' 
#' @example
#' 
HarmonyMatrix <- function(pc_mat, meta_data, vars_use, theta = NULL, lambda = NULL,
                          sigma = 0.1, alpha = .1, nclust = 100, tau = 0, 
                          block.size = 0.05, max.iter.harmony = 10, 
                          max.iter.cluster = 200, epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                          burn.in.time = 10, plot_convergence = FALSE, 
                          return_object = FALSE) {

    ## TODO: check for 
    ##    partially observed batch variables (WARNING)
    ##    batch variables with only 1 level (WARNING)
    ##    if lambda given, check correct length
    ##    if theta given, check correct length
    ##    very small batch size and tau=0: suggest tau>0
    ##    is PCA correct? 
  
    N <- nrow(meta_data)
    cells_as_cols <- TRUE
    if (ncol(pc_mat) != N) {
        if (nrow(pc_mat) == N) {
            pc_mat <- t(pc_mat)
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

    ## Pre-compute some useful statistics
    phi <- Reduce(rbind, lapply(vars_use, function(var_use) {t(onehot(meta_data[[var_use]]))}))
    N_b <- rowSums(phi)
    Pr_b <- N_b / N
    B_vec <- Reduce(c, lapply(vars_use, function(var_use) {length(unique(meta_data[[var_use]]))}))
    theta <- Reduce(c, lapply(1:length(B_vec), function(b) rep(theta[b], B_vec[b])))
    theta <- theta * (1 - exp(-(N_b / (nclust * tau)) ^ 2))

    lambda <- Reduce(c, lapply(1:length(B_vec), function(b) rep(lambda[b], B_vec[b])))
    lambda_mat <- diag(c(0, lambda))
                              
    ## RUN HARMONY
    harmonyObj <- new(harmony, 0) ## 0 is a dummy variable - will change later
    harmonyObj$setup(
        pc_mat, ## Z
        phi, ## Phi
        Pr_b, 
        sigma, ## sigma
        theta, ## theta
        max.iter.cluster, ## max.iter
        epsilon.cluster, ## kmeans converge.thresh
        epsilon.harmony, ## harmony epsilon
        TRUE, ## correct Z_orig only
        alpha, ## EXPERIMENTAL: alpha, strength of dirichlet prior for DKL penalty
        nclust, ## K
        tau, ## tau (desired cells/cluster)
        block.size, ## model$block.size
#        rep(1, N), ## EXPERIMENTAL FEATURE: each cell gets its own weight
        FALSE, ## do linear correction on Z_cos?
        burn.in.time, ## window size for kmeans convergence
        lambda_mat
    )

                               
    harmonyObj$harmonize(max.iter.harmony)
    if (plot_convergence) plot(HarmonyConvergencePlot(harmonyObj))
    
    ## Return either the R6 Harmony object or the corrected PCA matrix (default)
    if (return_object) {
      return(harmonyObj)
    } else {
      res <- as.matrix(harmonyObj$Z_corr)
      row.names(res) <- row.names(pc_mat)
      colnames(res) <- colnames(pc_mat)
      if (!cells_as_cols) 
          res <- t(res)
      return(res)      
    }
}

#' Encode vector of values from categorical variable as one-hot matrix 
#' 
#' @param
#' 
#' @return
#' 
#' @export
#' 
#' 
RunHarmony <- function(object, group.by.vars, dims.use, theta = NULL, lambda = NULL, sigma = 0.1, alpha = .1,
                       nclust = 100, tau = 0, block.size = 0.05, max.iter.harmony = 10, 
                       max.iter.cluster = 200, epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                       burn.in.time = 10, plot_convergence = FALSE) {
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
  
      
    message("Starting harmony")
    message(sprintf("Using top %d PCs", length(dims.use)))    
    
  
    harmonyEmbed <- HarmonyMatrix(object@dr$pca@cell.embeddings, object@meta.data, group.by.vars, 
                                   theta, lambda, sigma, alpha, nclust, tau, block.size, max.iter.harmony, 
                                   max.iter.cluster, epsilon.cluster, epsilon.harmony,
                                   burn.in.time, plot_convergence)
      
  
    rownames(harmonyEmbed) <- row.names(object@meta.data)
    colnames(harmonyEmbed) <- paste0("harmony_", 1:ncol(harmonyEmbed))
    
    object <- object %>%
                Seurat::SetDimReduction(reduction.type = "harmony", slot = "cell.embeddings", new.data = harmonyEmbed) %>%
                Seurat::SetDimReduction(reduction.type = "harmony", slot = "key", new.data = "Harmony") %>%
                Seurat::ProjectDim(reduction.type = "harmony", replace.dim = T, do.print = F)

    return(object)
    
}

















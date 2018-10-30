onehot <- function(x) {
	data.frame(x) %>%
		tibble::rowid_to_column("id") %>% 
		dplyr::mutate(dummy = 1) %>% 
		tidyr::spread(x, dummy, fill = 0) %>% 
		dplyr::select(-id) %>%
		as.matrix
}

HarmonyConvergencePlot <- function(harmonyObj) {
    ## ignore initial value (random initialization)
    ## break down kmeans objective into rounds
    obj_fxn <- data.frame(
        kmeans_idx = Reduce(c, lapply(harmonyObj$kmeans_rounds, function(rounds) {1:rounds})),
        harmony_idx = Reduce(c, lapply(1:length(harmonyObj$kmeans_rounds), function(i) {rep(i, harmonyObj$kmeans_rounds[i])})),
        val = tail(harmonyObj$objective_kmeans, -1)
    ) %>%
        tibble::rowid_to_column("idx")
##    data.table(obj_fxn)[, (tail(.SD$val, 1) - head(.SD$val, 1)) / head(.SD$val, 1), by = harmony_idx]
    obj_fxn %>% ggplot(aes(idx, val, col = harmony_idx)) + geom_point(shape = 21) + 
    labs(y = "Objective Function", x = "Iteration Number")
}

HarmonyMatrix <- function(pc_mat, batch_labels, batch_labels2 = NULL, theta = 1, theta2 = 1, sigma = 0.1, alpha = .1,
                          nclust = 100, tau = 0, block.size = 0.05, max.iter.harmony = 10, 
                          max.iter.cluster = 200, epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                          burn.in.time = 10, plot_convergence = FALSE) {
    
    cells_as_cols <- TRUE
    if (length(batch_labels) != ncol(pc_mat)) {
        if (length(batch_labels) == nrow(pc_mat)) {
            pc_mat <- t(pc_mat)
            cells_as_cols <- FALSE
        } else {
            stop("ERROR: Number of labels do not correspond to number of samples in PC matrix.")
        }
    }    

  
    ## RUN HARMONY
    harmonyObj <- new(harmony, 0) ## 0 is a dummy variable - will change later
    batch_mat <- t(onehot(batch_labels))
    harmonyObj$setup(
        pc_mat, ## Z
        batch_mat, ## Phi
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
        rep(1, length(batch_labels)), ## EXPERIMENTAL FEATURE: each cell gets its own weight
        FALSE, ## do linear correction on Z_cos?
        rep(1, nrow(batch_mat)), ## EXPERIMENTAL FEATURE: only correct certain batches
        burn.in.time ## window size for kmeans convergence
    )
  
    ## OPTIONAL: 2nd level of batch defined
    if (!is.null(batch_labels2)) {
      harmonyObj$setup_batch2(t(onehot(batch_labels2)), theta2, tau)
    }
  
    harmonyObj$harmonize(max.iter.harmony)
    if (plot_convergence) plot(HarmonyConvergencePlot(harmonyObj))
    
    res <- as.matrix(harmonyObj$Z_corr)
    row.names(res) <- row.names(pc_mat)
    colnames(res) <- colnames(pc_mat)
    if (!cells_as_cols) 
        res <- t(res)
    return(res)
}



RunHarmony <- function(object, group.by, dims.use, group.by.secondary = NULL, theta = 1, theta2 = 1, sigma = 0.1, alpha = .1,
                       nclust = 100, tau = 0, block.size = 0.05, max.iter.harmony = 10, 
                       max.iter.cluster = 200, epsilon.cluster = 1e-5, epsilon.harmony = 1e-4, 
                       burn.in.time = 10, plot_convergence = FALSE) {
    ## CHECK: PCs should be scaled. Unscaled PCs yield misleading results. 
    ##      sqrt(sum((apply(object@dr$pca@cell.embeddings, 2, sd) - object@dr$pca@sdev) ^ 2))  

    
    if (!"seurat" %in% class(object)) {
        stop("This Function is meant to be run on a Seurat object!")
    }    
    if (!"pca" %in% names(object@dr)) {
        stop("PCA must be computed before running Harmony.")
    }
    if (missing(dims.use)) {
        dims.use <- 1:ncol(object@dr$pca@cell.embeddings)        
    } else if (!all(dims.use %in% 1:ncol(object@dr$pca@cell.embeddings))) {
        stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")        
    }
    
    if (!group.by %in% colnames(object@meta.data)) {
      stop(sprintf("ERROR: Primary integration variable [%s] is not in meta.data"))
    }
    if (!is.null(group.by.secondary)) {
        if (!group.by.secondary %in% colnames(object@meta.data))
            stop(sprintf("ERROR: Secondary integration variable [%s] is not in meta.data"))
    }
  
  
    nbatches <- sum(table(object@meta.data[[group.by]]) > 10)
    if (nbatches < 2) {
        stop("Detected fewer than 2 batches with at least 10 cells each. Did you mean to use a different group.by variable?")
    }   
    
    message("Starting harmony")
    message(sprintf("Found %d datasets to integrate", nbatches))
    message(sprintf("Using top %d PCs", length(dims.use)))    
    
    if (!is.null(group.by.secondary)) {
      batches_secondary <- object@meta.data[[group.by.secondary]]
    } else {
      batches_secondary <- NULL
    }
  
    harmonyEmbed <- HarmonyMatrix(object@dr$pca@cell.embeddings, object@meta.data[[group.by]], batches_secondary, 
                                   theta, theta2, sigma, alpha, nclust, tau, block.size, max.iter.harmony, 
                                   max.iter.cluster, epsilon.cluster, epsilon.harmony,
                                   burn.in.time, plot_convergence)
      
  
    rownames(harmonyEmbed) <- row.names(object@meta.data)
    colnames(harmonyEmbed) <- paste0("harmony_", 1:ncol(harmonyEmbed))
    
    object %<>% SetDimReduction(reduction.type = "harmony", slot = "cell.embeddings", new.data = harmonyEmbed) %>%
                SetDimReduction(reduction.type = "harmony", slot = "key", new.data = "Harmony") %>%
                ProjectDim(reduction.type = "harmony", replace.dim = T, do.print = F)

    return(object)
    
}























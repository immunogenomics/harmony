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

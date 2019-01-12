#' @title RunHarmony
#'
#' @param object Seurat object to harmonize
#' @param group.by Metadata column to use as a grouping variable when harmonizing
#' @param dims.use PCA dimensions to use (Default: all available)
#' @param assay.use (Default: "RNA")
#' @param group.by.secondary Second metadata column to use as a grouping variable
#' @param theta 
#' @param theta2 
#' @param sigma 
#' @param alpha 
#' @param nclust 
#' @param tau 
#' @param block.size 
#' @param max.iter.harmony 
#' @param max.iter.cluster 
#' @param epsilon.cluster 
#' @param epsilon.harmony 
#' @param burn.in.time 
#' @param plot_convergence
#' @param reduction.save (Default: "harmony")
#'
#' @import Seurat
#' @importFrom glue glue
#' @export
#' @return
#' 
#' @examples
RunHarmony <- function(object, ...) {
  UseMethod("RunHarmony")
}

#' @rdname RunHarmony
#' @method RunHarmony Seurat
#' @return
#' @export
RunHarmony.Seurat <- function(object, 
                              group.by, 
                              dims.use,
                              assay.use = "RNA",
                              group.by.secondary = NULL, 
                              theta = 1, 
                              theta2 = 1, 
                              sigma = 0.1, 
                              alpha = .1,
                              nclust = 100, 
                              tau = 0, 
                              block.size = 0.05, 
                              max.iter.harmony = 10, 
                              max.iter.cluster = 200, 
                              epsilon.cluster = 1e-5, 
                              epsilon.harmony = 1e-4, 
                              burn.in.time = 10, 
                              plot_convergence = FALSE,
                              reduction.save = "harmony") {
  
  if (!"pca" %in% names(object)) {
    message("PCA must be computed before running Harmony.  Since it is missing, will not compute PCA...")
    object <- RunPCA(object, assay = assay.use)
  }
  if (missing(dims.use)) {
    dims.use <- 1:length(object[["pca"]])        
  } else if (!all(dims.use %in% 1:length(object[['pca']]))) {
    stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")        
  }
  
  if (!(group.by) %in% colnames(object@meta.data)) {
    stop(glue("ERROR: Primary integration variable {group.by} is not in meta.data"))
  }
  if (!is.null(group.by.secondary)) {
    if (!group.by.secondary %in% colnames(object@meta.data))
      stop(glue("ERROR: Secondary integration variable {group.by.secondary} is not in meta.data"))
  }
  
  nbatches <- sum(table(object[[group.by]]) > 10)
  if (nbatches < 2) {
    stop("Detected fewer than 2 batches with at least 10 cells each. Did you mean to use a different group.by variable?")
  }   
  
  message(glue("Starting harmony\nFound {nbatches} datasets to integrate\nUsing top {length(dims.use)} PCs"))
  
  if (!is.null(group.by.secondary)) {
    
    batches_secondary <- object[[group.by.secondary]][,1]
  } else {
    batches_secondary <- NULL
  }
  
  ce <- Embeddings(object = object, 
                   reduction = "pca", 
                   assay = assay.use)
  harmonyEmbed <- HarmonyMatrix(pc_mat = ce, 
                                batch_labels = object[[group.by]][,1], 
                                batch_labels2 = batches_secondary, 
                                theta = theta, 
                                theta2 = theta2, 
                                sigma = sigma, 
                                alpha = alpha, 
                                nclust = nclust, 
                                tau = tau, 
                                block.size = block.size, 
                                max.iter.harmony = max.iter.harmony, 
                                max.iter.cluster = max.iter.cluster, 
                                epsilon.cluster = epsilon.cluster, 
                                epsilon.harmony = epsilon.harmony,
                                burn.in.time = burn.in.time, 
                                plot_convergence = plot_convergence)
  
  rownames(harmonyEmbed) <- colnames(object)
  colnames(harmonyEmbed) <- glue("harmony_{1:ncol(harmonyEmbed)}") %>% as.character()
  
  harmonydata <- CreateDimReducObject(embeddings = harmonyEmbed, assay = assay.use, key="harmony")
  object[[reduction.save]] <- harmonydata
  
  return(object)
}

#' @rdname RunHarmony
#' @method RunHarmony seurat
#' @return
#' @export
RunHarmony.seurat <- function(object,
                              group.by,
                              dims.use,
                              group.by.secondary = NULL,
                              theta = 1,
                              theta2 = 1,
                              sigma = 0.1,
                              alpha = .1,
                              nclust = 100,
                              tau = 0,
                              block.size = 0.05,
                              max.iter.harmony = 10,
                              max.iter.cluster = 200,
                              epsilon.cluster = 1e-5,
                              epsilon.harmony = 1e-4,
                              burn.in.time = 10,
                              plot_convergence = FALSE) {
  ## CHECK: PCs should be scaled. Unscaled PCs yield misleading results.
  ##      sqrt(sum((apply(object@dr$pca@cell.embeddings, 2, sd) - object@dr$pca@sdev) ^ 2))


  if (!"pca" %in% names(object@dr)) {
    stop("PCA must be computed before running Harmony.")
  }
  if (missing(dims.use)) {
    dims.use <- 1:ncol(PCAEmbed(object))
  } else if (!all(dims.use %in% 1:ncol(PCAEmbed(object)))) {
    stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")
  }

  if (!group.by %in% colnames(object@meta.data)) {
    stop(glue("ERROR: Primary integration variable {group.by} is not in meta.data"))
  }
  if (!is.null(group.by.secondary)) {
    if (!group.by.secondary %in% colnames(object@meta.data)) {
      stop(glue("ERROR: Secondary integration variable {group.by.secondary} is not in meta.data"))
    }
  }


  nbatches <- sum(table(object@meta.data[[group.by]]) > 10)
  if (nbatches < 2) {
    stop("Detected fewer than 2 batches with at least 10 cells each. Did you mean to use a different group.by variable?")
  }

  message("Starting harmony")
  message(glue("Found {nbatches} datasets to integrate"))
  message(glue("Using top {length(dims.use)} PCs"))

  if (!is.null(group.by.secondary)) {
    batches_secondary <- object@meta.data[[group.by.secondary]]
  } else {
    batches_secondary <- NULL
  }

  ce <- PCAEmbed(object = object)
  
  harmonyEmbed <- HarmonyMatrix(ce, 
                                object@meta.data[[group.by]], 
                                batches_secondary,
                                theta, 
                                theta2, 
                                sigma, 
                                alpha, 
                                nclust, 
                                tau, 
                                block.size, 
                                max.iter.harmony,
                                max.iter.cluster, 
                                epsilon.cluster, 
                                epsilon.harmony,
                                burn.in.time, 
                                plot_convergence)

  rownames(harmonyEmbed) <- row.names(object@meta.data)
  colnames(harmonyEmbed) <- glue("harmony_{1:ncol(harmonyEmbed)}")

  object %<>% SetDimReduction(reduction.type = "harmony", slot = "cell.embeddings", new.data = harmonyEmbed) %>%
    SetDimReduction(reduction.type = "harmony", slot = "key", new.data = "Harmony") %>%
    ProjectDim(reduction.type = "harmony", replace.dim = T, do.print = F)

  return(object)
}
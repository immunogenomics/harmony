#' Title
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
#' @importFrom glue glue
#' @importFrom Seurat Embeddings CreateDimReducObject 
#' @return
#' @export
#'
#' @examples
RunHarmony <- function(object, 
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
  ## CHECK: PCs should be scaled. Unscaled PCs yield misleading results. 
  ##      sqrt(sum((apply(object@dr$pca@cell.embeddings, 2, sd) - object@dr$pca@sdev) ^ 2))  
  
  
  group.by
  if (!"Seurat" %in% class(object)) {
    stop("This Function is meant to be run on a Seurat object!")
  }    
  if (!"pca" %in% names(object)) {
    message("PCA must be computed before running Harmony.  Since it is missing, will not compute PCA...")
    object <- RunPCA(object, assay = assay.use)
  }
  if (missing(dims.use)) {
    dims.use <- 1:ncol(object[["pca"]])        
  } else if (!all(dims.use %in% 1:length(names(object[['pca']])))) {
    stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")        
  }
  
  if (!(group.by) %in% colnames(object@meta.data)) {
    stop(glue("ERROR: Primary integration variable {group.by} is not in meta.data"))
  }
  if (!is.null(group.by.secondary)) {
    if (!group.by.secondary %in% colnames(object@meta.data))
      stop(glue("ERROR: Secondary integration variable {group.by.secondary} is not in meta.data"))
  }
  
  nbatches <- sum(table(object@meta.data[,group.by]) > 10)
  if (nbatches < 2) {
    stop("Detected fewer than 2 batches with at least 10 cells each. Did you mean to use a different group.by variable?")
  }   
  
  message("Starting harmony")
  message(glue("Found {nbatches} datasets to integrate"))
  message(glue("Using top {length(dims.use)} PCs"))    
  
  if (!is.null(group.by.secondary)) {
    
    batches_secondary <- object@meta.data[,group.by.secondary]
  } else {
    batches_secondary <- NULL
  }
  
  ce <- Embeddings(object = object, reduction = "pca", assay = assay.use)
  harmonyEmbed <- HarmonyMatrix(ce, object@meta.data[,group.by], batches_secondary, 
                                theta, theta2, sigma, alpha, nclust, tau, block.size, max.iter.harmony, 
                                max.iter.cluster, epsilon.cluster, epsilon.harmony,
                                burn.in.time, plot_convergence)
  
  rownames(harmonyEmbed) <- colnames(object)
  colnames(harmonyEmbed) <- glue("harmony_{1:ncol(harmonyEmbed)}") %>% as.character()
  
  harmonydata <- CreateDimReducObject(embeddings = harmonyEmbed, assay = assay.use, key="harmony")
  object[[reduction.save]] <- harmonydata
  
  return(object)
  
}

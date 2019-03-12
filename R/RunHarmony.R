#' @title RunHarmony
#'
#' @param object Data object to harmonize
#' @param grouping.var Metadata column(s) to use as a grouping variable when harmonizing
#' @param dims.use PCA dimensions to use (Default: all available)
#' @param assay.use (Default: "RNA")
#' @param theta
#' @param lambda
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
#' @importFrom Hmisc %nin%
#' @export
#' @return
#'
#' @examples
RunHarmony <- function(object, ...) {
  UseMethod("RunHarmony")
}

#' @rdname RunHarmony
#' @method RunHarmony Seurat
#' @import Seurat
#' @return
#' @export
RunHarmony.Seurat <- function(object,
                              grouping.vars,
                              dims.use,
                              assay.use = "RNA",
                              theta = NULL,
                              lambda = NULL,
                              sigma = 0.1,
                              alpha = 0.1,
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
  } else if (!all(dims.use %in% 1:length(object[["pca"]]))) {
    stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")
  }

  if (!all(grouping.vars %in% colnames(object@meta.data))) {
    stop(glue("ERROR: The integration variable(s) {grouping.vars[which(grouping.vars %nin% colnames(object@meta.data))]} is not in meta.data"))
  }

  message(glue("Starting harmony\nUsing top {length(dims.use)} PCs"))

  ce <- Embeddings(
    object = object,
    reduction = "pca",
    assay = assay.use
  )[,dims.use]

  harmonyEmbed <- HarmonyMatrix(
    pc_mat = ce,
    meta_data = object@meta.data,
    vars_use = grouping.vars,
    theta = theta,
    lambda = lambda,
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
    plot_convergence = plot_convergence,
    return_object = FALSE,
    init_mode = 'kmeans'
  )

  rownames(harmonyEmbed) <- colnames(object)
  colnames(harmonyEmbed) <- glue("harmony_{1:ncol(harmonyEmbed)}") %>% as.character()

  harmonydata <- CreateDimReducObject(embeddings = harmonyEmbed, assay = assay.use, key = "harmony")
  object[[reduction.save]] <- harmonydata

  return(object)
}

#' @rdname RunHarmony
#' @method RunHarmony seurat
#' @import Seurat
#' @return
#' @export
RunHarmony.seurat <- function(object,
                              grouping.vars,
                              dims.use,
                              theta = NULL,
                              lambda = NULL,
                              sigma = 0.1,
                              alpha = 0.1,
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


  if (!"pca" %in% names(object@dr)) {
    message("PCA must be computed before running Harmony.  Since it is missing, will not compute PCA...")
    object <- RunPCA(object)
  }
  if (missing(dims.use)) {
    dims.use <- 1:ncol(PCAEmbed(object))
  } else if (!all(dims.use %in% 1:ncol(PCAEmbed(object)))) {
    stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")
  }

  if (!all(grouping.vars %in% colnames(object@meta.data))) {
    stop(glue("ERROR: The integration variable(s) {grouping.vars[which(grouping.vars %nin% colnames(object@meta.data))]} is not in meta.data"))
  }

  message(glue("Starting harmony\nUsing top {length(dims.use)} PCs"))

  ce <- PCAEmbed(object = object)

  harmonyEmbed <- HarmonyMatrix(
    pc_mat = ce,
    meta_data = object@meta.data,
    vars_use = grouping.vars,
    theta = theta,
    lambda = lambda,
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
    plot_convergence = plot_convergence,
    return_object = FALSE,
    init_mode = 'kmeans'
  )

  rownames(harmonyEmbed) <- row.names(object@meta.data)
  colnames(harmonyEmbed) <- glue("harmony_{1:ncol(harmonyEmbed)}")

  object %<>% SetDimReduction(reduction.type = reduction.save, slot = "cell.embeddings", new.data = harmonyEmbed) %>%
    SetDimReduction(reduction.type = reduction.save, slot = "key", new.data = "Harmony") %>%
    ProjectDim(reduction.type = reduction.save, replace.dim = T, do.print = F)

  return(object)
}

#' @rdname RunHarmony
#' @method RunHarmony SingleCellExperiment
#' @import SingleCellExperiment
#' @return
#' @export
#'
RunHarmony.SingleCellExperiment <- function(object,
                                            grouping.vars,
                                            dims.use,
                                            theta = NULL,
                                            lambda = NULL,
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


  if (!"PCA" %in% reducedDimNames(object)) {
    stop("PCA must be computed before running Harmony.")
  }
  if (missing(dims.use)) {
    dims.use <- 1:ncol(reducedDim(object, "PCA"))
  } else if (!all(dims.use %in% 1:ncol(reducedDim(object, "PCA")))) {
    stop("Trying to use more dimensions than computed with PCA. Rereun PCA with more dimensions or use fewer PCs.")
  }

  if (!all(grouping.vars %in% colnames(colData(object)))) {
    stop(glue("ERROR: The integration variable(s) {grouping.vars[which(grouping.vars %nin% colnames(colData(object)))]} is not in meta.data"))
  }

  message(glue("Starting harmony\nUsing top {length(dims.use)} PCs"))

  ce <- reducedDim(
    x = object,
    type = "PCA"
  )

  harmonyEmbed <- HarmonyMatrix(
    pc_mat = ce,
    meta_data = colData(object),
    vars_use = grouping.vars,
    theta = theta,
    lambda = lambda,
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
    plot_convergence = plot_convergence,
    return_object = FALSE,
    init_mode = 'kmeans'
  )

  rownames(harmonyEmbed) <- row.names(colData(object))
  colnames(harmonyEmbed) <- glue("harmony_{1:ncol(harmonyEmbed)}")

  reducedDim(object, "HARMONY") <- harmonyEmbed

  return(object)
}

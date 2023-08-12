#' Harmony single cell integration
#'
#' Run Harmony algorithm with Seurat and SingleCellAnalysis pipelines.
#'
#' @param object Pipeline object. Must have PCA computed.
#' @param group.by.vars Which variable(s) to remove (character vector).
#' @param dims.use Which PCA dimensions to use for Harmony. By default, use all
#' @param theta Diversity clustering penalty parameter. Specify for each
#' variable in group.by.vars. Default theta=2. theta=0 does not encourage any
#'  diversity. Larger values of theta result in more diverse clusters.
#' Smaller values result in more aggressive correction.
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales
#' the distance from a cell to cluster centroids. Larger values of sigma result
#'  in cells assigned to more clusters. Smaller values of sigma make soft
#'  kmeans cluster approach hard clustering.
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple
#'  linear regression.
#' @param max_iter Maximum number of rounds to run Harmony. One round
#'  of Harmony involves one clustering and one correction step.
#' @param early_stop When to stop Harmony iteration before reaching max_iter 
#' when the change in objectie function is small enough (< 1e-4)
#' @param plot_convergence Whether to print the convergence plot of the
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging.
#' @param verbose Whether to print progress messages. TRUE to print, FALSE to
#'  suppress.
#' @param reduction.save Keyword to save Harmony reduction. Useful if you want
#' to try Harmony with multiple parameters and save them as e.g.
#' 'harmony_theta0', 'harmony_theta1', 'harmony_theta2'
#' @param ... other parameters
#'
#'
#' @rdname RunHarmony
#' @export
RunHarmony <- function(object, group.by.vars, ...) {
    UseMethod("RunHarmony")
}



#' @rdname RunHarmony
#' @param reduction Name of dimension reduction to use. Default is PCA.
#' @param project.dim Project dimension reduction loadings. Default TRUE.
#' @return Seurat (version 3) object. Harmony dimensions placed into
#' dimensional reduction object harmony. For downstream Seurat analyses,
#' use reduction='harmony'.
#' @export
RunHarmony.Seurat <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  sigma = 0.1,
  nclust = NULL,
  max_iter = 10,
  early_stop = TRUE,
  plot_convergence = FALSE,
  verbose = TRUE,
  reduction.save = "harmony",
  project.dim = TRUE,
  .options = harmony_options(),
  ...
) {
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Harmony on a Seurat object requires Seurat")
  }
  if (!reduction %in% Seurat::Reductions(object = object)) {
      stop(paste(reduction, "cell embeddings not found in Seurat object.",
                 "For a Seurat preprocessing walkthrough, please refer to the vignette"))
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rerun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = Seurat::Cells(x = object[[reduction]])
  )

  harmonyEmbed <- RunHarmony(
    data_mat = embedding[, dims.use],
    meta_data = metavars_df,
    vars_use = group.by.vars,
    theta = theta,
    sigma = sigma,
    nclust = nclust,
    max_iter = max_iter,
    early_stop = early_stop,
    plot_convergence= plot_convergence,
    return_object = FALSE,
    verbose = verbose,
    .options = .options,
    ...
  )

  reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
  rownames(harmonyEmbed) <- rownames(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))

  object[[reduction.save]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = Seurat::DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}



#' @rdname RunHarmony
#' @return SingleCellExperiment object. After running RunHarmony, the corrected
#' cell embeddings can be accessed with reducedDim(object, "Harmony").
#' @export
RunHarmony.SingleCellExperiment <- function(
    object,
    group.by.vars,
    dims.use = NULL,
    theta = NULL,
    sigma = 0.1,
    nclust = NULL,
    max_iter = 10,
    early_stop = TRUE,
    plot_convergence = FALSE,
    verbose = TRUE,
    reduction.save = "HARMONY",
    .options = harmony_options(),
    ...
) {

    ## Get PCA embeddings
    if (!"PCA" %in% SingleCellExperiment::reducedDimNames(object)) {
        stop("PCA must be computed before running Harmony.")
    }
    pca_embedding <- SingleCellExperiment::reducedDim(object, "PCA")
    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    }

    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    }
    dims_avail <- seq_len(ncol(pca_embedding))
    if (!all(dims.use %in% dims_avail)) {
        stop("trying to use more dimensions than computed with PCA. Rerun
            PCA with more dimensions or use fewer PCs")
    }

    metavars_df <- SingleCellExperiment::colData(object)
    if (!all(group.by.vars %in% colnames(metavars_df))) {
        stop('Trying to integrate over variables missing in colData')
    }

    harmonyEmbed <- RunHarmony(
        data_mat = pca_embedding[, dims.use], # is here an error? quick fix 
        meta_data = metavars_df,
        vars_use = group.by.vars,
        theta = theta,
        sigma = sigma,
        nclust = nclust,
        max_iter = max_iter,
        early_stop = early_stop,
        plot_convergence= plot_convergence,
        return_object = FALSE,
        verbose = verbose,
        .options = .options,
        ...
    )
   

    rownames(harmonyEmbed) <- row.names(metavars_df)
    colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))
    SingleCellExperiment::reducedDim(object, reduction.save) <- harmonyEmbed

    return(object)
}

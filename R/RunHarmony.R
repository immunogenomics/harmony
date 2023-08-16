
#' Run harmony algorithm generic function
#'
#' This is a generic that provides wrappers for Seurat and SingleCellExperiment
#' objects. Also, it allows harmony standalone with a matrix and a metadata
#' dataframe.
#' 
#' @rdname RunHarmony
#' @export
RunHarmony <- function(...) {
    UseMethod("RunHarmony")
}



#' Applies harmony on a Seurat object cell embedding.
#'
#' @rdname RunHarmony
#' @param object the Seurat object. It needs to have the appropriate slot
#'     of cell embeddings precomputed.
#' @param group.by.vars the name(s) of covariates that harmony will remove
#'     its effect on the data.
#' @param reduction.use Name of dimension reduction to use. Default is pca.
#' @param dims.use indices of the cell embedding features to be used
#' @param reduction.save the name of the new slot that is going to be created by
#'     harmony. By default, harmony.
#' @param project.dim Project dimension reduction loadings. Default TRUE.
#' @param ... harmony algorithm parameters to be passed on RunHarmony.default
#' @return Seurat object. Harmony dimensions placed into a new slot in the Seurat
#' object according to the reduction.save. For downstream Seurat analyses,
#' use reduction='harmony'.
#' 
#' @export
RunHarmony.Seurat <- function(
  object,
  group.by.vars,
  reduction.use = 'pca',
  dims.use = NULL,
  verbose = TRUE,
  reduction.save = "harmony",
  project.dim = TRUE,
  ...
) {
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Harmony on a Seurat object requires Seurat")
  }
  if (!reduction.use %in% Seurat::Reductions(object = object)) {
      stop(paste(reduction.use, "cell embeddings not found in Seurat object.",
                 "For a Seurat preprocessing walkthrough, please refer to the vignette"))
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction.use)
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
    cells = Seurat::Cells(x = object[[reduction.use]])
  )

  harmonyEmbed <- RunHarmony(
    data_mat = embedding[, dims.use],
    meta_data = metavars_df,
    vars_use = group.by.vars,
    return_object = FALSE,
    ...
  )

  reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
  rownames(harmonyEmbed) <- rownames(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))

  object[[reduction.save]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = Seurat::DefaultAssay(object = object[[reduction.use]]),
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



#' Applies harmony on PCA cell embeddings of a SingleCellExperiment.
#' 
#' @param object SingleCellExperiment with the PCA reducedDim cell embeddings populated 
#' @param group.by.vars the name(s) of covariates that harmony will remove
#'     its effect on the data.
#' @param dims.use a vector of indices that allows only selected cell embeddings
#'     features to be used.
#' @param verbose enable verbosity 
#' @param reduction.save the name of the new slot that is going to be created by
#'     harmony. By default, HARMONY.
#' @param ... harmony algorithm parameters to be passed on RunHarmony.default 
#' @rdname RunHarmony
#' @return SingleCellExperiment object. After running RunHarmony, the corrected
#' cell embeddings can be accessed with reducedDim(object, "Harmony").
#' @export
RunHarmony.SingleCellExperiment <- function(
    object,
    group.by.vars,
    dims.use = NULL,
    verbose = TRUE,
    reduction.save = "HARMONY",
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
        return_object = FALSE,
        verbose = verbose,
        ...
    )
   

    rownames(harmonyEmbed) <- row.names(metavars_df)
    colnames(harmonyEmbed) <- paste0(reduction.save, "_", seq_len(ncol(harmonyEmbed)))
    SingleCellExperiment::reducedDim(object, reduction.save) <- harmonyEmbed

    return(object)
}

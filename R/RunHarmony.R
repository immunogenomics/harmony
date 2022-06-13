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
#' @param lambda Ridge regression penalty parameter. Specify for each variable 
#' in group.by.vars. Default lambda=1. Lambda must be strictly positive. 
#' Smaller values result in more aggressive correction. 
#' @param sigma Width of soft kmeans clusters. Default sigma=0.1. Sigma scales 
#' the distance from a cell to cluster centroids. Larger values of sigma result
#'  in cells assigned to more clusters. Smaller values of sigma make soft 
#'  kmeans cluster approach hard clustering. 
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple
#'  linear regression. 
#' @param tau Protection against overclustering small datasets with large ones.
#'  tau is the expected number of cells per cluster. 
#' @param block.size What proportion of cells to update during clustering.
#'  Between 0 to 1, default 0.05. Larger values may be faster but less accurate 
#' @param max.iter.cluster Maximum number of rounds to run clustering at each 
#' round of Harmony. 
#' @param epsilon.cluster Convergence tolerance for clustering round of Harmony
#'  Set to -Inf to never stop early. 
#' @param max.iter.harmony Maximum number of rounds to run Harmony. One round
#'  of Harmony involves one clustering and one correction step. 
#' @param epsilon.harmony Convergence tolerance for Harmony. Set to -Inf to
#'  never stop early. 
#' @param plot_convergence Whether to print the convergence plot of the 
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging. 
#' @param verbose Whether to print progress messages. TRUE to print, FALSE to
#'  suppress.
#' @param reference_values (Advanced Usage) Defines reference dataset(s). Cells
#' that have batch variables values matching reference_values will not be moved
#' @param reduction.save Keyword to save Harmony reduction. Useful if you want
#' to try Harmony with multiple parameters and save them as e.g. 
#' 'harmony_theta0', 'harmony_theta1', 'harmony_theta2'
#' @param assay.use (Seurat V3 only) Which assay to Harmonize with (RNA 
#' by default). 
#' @param ... other parameters 
#' 
#' @examples
#' 
#' ## Seurat Version 2
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'     pkg_version <- packageVersion('Seurat')
#'     if (pkg_version >= "2.0" & pkg_version < "3.0") {
#'         data(cell_lines_small_seurat_v2)
#'         seuratObject <- RunHarmony(cell_lines_small_seurat_v2, 'dataset', 
#'                                    lambda = .1, verbose = FALSE)
#' 
#'         ## Harmony cell embeddings
#'         harmony_embedding <- Seurat::GetCellEmbeddings(
#'             seuratObject, 'harmony'
#'         )
#'         harmony_embedding[seq_len(5), seq_len(5)] 
#' 
#'         ## Harmony gene loadings
#'         harmony_loadings <- Seurat::GetGeneLoadings(
#'             seuratObject, 'harmony'
#'         )
#'         harmony_loadings[seq_len(5), seq_len(5)] 
#' 
#'         p1 <- Seurat::DimPlot(seuratObject, reduction.use = 'harmony', 
#'                       group.by = 'dataset', do.return = TRUE)
#'         p2 <- Seurat::VlnPlot(seuratObject, features.plot = 'Harmony1', 
#'                       group.by = 'dataset', do.return = TRUE)
#'         cowplot::plot_grid(p1,p2)
#'     }
#' }
#' 
#' ## Seurat Version 3
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'     pkg_version <- packageVersion('Seurat')
#'     if (pkg_version >= "3.0" & pkg_version < "4.0") {
#'         data(cell_lines_small_seurat_v3)
#'         seuratObject <- RunHarmony(cell_lines_small_seurat_v3, 'dataset', 
#'                                    lambda = .1, verbose = FALSE)
#'         ## Harmony cell embeddings
#'         harmony_embedding <- Seurat::Embeddings(seuratObject, 'harmony')
#'         harmony_embedding[seq_len(5), seq_len(5)] 
#'         ## Harmony gene loadings
#'         harmony_loadings <- Seurat::Loadings(seuratObject, 'harmony')
#'         harmony_loadings[seq_len(5), seq_len(5)] 
#' 
#'         p1 <- Seurat::DimPlot(seuratObject, reduction = 'harmony', 
#'                               group.by = 'dataset', do.return = TRUE)
#'         p2 <- Seurat::VlnPlot(seuratObject, features = 'harmony_1', 
#'                               group.by = 'dataset', do.return = TRUE)
#'         cowplot::plot_grid(p1, p2)
#'     }
#' }
#' 
#' ## SingleCellExperiment 
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#' 
#'     data(cell_lines_small_sce)
#'     sceObject <- RunHarmony(cell_lines_small_sce, 'dataset',
#'                             lambda = .1, verbose = FALSE)
#' 
#'     ## Harmony cell embeddings
#'     harmony_embedding <- SingleCellExperiment::reducedDim(
#'         sceObject, 'HARMONY'
#'     )
#'     harmony_embedding[seq_len(5), seq_len(5)]
#' 
#'     ## Plot the Harmonized embeddings
#'     ## Colored by batch and cell type
#'     SingleCellExperiment::reducedDim(sceObject, 'HARMONY') %>% 
#'         cbind(SingleCellExperiment::colData(sceObject)) %>% 
#'         data.frame() %>% 
#'     tidyr::gather(key, val, dataset, cell_type) %>%
#'     dplyr::mutate(key = dplyr::case_when(
#'         key == 'dataset' ~ 'Dataset batch', 
#'         key == 'cell_type' ~ 'Known cell types'
#'     )) %>% 
#'     dplyr::sample_frac(1L, FALSE) %>% 
#'     ggplot2::ggplot(ggplot2::aes(x = harmony_1, 
#'                                  y = harmony_2, 
#'                                  color = val)) + 
#'         ggplot2::geom_point() + 
#'         ggplot2::facet_wrap(~key) + 
#'         ggplot2::theme_test(base_size = 12) + 
#'         ggplot2::labs(title = 'Cell embeddings after Harmony')
#' }
#' 
#' @rdname RunHarmony
#' @export 
RunHarmony <- function(object, group.by.vars, ...) {
    UseMethod("RunHarmony")
}

#' @rdname RunHarmony
#' @return Seurat (version 2) object. Harmony dimensions placed into 
#' dimensional reduction object harmony. For downstream Seurat analyses, 
#' use reduction.use='harmony' and reduction.type='harmony'.
#' @export
RunHarmony.seurat <- function(
    object, 
    group.by.vars, 
    dims.use = NULL, 
    theta = NULL, 
    lambda = NULL, 
    sigma = 0.1, 
    nclust = NULL, 
    tau = 0, 
    block.size = 0.05, 
    max.iter.harmony = 10, 
    max.iter.cluster = 20, 
    epsilon.cluster = 1e-5, 
    epsilon.harmony = 1e-4, 
    plot_convergence = FALSE, 
    verbose = TRUE, 
    reference_values = NULL,
    reduction.save = "harmony",
    ...
) {
    
    tryCatch(
        pca_embedding <- Seurat::GetCellEmbeddings(
            object, reduction.type = 'pca'
        ), 
        error = function(e) {
            if (verbose) {
                message('Harmony needs PCA. Trying to run PCA now.')
            }            
            tryCatch(
                object <- Seurat::RunPCA(object, do.print = verbose), 
                error = function(e) {
                    stop('Harmony needs PCA. Tried to run PCA and failed.')
                }
            )
        }
    )    
    if (is.null(dims.use)) {
        dims.use <- seq_len(ncol(pca_embedding))
    } 
    dims_avail <- seq_len(ncol(pca_embedding))
    if (!all(dims.use %in% dims_avail)) {
        stop("trying to use more dimensions than computed with PCA. Rereun 
            PCA with more dimensions or use fewer PCs")
    }
    if (length(dims.use) == 1) {
        stop("only specified one dimension in dims.use")
    }
    metavars_df <- Seurat::FetchData(object, group.by.vars)
    
    harmonyEmbed <- HarmonyMatrix(
        pca_embedding,
        metavars_df,
        group.by.vars, 
        FALSE, 
        0, 
        theta, 
        lambda, 
        sigma, 
        nclust, 
        tau, 
        block.size, 
        max.iter.harmony, 
        max.iter.cluster,
        epsilon.cluster, 
        epsilon.harmony,
        plot_convergence, 
        FALSE, 
        verbose, 
        reference_values)
    
    rownames(harmonyEmbed) <- row.names(pca_embedding)
    colnames(harmonyEmbed) <- paste0("harmony_", seq_len(ncol(harmonyEmbed)))
    
    object <- object %>%
        Seurat::SetDimReduction(reduction.type = "harmony", 
                                slot = "cell.embeddings", 
                                new.data = harmonyEmbed) %>%
        Seurat::SetDimReduction(reduction.type = "harmony", 
                                slot = "key", 
                                new.data = "Harmony") %>%
        Seurat::ProjectDim(reduction.type = "harmony", 
                            replace.dim = TRUE, do.print = FALSE)
    
    return(object)
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
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = 'RNA',
  project.dim = TRUE,  
  weight.by = NULL,
  ...
) {
  if (reduction == 'pca') {
    tryCatch(
      embedding <- Seurat::Embeddings(object, reduction = 'pca'),
      error = function(e) {
        if (verbose) {
          message('Harmony needs PCA. Trying to run PCA now.')
        }
        tryCatch(
          object <- Seurat::RunPCA(
            object, assay = assay.use, verbose = verbose
          ),
          error = function(e) {
            stop('Harmony needs PCA. Tried to run PCA and failed.')
          }
        )
      }
    )
  } else {
    available.dimreduc <- names(methods::slot(object = object, name = 'reductions'))
    if (!(reduction %in% available.dimreduc)) {
      stop("Requested dimension reduction is not present in the Seurat object")
    }
    embedding <- Seurat::Embeddings(object, reduction = reduction)
  }
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(object, group.by.vars)

  ## if using weights, set them up here
  if (is.null(weight.by)) {
    weights <- NULL
  } else {
    if (!weight.by %in% colnames(object@meta.data)) {
        stop('Weighting variable not in meta.data')
    }
    y <- object@meta.data[[weight.by]]
    if (!is(y, 'character') & !is(y, 'factor')) {
        stop('Weighting variable must either encode factor or character type. Numerical not permitted.')
    }
    weights <- as.numeric(((1 / table(y))[y]) * (length(y) / length(unique(y))))
  }
    
  harmonyEmbed <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    FALSE,
    verbose,
    reference_values,
    weights = weights
  )

  rownames(harmonyEmbed) <- row.names(embedding)
  colnames(harmonyEmbed) <- paste0("harmony_", seq_len(ncol(harmonyEmbed)))

  suppressWarnings({
    harmonydata <- Seurat::CreateDimReducObject(
      embeddings = harmonyEmbed,
      assay = assay.use,
      key = "harmony"
    )
  })

  object[[reduction.save]] <- harmonydata
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = "harmony",
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
    lambda = NULL, 
    sigma = 0.1, 
    nclust = NULL, 
    tau = 0, 
    block.size = 0.05, 
    max.iter.harmony = 10, 
    max.iter.cluster = 20, 
    epsilon.cluster = 1e-5, 
    epsilon.harmony = 1e-4, 
    plot_convergence = FALSE, 
    verbose = TRUE, 
    reference_values = NULL,
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
        stop("trying to use more dimensions than computed with PCA. Rereun 
            PCA with more dimensions or use fewer PCs")
    }
    
    metavars_df <- SingleCellExperiment::colData(object)
    if (!all(group.by.vars %in% colnames(metavars_df))) {
        stop('Trying to integrate over variables missing in colData')
    }
    
    harmonyEmbed <- HarmonyMatrix(
        pca_embedding,
        metavars_df,
        group.by.vars, 
        FALSE, 
        0, 
        theta, 
        lambda, 
        sigma, 
        nclust, 
        tau, 
        block.size, 
        max.iter.harmony, 
        max.iter.cluster,
        epsilon.cluster, 
        epsilon.harmony,
        plot_convergence, 
        FALSE, 
        verbose, 
        reference_values
    )
    
    rownames(harmonyEmbed) <- row.names(metavars_df)
    colnames(harmonyEmbed) <- paste0("harmony_", seq_len(ncol(harmonyEmbed)))
    SingleCellExperiment::reducedDim(object, reduction.save) <- harmonyEmbed
    
    return(object)
}

context('Test Harmony interfaces to external packages')

library(harmony)

check_seurat2 <- function() {
    if (!requireNamespace('Seurat', quietly = TRUE)) {
        skip('Seurat V2 not available')
    } else {
        pkg_version <- packageVersion('Seurat')
        if (pkg_version >= "3.0" | pkg_version < "2.0") {
            skip('Seurat V2 not available')
        }   
    }
}

test_that('Seurat V2 interface works', {
    check_seurat2()
    data(cell_lines_small_seurat_v2)
    obj <- cell_lines_small_seurat_v2
    library(Seurat)
    obj <- RunHarmony(obj, "dataset", theta = 1, nclust = 50, lambda = .1,
                       max.iter.cluster = 5, max.iter.harmony = 2, 
                       verbose = FALSE)
    
    expect_true(
      tryCatch({
        V <- Seurat::GetCellEmbeddings(obj, 'harmony')
        return(TRUE)
      }, error = function(e) {
        return(FALSE)
      })
    )

    V <- Seurat::GetCellEmbeddings(obj, 'harmony')
    V_pca <- Seurat::GetCellEmbeddings(obj, 'pca')
    expect_equal(sum(is.na(V)), 0)
    expect_equal(dim(V), dim(V_pca))

})

check_seurat3 <- function() {
    if (!requireNamespace('Seurat', quietly = TRUE)) {
        skip('Seurat V3 not available')
    } else {
        pkg_version <- packageVersion('Seurat')
        if (pkg_version >= "4.0" | pkg_version < "3.0") {
            skip('Seurat V3 not available')
        }   
    }
}

test_that('Seurat V3 interface works', {
    check_seurat3()
    data(cell_lines_small_seurat_v3)
    obj <- cell_lines_small_seurat_v3
    library(Seurat)
    obj <- RunHarmony(obj, "dataset", theta = 1, nclust = 50, lambda = .1,
                       max.iter.cluster = 5, max.iter.harmony = 2, 
                       verbose = FALSE)
    
    expect_true(
      tryCatch({
        V <- Seurat::Embeddings(obj, 'harmony')
        return(TRUE)
      }, error = function(e) {
        return(FALSE)
      })
    )

    V <- Seurat::Embeddings(obj, 'harmony')
    V_pca <- Seurat::Embeddings(obj, 'pca')
    expect_equal(sum(is.na(V)), 0)
    expect_equal(dim(V), dim(V_pca))
    
})



test_that('SingleCellExperiment interface works', {
    if (!requireNamespace('SingleCellExperiment', quietly = TRUE)) {
        skip('SingleCellExperiment not available')
    } 

    data(cell_lines_small_sce)
    obj <- cell_lines_small_sce
    library(SingleCellExperiment)
    obj <- RunHarmony(obj, "dataset", theta = 1, nclust = 50, lambda = .1,
                       max.iter.cluster = 5, max.iter.harmony = 2, 
                       verbose = FALSE)

    expect_true('HARMONY' %in% SingleCellExperiment::reducedDimNames(obj))
    
    V <- SingleCellExperiment::reducedDim(obj, 'HARMONY')
    V_pca <- SingleCellExperiment::reducedDim(obj, 'PCA')
    expect_equal(sum(is.na(V)), 0)
    expect_equal(dim(V), dim(V_pca))
    
})




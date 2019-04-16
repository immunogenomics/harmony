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
    data(cell_lines_small_seurat)
    obj <- cell_lines_small_seurat
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


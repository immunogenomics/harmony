context('Test Harmony interfaces to external packages')

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
    library(Seurat)
    cell_lines_small_seurat <- RunHarmony(cell_lines_small_seurat, "dataset", theta = 1, nclust = 50, lambda = .1,
                       max.iter.cluster = 5, max.iter.harmony = 2, verbose = FALSE)
    expect_true('harmony' %in% names(cell_lines_small_seurat@dr))
    expect_equal(sum(is.na(cell_lines_small_seurat@dr$harmony@cell.embeddings)), 0)
    expect_equal(dim(cell_lines_small_seurat@dr$harmony@cell.embeddings), dim(cell_lines_small_seurat@dr$pca@cell.embeddings))
})


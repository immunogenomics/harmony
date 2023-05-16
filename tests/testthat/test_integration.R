context('Test main Harmony integration function: HarmonyMatrix')
library(harmony)
data(cell_lines_small)

obj <- HarmonyMatrix(cell_lines_small$scaled_pcs, cell_lines_small$meta_data, 
                     theta = 1, nclust = 50, lambda = .1,
                     max.iter.cluster = 10, max.iter.harmony = 5,
                     'dataset', do_pca = FALSE, return_object = TRUE, 
                     verbose = FALSE)

test_that('dimensions match in Harmony object data structures', {
    expect_equal(dim(obj$Y), c(obj$d, obj$K))
    expect_equal(dim(obj$Z_corr), c(obj$d, obj$N))
    expect_equal(dim(obj$Z_cos), c(obj$d, obj$N))
    expect_equal(dim(obj$R), c(obj$K, obj$N))    
})

test_that('R defines proper probability distributions', {
    expect_gte(min(obj$R), 0)
    expect_lte(max(obj$R), 1)
    expect_equal(colSums(obj$R), rep(1, obj$N))
})

test_that('there are no null values in the corrected embedding', {
    expect_true(all(!is.infinite(obj$Z_corr)))
    expect_true(all(!is.na(obj$Z_corr)))
    expect_true(all(!is.infinite(obj$Z_cos)))
    expect_true(all(!is.na(obj$Z_cos)))
})


test_that('increasing theta decreases chi2 between Cluster and Batch assign', {
    obj0 <- HarmonyMatrix(cell_lines_small$scaled_pcs, 
                          cell_lines_small$meta_data,
                         theta = 0, nclust = 20, lambda = .1,
                         max.iter.cluster = 5, max.iter.harmony = 2,
                         'dataset', do_pca = FALSE, return_object = TRUE,
                         verbose = FALSE)    
    obj1 <- HarmonyMatrix(cell_lines_small$scaled_pcs, 
                          cell_lines_small$meta_data, 
                         theta = 1, nclust = 20, lambda = .1,
                         max.iter.cluster = 5, max.iter.harmony = 2,
                         'dataset', do_pca = FALSE, return_object = TRUE, 
                         verbose = FALSE)    
    
    expect_gt(
        sum(((obj0$O - obj0$E) ^ 2) / obj0$E), 
        sum(((obj1$O - obj1$E) ^ 2) / obj1$E)
    )
})

test_that('error messages work', {
    expect_error(
        HarmonyMatrix(cell_lines_small$scaled_pcs, cell_lines_small$meta_data,
                      do_pca = FALSE)
    )
    expect_error(
        HarmonyMatrix(cell_lines_small$scaled_pcs, cell_lines_small$meta_data,
                      'fake_variable', do_pca = FALSE)
    )
    expect_error(
        HarmonyMatrix(cell_lines_small$scaled_pcs, 
                      head(cell_lines_small$meta_data, -1), 'dataset', 
                      do_pca = FALSE)
    )
})

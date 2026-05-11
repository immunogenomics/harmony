context('Test 2-level variable correction: cell_lines with cell_type and dataset')
library(harmony)
data(cell_lines)

obj <- RunHarmony(
    cell_lines$scaled_pcs, cell_lines$meta_data,
    vars_use = c('cell_type', 'dataset'),
    theta = c(1, 1), nclust = 50, max_iter = 10,
    return_object = TRUE, verbose = FALSE,
    .options = harmony_options(max.iter.cluster = 10)
)

test_that('two-variable run: core dimensions are consistent', {
    expect_equal(dim(obj$Y),        c(obj$d, obj$K))
    expect_equal(dim(obj$getZcorr()), c(obj$d, obj$N))
    expect_equal(dim(obj$getZorig()), c(obj$d, obj$N))
    expect_equal(dim(obj$R),        c(obj$K, obj$N))
})

test_that('two-variable run: O and E columns span all levels of both covariates', {
    n_levels <- length(unique(cell_lines$meta_data$cell_type)) +
                length(unique(cell_lines$meta_data$dataset))
    expect_equal(ncol(obj$O), n_levels)
    expect_equal(ncol(obj$E), n_levels)
})

test_that('two-variable run: R defines proper probability distributions', {
    expect_gte(min(obj$R), 0)
    expect_lte(max(obj$R), 1)
    expect_equal(colSums(obj$R), rep(1, obj$N), tolerance = 1e-5)
})

test_that('two-variable run: corrected embedding has no NA or infinite values', {
    Z_corr <- obj$getZcorr()
    expect_true(all(!is.na(Z_corr)))
    expect_true(all(!is.infinite(Z_corr)))
})

test_that('two-variable run: higher theta reduces batch/cluster association for both covariates', {
    obj_lo <- RunHarmony(
        cell_lines$scaled_pcs, cell_lines$meta_data,
        vars_use = c('cell_type', 'dataset'),
        theta = c(0, 0), nclust = 20, max_iter = 2,
        return_object = TRUE, verbose = FALSE
    )
    obj_hi <- RunHarmony(
        cell_lines$scaled_pcs, cell_lines$meta_data,
        vars_use = c('cell_type', 'dataset'),
        theta = c(2, 2), nclust = 20, max_iter = 2,
        return_object = TRUE, verbose = FALSE
    )
    chi2_lo <- sum(((obj_lo$O - obj_lo$E) ^ 2) / obj_lo$E)
    chi2_hi <- sum(((obj_hi$O - obj_hi$E) ^ 2) / obj_hi$E)
    expect_gt(chi2_lo, chi2_hi)
})

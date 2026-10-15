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

test_that('a covariate with theta zero does not change cluster assignments', {
    data(cell_lines_small)
    options <- harmony_options(
        block.size = 0.2, max.iter.cluster = 1, epsilon.cluster = -Inf
    )

    set.seed(11)
    one_variable <- RunHarmony(
        cell_lines_small$scaled_pcs, cell_lines_small$meta_data,
        vars_use = 'dataset', theta = 1, nclust = 5, sigma = 0.1,
        max_iter = 0, return_object = TRUE, verbose = FALSE,
        .options = options
    )
    set.seed(11)
    two_variables <- RunHarmony(
        cell_lines_small$scaled_pcs, cell_lines_small$meta_data,
        vars_use = c('dataset', 'cell_type'), theta = c(1, 0),
        nclust = 5, sigma = 0.1, max_iter = 0,
        return_object = TRUE, verbose = FALSE, .options = options
    )

    set.seed(29)
    one_variable$cluster_cpp()
    set.seed(29)
    two_variables$cluster_cpp()

    expect_equal(two_variables$R, one_variable$R, tolerance = 1e-6)
})

context('Test main Harmony integration function: RunHarmony')
library(harmony)
data(cell_lines_small)

obj <- RunHarmony(cell_lines_small$scaled_pcs, cell_lines_small$meta_data, 'dataset',
                  theta = 1, nclust = 50, max_iter = 5, return_object = TRUE,
                  verbose = FALSE, .options = harmony_options(max.iter.cluster = 10))

test_that('dimensions match in Harmony object data structures', {
    expect_equal(dim(obj$Y), c(obj$d, obj$K))
    expect_equal(dim(obj$getZcorr()), c(obj$d, obj$N))
    expect_equal(dim(obj$getZorig()), c(obj$d, obj$N))
    expect_equal(dim(obj$R), c(obj$K, obj$N))
})

test_that('R defines proper probability distributions', {
    expect_gte(min(obj$R), 0)
    expect_lte(max(obj$R), 1)
    expect_equal(colSums(obj$R), rep(1, obj$N), tolerance = 1e-5)
})

test_that('there are no null values in the corrected embedding', {
    Z_corr <- obj$getZcorr()
    expect_true(all(!is.infinite(Z_corr)))
    expect_true(all(!is.na(Z_corr)))
})


test_that('increasing theta decreases chi2 between Cluster and Batch assign', {
    obj0 <- RunHarmony(cell_lines_small$scaled_pcs, cell_lines_small$meta_data, 'dataset',
                       theta = 0, nclust = 20, max_iter = 2, return_object = TRUE,
                       verbose = FALSE)
    obj1 <- RunHarmony(cell_lines_small$scaled_pcs, cell_lines_small$meta_data, 'dataset',
                       theta = 1, nclust = 5, max_iter = 2, return_object = TRUE,
                       verbose = FALSE)

    expect_gt(
        sum(((obj0$O - obj0$E) ^ 2) / obj0$E),
        sum(((obj1$O - obj1$E) ^ 2) / obj1$E)
    )
})

test_that('error messages work', {
    expect_error(
        RunHarmony(cell_lines_small$scaled_pcs, cell_lines_small$meta_data, 'fake_variable')
    )

    expect_error(
        RunHarmony(cell_lines_small$scaled_pcs, cell_lines_small$meta_data, 'dataset', lambda = c(1,2))
    )

    expect_error(
        RunHarmony(cell_lines_small$scaled_pcs, head(cell_lines_small$meta_data, -1), 'dataset')
    )

})

test_that('outer convergence requires a small nonnegative decrease', {
    original_objectives <- obj$objective_harmony
    transitions <- list(
        c(100, 99.5),
        c(100, 99),
        c(100, 90),
        c(100, 101),
        c(0, 0),
        c(0, 1),
        c(0, -1),
        c(NaN, 1),
        c(1, NaN),
        c(Inf, 1),
        c(1, Inf),
        c(-Inf, 1),
        c(1, -Inf)
    )
    converged <- vapply(transitions, function(values) {
        obj$objective_harmony <- values
        obj$check_convergence(1)
    }, logical(1))
    obj$objective_harmony <- original_objectives

    expect_identical(
        converged,
        c(TRUE, FALSE, FALSE, FALSE, TRUE, rep(FALSE, 8))
    )
})

test_that('Harmony continues after an outer objective increase', {
    args <- list(
        data_mat = cell_lines_small$scaled_pcs,
        meta_data = cell_lines_small$meta_data,
        vars_use = 'dataset',
        theta = 8,
        sigma = 0.03,
        lambda = 0.1,
        nclust = 5,
        max_iter = 4,
        return_object = TRUE,
        verbose = FALSE,
        .options = harmony_options(
            max.iter.cluster = 10,
            epsilon.harmony = 0.01
        )
    )

    set.seed(69)
    early_stop_run <- do.call(RunHarmony, args)
    args$early_stop <- FALSE
    set.seed(69)
    full_run <- do.call(RunHarmony, args)

    expect_gt(early_stop_run$objective_harmony[4],
              early_stop_run$objective_harmony[3])
    expect_lt(early_stop_run$objective_harmony[5],
              early_stop_run$objective_harmony[4])
    expect_length(early_stop_run$objective_harmony, 5)
    expect_equal(
        early_stop_run$objective_harmony,
        full_run$objective_harmony,
        tolerance = 0
    )
    expect_equal(
        early_stop_run$getZcorr(),
        full_run$getZcorr(),
        tolerance = 0
    )
})

test_that('ordinary one-covariate output matches the pre-change result', {
    set.seed(42)
    ordinary_run <- RunHarmony(
        cell_lines_small$scaled_pcs,
        cell_lines_small$meta_data,
        vars_use = 'dataset',
        theta = 1,
        sigma = 0.1,
        lambda = 1,
        nclust = 5,
        max_iter = 5,
        early_stop = TRUE,
        return_object = TRUE,
        verbose = FALSE,
        .options = harmony_options(
            block.size = 0.05,
            max.iter.cluster = 10,
            epsilon.cluster = 1e-3,
            epsilon.harmony = 1e-2,
            batch.prop.cutoff = 1e-5
        )
    )
    Z_corr <- ordinary_run$getZcorr()
    cells <- c(1, 73, 151, 227, 300)
    z_index <- rbind(
        c(4, 1),
        c(1, 73),
        c(6, 151),
        c(1, 227),
        c(1, 300)
    )

    # Reference from da8de4d, before the optimizer reliability changes.
    expected_objective <- c(
        723.603759765625,
        722.236083984375
    )
    expected_R_max <- c(
        0.99899065494537354,
        0.91011166572570801,
        0.95751070976257324,
        0.75582998991012573,
        0.89015030860900879
    )
    expected_Z <- c(
        -0.0088489148765802383,
        -0.0078113484196364880,
        0.0070061264559626579,
        -0.0132178822532296181,
        -0.0114015899598598480
    )

    expect_length(ordinary_run$objective_harmony, 2)
    expect_equal(
        ordinary_run$objective_harmony,
        expected_objective,
        tolerance = 5e-4
    )
    expect_equal(sum(ordinary_run$R^2), 251.39820952111177,
                 tolerance = 1e-4)
    expect_equal(sum(Z_corr^2), 0.050813484665478796,
                 tolerance = 1e-6)
    expect_equal(
        apply(ordinary_run$R[, cells, drop = FALSE], 2, max),
        expected_R_max,
        tolerance = 2e-6
    )
    expect_equal(
        Z_corr[z_index],
        expected_Z,
        tolerance = 2e-7
    )
})

test_that('small sigma stays finite and invalid assignments fail closed', {
    donor_meta <- data.frame(
        donor = factor(cell_lines_small$meta_data$cell_type)
    )
    set.seed(1)
    small_sigma <- RunHarmony(
        cell_lines_small$scaled_pcs,
        donor_meta,
        vars_use = 'donor',
        theta = 1,
        sigma = 1e-8,
        lambda = 1,
        nclust = 5,
        max_iter = 0,
        return_object = TRUE,
        verbose = FALSE,
        .options = harmony_options(max.iter.cluster = 1)
    )

    distance <- 2 * (1 - crossprod(
        small_sigma$Y,
        small_sigma$getZcorr()
    ))
    legacy_R <- exp(-distance / 1e-8)
    legacy_R <- sweep(legacy_R, 2, colSums(legacy_R), '/')
    expect_false(all(is.finite(legacy_R)))

    invisible(small_sigma$cluster_cpp())
    small_sigma$moe_correct_ridge_cpp()
    state <- c(
        small_sigma$R,
        small_sigma$O,
        small_sigma$E,
        small_sigma$objective_kmeans,
        small_sigma$objective_harmony,
        small_sigma$getZcorr()
    )
    expect_true(all(is.finite(state)))
    expect_equal(
        colSums(small_sigma$R),
        rep(1, small_sigma$N),
        tolerance = 1e-6
    )

    small_sigma$R <- legacy_R
    expect_error(
        small_sigma$moe_correct_ridge_cpp(),
        paste0(
            'nonfinite state during correction input.*',
            'sigma=\\[1e-08, 1e-08\\].*',
            'theta=\\[1, 1\\].*',
            'block\\.size=0.05, N=300, K=5'
        )
    )

    set.seed(1)
    tiny_sigma <- RunHarmony(
        cell_lines_small$scaled_pcs,
        donor_meta,
        vars_use = 'donor',
        theta = 1,
        sigma = 1e-40,
        lambda = 1,
        nclust = 5,
        max_iter = 1,
        early_stop = FALSE,
        return_object = TRUE,
        verbose = FALSE,
        .options = harmony_options(max.iter.cluster = 1)
    )
    tiny_state <- c(
        tiny_sigma$R,
        tiny_sigma$O,
        tiny_sigma$E,
        tiny_sigma$objective_kmeans,
        tiny_sigma$objective_harmony,
        tiny_sigma$getZcorr()
    )
    expect_true(all(is.finite(tiny_state)))
    expect_equal(
        colSums(tiny_sigma$R),
        rep(1, tiny_sigma$N),
        tolerance = 1e-6
    )
})

test_that('stable assignments match moderate logits and ignore constant shifts', {
    set.seed(11)
    assignment_obj <- RunHarmony(
        cell_lines_small$scaled_pcs,
        cell_lines_small$meta_data,
        vars_use = 'dataset',
        theta = 0,
        sigma = 0.3,
        nclust = 5,
        max_iter = 0,
        return_object = TRUE,
        verbose = FALSE,
        .options = harmony_options(max.iter.cluster = 0)
    )

    distance <- 2 * (1 - crossprod(
        assignment_obj$Y,
        assignment_obj$getZcorr()
    ))
    legacy_R <- exp(-distance / 0.3)
    legacy_R <- sweep(legacy_R, 2, colSums(legacy_R), '/')
    expect_equal(assignment_obj$R, legacy_R, tolerance = 1e-5)

    invisible(assignment_obj$cluster_cpp())
    invisible(assignment_obj$cluster_cpp())
    unshifted_R <- assignment_obj$R
    shift <- 15 * assignment_obj$getZcorr()[, 1]
    assignment_obj$Y <- assignment_obj$Y + matrix(
        shift,
        nrow = assignment_obj$d,
        ncol = assignment_obj$K
    )
    invisible(assignment_obj$cluster_cpp())

    expect_true(all(is.finite(assignment_obj$R)))
    expect_equal(assignment_obj$R, unshifted_R, tolerance = 1e-4)
    expect_equal(
        colSums(assignment_obj$R),
        rep(1, assignment_obj$N),
        tolerance = 1e-6
    )
})

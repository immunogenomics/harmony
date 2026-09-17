context('Test 2-level variable correction: cell_lines with cell_type and dataset')
library(harmony)
data(cell_lines)
data(cell_lines_small)

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

test_that('log-space updates preserve the multiplicative assignment formula', {
    balanced_cells <- unlist(lapply(c('jurkat', 't293'), function(cell_type) {
        utils::head(which(cell_lines_small$meta_data$cell_type == cell_type), 24)
    }))
    assignment_data <- t(cell_lines_small$scaled_pcs[
        balanced_cells,
        seq_len(4)
    ])
    assignment_meta <- data.frame(
        donor = factor(rep(c('donor_male', 'donor_female'), each = 24)),
        chemistry = factor(rep(
            rep(c('10x_3prime_v2', '10x_3prime_v3'), each = 12),
            2
        ))
    )
    vars <- c('donor', 'chemistry')
    phi <- Reduce(rbind, lapply(vars, function(var) {
        Matrix::t(Matrix::sparse.model.matrix(
            ~0 + as.factor(assignment_meta[[var]])
        ))
    }))
    B_vec <- vapply(vars, function(var) {
        nlevels(assignment_meta[[var]])
    }, integer(1))
    level_theta <- rep(c(0.4, 0.8), times = B_vec)

    set.seed(11)
    assignment_obj <- methods::new(harmony:::harmony)
    assignment_obj$setup(
        assignment_data,
        phi,
        rep(0.3, 3),
        level_theta,
        -1,
        0.2,
        1L,
        -Inf,
        -Inf,
        3L,
        1,
        B_vec,
        1e-5,
        FALSE
    )
    assignment_obj$init_cluster_cpp()
    base_R <- assignment_obj$R

    count_index <- matrix(seq_len(3 * sum(B_vec)), nrow = 3)
    frozen_E <- 0.5 + count_index / 5
    frozen_O <- 0.25 + ((7 * count_index) %% 11) / 3
    assignment_obj$E <- assignment_obj$E + frozen_E
    assignment_obj$O <- assignment_obj$O + frozen_O
    invisible(assignment_obj$cluster_cpp())

    level_factor <- sweep(
        (2 * frozen_E + 1) / (frozen_O + frozen_E + 1),
        2,
        level_theta,
        '^'
    )
    block_end <- cumsum(B_vec)
    block_start <- c(1, utils::head(block_end, -1) + 1)
    selected_factor <- Map(function(first, last) {
        rows <- seq.int(first, last)
        as.matrix(level_factor[, rows, drop = FALSE] %*%
                  phi[rows, , drop = FALSE])
    }, block_start, block_end)
    expected_R <- base_R * Reduce(`*`, selected_factor)
    expected_R <- unname(sweep(expected_R, 2, colSums(expected_R), '/'))

    expect_equal(assignment_obj$R, expected_R, tolerance = 1e-6)
})

test_that('extreme theta and unbalanced chemistry stay finite', {
    donor_chemistry_meta <- data.frame(
        donor = factor(
            cell_lines_small$meta_data$cell_type,
            levels = c('jurkat', 't293'),
            labels = c('donor_male', 'donor_female')
        ),
        chemistry = factor(c(
            rep('10x_3prime_v2', nrow(cell_lines_small$meta_data) - 1),
            '10x_3prime_v3'
        ))
    )

    set.seed(1)
    extreme_obj <- RunHarmony(
        cell_lines_small$scaled_pcs,
        donor_chemistry_meta,
        vars_use = c('donor', 'chemistry'),
        theta = c(1000, 1000),
        sigma = 0.1,
        lambda = 1,
        nclust = 5,
        max_iter = 1,
        early_stop = FALSE,
        return_object = TRUE,
        verbose = FALSE,
        .options = harmony_options(
            block.size = 0.07,
            max.iter.cluster = 1
        )
    )

    state <- c(
        extreme_obj$R,
        extreme_obj$O,
        extreme_obj$E,
        extreme_obj$objective_kmeans,
        extreme_obj$objective_harmony,
        extreme_obj$getZcorr()
    )
    expect_true(all(is.finite(state)))
    expect_equal(
        colSums(extreme_obj$R),
        rep(1, extreme_obj$N),
        tolerance = 1e-6
    )
})

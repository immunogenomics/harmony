# Standalone benchmark worker. Each process loads exactly one package build.
args <- commandArgs(TRUE)
library(harmony, lib.loc = args[1])
RhpcBLASctl::blas_set_num_threads(1)
data(cell_lines_small)
data(cell_lines)
mode <- args[2]
output <- args[3]

if (mode == 'resource') {
    n <- as.integer(args[4])
    block <- as.numeric(args[5])
    idx <- rep(seq_len(nrow(cell_lines$scaled_pcs)), length.out = n)
    set.seed(42)
    elapsed <- system.time(object <- RunHarmony(
        cell_lines$scaled_pcs[idx, ], cell_lines$meta_data[idx, ], 'dataset',
        theta = 1, sigma = 0.1, lambda = 1, nclust = 100,
        max_iter = 2, early_stop = FALSE, return_object = TRUE,
        verbose = FALSE,
        .options = harmony_options(max.iter.cluster = 5, block.size = block)
    ))[['elapsed']]
    write.csv(data.frame(
        cells = n, block = block, elapsed = elapsed,
        objective_finite = all(is.finite(object$objective_harmony))
    ), output, row.names = FALSE)
} else {
    rows <- list()
    states <- list()
    for (dataset in c('small', 'full')) {
        dat <- if (dataset == 'small') cell_lines_small else cell_lines
        for (seed in c(1, 42, 69)) {
            for (sigma in c(0.1, 0.03, 0.02, 0.01, 1e-8)) {
                for (covariate in c('dataset', 'cell_type')) {
                    key <- paste(dataset, seed, sigma, covariate, sep = '_')
                    set.seed(seed)
                    object <- RunHarmony(
                        dat$scaled_pcs, dat$meta_data, covariate,
                        theta = 1, sigma = sigma, lambda = 1, nclust = 5,
                        max_iter = 10, return_object = TRUE, verbose = FALSE,
                        .options = harmony_options(max.iter.cluster = 10)
                    )
                    state <- list(R = object$R, Z = object$getZcorr(),
                                  objective = object$objective_harmony)
                    rows[[key]] <- data.frame(
                        key, dataset, seed, sigma, covariate,
                        finite = all(is.finite(c(object$R, object$O, object$E,
                            object$objective_kmeans, object$objective_harmony))),
                        corrected_finite = all(is.finite(state$Z)),
                        rounds = length(state$objective) - 1
                    )
                    states[[key]] <- state
                }
            }
        }
    }
    write.csv(do.call(rbind, rows), output, row.names = FALSE)
    saveRDS(states, paste0(output, '.rds'))
}

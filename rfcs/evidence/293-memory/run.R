# Run from the repository root inside an R environment with package dependencies.
# Arguments: scratch directory, upstream library, PR293 library, revised library.
# Requires GNU /usr/bin/time. See session.txt for versions.
# Do not run CPU-heavy jobs alongside the timed comparison.
# Reproduction (shell commands from the repository root):
# task_dir=$(mktemp -d)
# mkdir "$task_dir/upstream" "$task_dir/pr293" "$task_dir/revised"
# mkdir "$task_dir/lib-upstream" "$task_dir/lib-pr293" "$task_dir/lib-revised"
# git archive df19af23ae0639bd6ea2da63898f973f08c85862 | tar -x -C "$task_dir/upstream"
# git archive 6470e20ffa700fbd516d4bdbba72d22cb97518eb | tar -x -C "$task_dir/pr293"
# git archive HEAD | tar -x -C "$task_dir/revised"
# R CMD INSTALL -l "$task_dir/lib-upstream" "$task_dir/upstream"
# R CMD INSTALL -l "$task_dir/lib-pr293" "$task_dir/pr293"
# R CMD INSTALL -l "$task_dir/lib-revised" "$task_dir/revised"
# Rscript rfcs/evidence/293-memory/run.R "$task_dir" \
#   "$task_dir/lib-upstream" "$task_dir/lib-pr293" "$task_dir/lib-revised"
# The runner prints Markdown tables for the numerical results.
# The runner replaces CSV results here and keeps logs and matrices in scratch.
args <- commandArgs(TRUE)
scratch <- normalizePath(args[1])
libraries <- setNames(normalizePath(args[2:4]), c('upstream', 'pr293', 'revised'))
worker <- normalizePath('rfcs/evidence/293-memory/benchmark.R')
destination <- 'rfcs/evidence/293-memory'
rscript <- file.path(R.home('bin'), 'Rscript')
run_numerical <- function() {
    numerical <- list()
    states <- list()
    for (version in names(libraries)) {
        output <- file.path(scratch, paste0(version, '-numerical.csv'))
        status <- system2(rscript, shQuote(c(worker, libraries[[version]],
            'numerical', output)), stdout = paste0(output, '.log'),
            stderr = paste0(output, '.log'))
        stopifnot(status == 0)
        numerical[[version]] <- cbind(version, read.csv(output))
        states[[version]] <- readRDS(paste0(output, '.rds'))
    }
    write.csv(do.call(rbind, numerical), file.path(destination, 'numerical.csv'),
              row.names = FALSE)
    agreement <- do.call(rbind, lapply(names(states$revised), function(key) {
        old <- states$pr293[[key]]
        new <- states$revised[[key]]
        data.frame(key, assignments = max(abs(old$R - new$R)),
                   coordinates = max(abs(old$Z - new$Z)),
                   same_rounds = length(old$objective) == length(new$objective),
                   identical_objective = identical(old$objective, new$objective))
    }))
    write.csv(agreement, file.path(destination, 'agreement.csv'), row.names = FALSE)
    source(file.path(destination, 'tables.R'))
    print_tables(do.call(rbind, numerical), agreement)
}
resources <- list()
for (n in c(10000, 50000, 100000)) {
    for (block in c(0.05, 1)) {
        for (replicate in 1:3) {
            # Alternate process order to limit systematic order effects.
            versions <- if (replicate %% 2) c('pr293', 'revised') else
                c('revised', 'pr293')
            for (version in versions) {
                key <- paste(version, n, block, replicate, sep = '-')
                output <- file.path(scratch, paste0(key, '.csv'))
                rss <- paste0(output, '.rss')
                status <- system2('/usr/bin/time', c('-f', '%M', '-o',
                    shQuote(rss), shQuote(rscript),
                    shQuote(c(worker, libraries[[version]], 'resource',
                              output, n, block))),
                    stdout = paste0(output, '.log'),
                    stderr = paste0(output, '.log'))
                stopifnot(status == 0)
                resources[[key]] <- cbind(version, replicate,
                    read.csv(output), peak_mib = scan(rss, quiet = TRUE) / 1024)
                write.csv(do.call(rbind, resources),
                    file.path(destination, 'resources.csv'), row.names = FALSE)
            }
        }
    }
}
# Load the large result matrices only after timed child processes finish.
# Otherwise Linux child peak RSS can inherit the controller's larger footprint.
run_numerical()
writeLines(sub('[[:space:]]+$', '', capture.output(sessionInfo())),
           file.path(destination, 'session.txt'))

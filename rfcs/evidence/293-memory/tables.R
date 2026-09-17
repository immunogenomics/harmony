print_tables <- function(numerical, agreement) {
    cat('| Sigma | Original failures | Full-matrix failures | Cell-by-cell failures |\n',
        '|---|---:|---:|---:|\n', sep = '')
    for (sigma in c(0.1, 0.03, 0.02, 0.01, 1e-8)) {
        counts <- vapply(c('upstream', 'pr293', 'revised'), function(version) {
            finite <- numerical$finite[numerical$version == version &
                                       numerical$sigma == sigma]
            paste(sum(!finite), 'of', length(finite))
        }, character(1))
        cat('|', format(sigma), '|', paste(counts, collapse = ' | '), '|\n')
    }
    cat('\n| Full matrix versus cell by cell | Exact matches |\n',
        '|---|---:|\n', sep = '')
    matches <- c(Assignments = sum(agreement$assignments == 0),
                 'Corrected coordinates' = sum(agreement$coordinates == 0),
                 'Objective traces' = sum(agreement$identical_objective &
                                          agreement$same_rounds))
    for (label in names(matches)) {
        cat('|', label, '|', matches[[label]], 'of', nrow(agreement), '|\n')
    }
}

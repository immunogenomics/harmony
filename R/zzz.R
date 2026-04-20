#' @importFrom utils packageVersion vignette
.onAttach <- function(libname, pkgname) {
    ver = packageVersion(pkgname)
    vignette_built = nrow(vignette(package = 'harmony')$results) != 0

    if (vignette_built) {
        rebuild_advice = ""
    } else {
        rebuild_advice = "Reinstall with {.code build_vignettes = TRUE}, then "
    }

    vign_emph = cli::combine_ansi_styles("bold", "red2")

    vignette_msg = paste0('{cli::symbol$bullet} {.strong Read the guide}: ',
                          rebuild_advice,
                          'run {vign_emph("vignette()")}')

    cli::cli_inform(paste0("{cli::symbol$bullet} This is Harmony2 version {.strong ", ver, "}"), class = "packageStartupMessage")
    cli::cli_inform(vignette_msg, class = "packageStartupMessage")
    cli::cli_inform("{cli::symbol$bullet} {.strong Get help}: Visit the website at https://korsunskylab.github.io/harmony2/ and report issues on https://github.com/immunogenomics/harmony/issues", class = "packageStartupMessage")
}

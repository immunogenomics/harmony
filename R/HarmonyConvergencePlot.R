#' @title HarmonyConvergencePlot
#'
#' @param harmony_obj
#'
#' @importFrom magrittr "%>%"
#' @importFrom tibble rowid_to_column
#' @importFrom purrr reduce
#' @importFrom utils tail
#' @import ggplot2
#'
#' @return
#' @export
#'
#' @examples
HarmonyConvergencePlot <- function(harmony_obj) {
  ## ignore initial value (random initialization)
  ## break down kmeans objective into rounds
  obj_fxn <- data.frame(
    kmeans_idx = reduce(c,
                        lapply(harmony_obj$kmeans_rounds,
                               function(rounds) {
                                   1:rounds
                                   }
                               )
                        ),
    harmony_idx = reduce(c, lapply(1:length(harmony_obj$kmeans_rounds),
                                   function(i) {
                                       rep(i, harmony_obj$kmeans_rounds[i])
                                       }
                                   )
                         ),
    val = tail(harmony_obj$objective_kmeans, -1)) %>%
    rowid_to_column("idx")
  ##    data.table(obj_fxn)[, (tail(.SD$val, 1) - head(.SD$val, 1)) / head(.SD$val, 1), by = harmony_idx]
  obj_fxn %>%
    ggplot(aes(x = idx,
               y = val,
               col = harmony_idx)) +
    geom_point(shape = 21) +
    labs(y = "Objective Function",
         x = "Iteration Number"
         )
}
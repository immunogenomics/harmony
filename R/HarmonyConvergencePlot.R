#' HarmonyConvergencePlot
#'
#' @param harmonyObj 
#'
#' @importFrom magrittr "%>%"
#' @importFrom tibble rowid_to_column
#' @import ggplot2
#' 
#' @return
#' @export
#'
#' @examples
HarmonyConvergencePlot <- function(harmonyObj) {
    ## ignore initial value (random initialization)
    ## break down kmeans objective into rounds
    obj_fxn <- data.frame(
        kmeans_idx = Reduce(c, 
                            lapply(harmonyObj$kmeans_rounds, 
                                   function(rounds) {
                                       1:rounds
                                       }
                                   )
                            ),
        harmony_idx = Reduce(c, 
                             lapply(1:length(harmonyObj$kmeans_rounds), 
                                    function(i) {
                                        rep(i, 
                                            harmonyObj$kmeans_rounds[i])
                                        }
                                    )
                             ),
        val = tail(harmonyObj$objective_kmeans, -1)
    ) %>%
        rowid_to_column("idx")
##    data.table(obj_fxn)[, (tail(.SD$val, 1) - head(.SD$val, 1)) / head(.SD$val, 1), by = harmony_idx]
    obj_fxn %>% 
        ggplot(aes(idx, 
                   val, 
                   col = harmony_idx)) + 
        geom_point(shape = 21) + 
        labs(y = "Objective Function", 
             x = "Iteration Number")
}

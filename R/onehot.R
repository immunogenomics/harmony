#' onehot
#'
#' @param x 
#'
#' @importFrom tibble rowid_to_column
#' @importFrom dplyr mutate select
#' @importFrom tidyr spread
#'
#' @return
#' @export
#'
#' @examples
onehot <- function(x) {
	data.frame(x) %>%
		rowid_to_column("id") %>% 
		mutate(dummy = 1) %>% 
		spread(x, dummy, fill = 0) %>% 
		select(-id) %>%
		as.matrix()
}

onehot <- function(x) {
	data.frame(x) %>%
		tibble::rowid_to_column("id") %>% 
		dplyr::mutate(dummy = 1) %>% 
		tidyr::spread(x, dummy, fill = 0) %>% 
		dplyr::select(-id) %>%
		as.matrix
}

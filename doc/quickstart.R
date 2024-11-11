## ----eval=FALSE---------------------------------------------------------------
#  install.packages('harmony')

## -----------------------------------------------------------------------------
library(harmony)

## -----------------------------------------------------------------------------
data(cell_lines)
V <- cell_lines$scaled_pcs
meta_data <- cell_lines$meta_data


## ----class.source='fold-hide', fig.width=5, fig.height=3, fig.align="center"----

library(ggplot2)

do_scatter <- function(xy, meta_data, label_name, base_size = 12) {    
    palette_use <- c(`jurkat` = '#810F7C', `t293` = '#D09E2D',`half` = '#006D2C')
    xy <- xy[, 1:2]
    colnames(xy) <- c('X1', 'X2')
    plt_df <- xy %>% data.frame() %>% cbind(meta_data)
    plt <- ggplot(plt_df, aes(X1, X2, col = !!rlang::sym(label_name), fill = !!rlang::sym(label_name))) + 
        theme_test(base_size = base_size) +
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1,
                                                        shape = 16, size = 4))) +
        scale_color_manual(values = palette_use) +
        scale_fill_manual(values = palette_use) +
        theme(plot.title = element_text(hjust = .5)) +
        labs(x = "PC 1", y = "PC 2") +
        theme(legend.position = "none") +
        geom_point(shape = '.')
    
    ## Add labels
    data_labels <- plt_df %>%
        dplyr::group_by(!!rlang::sym(label_name)) %>%
        dplyr::summarise(X1 = mean(X1), X2 = mean(X2)) %>%
        dplyr::ungroup()
    plt + geom_label(data = data_labels, aes(label = !!rlang::sym(label_name)), 
                            color = "white", size = 4)
}
p1 <- do_scatter(V, meta_data, 'dataset') + 
    labs(title = 'Colored by dataset')
p2 <- do_scatter(V, meta_data, 'cell_type') + 
    labs(title = 'Colored by cell type')

cowplot::plot_grid(p1, p2)


## -----------------------------------------------------------------------------
harmony_embeddings <- harmony::RunHarmony(
    V, meta_data, 'dataset', verbose=FALSE
)


## ---- fig.width=5, fig.height=3, fig.align="center"---------------------------
p1 <- do_scatter(harmony_embeddings, meta_data, 'dataset') + 
    labs(title = 'Colored by dataset')
p2 <- do_scatter(harmony_embeddings, meta_data, 'cell_type') + 
    labs(title = 'Colored by cell type')
cowplot::plot_grid(p1, p2, nrow = 1)


## -----------------------------------------------------------------------------
sessionInfo()



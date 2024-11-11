## ---- message=FALSE, warning=FALSE, class.source = 'fold-hide'----------------

## Source required libraries
library(data.table)
library(tidyverse)
library(ggthemes)
library(ggrepel)
library(harmony)
library(patchwork)
library(tidyr)

## Useful util functions

cosine_normalize <- function(X, margin) {
    if (margin == 1) {
        res <- sweep(as.matrix(X), 1, sqrt(rowSums(X ^ 2)), '/')
        row.names(res) <- row.names(X)
        colnames(res) <- colnames(X)        
    } else {
        res <- sweep(as.matrix(X), 2, sqrt(colSums(X ^ 2)), '/')
        row.names(res) <- row.names(X)
        colnames(res) <- colnames(X)
    }
    return(res)
}

onehot <- function(vals) {
    t(model.matrix(~0 + as.factor(vals)))
}


colors_use <- c(`jurkat` = rgb(129, 15, 124, maxColorValue=255),
                `t293` = rgb(208, 158, 45, maxColorValue=255),
                `half` = rgb(0, 109, 44, maxColorValue=255))


do_scatter <- function(umap_use, meta_data, label_name, no_guides = TRUE, do_labels = TRUE, nice_names, 
                       palette_use = colors_use,
                       pt_size = 4, point_size = .5, base_size = 10, do_points = TRUE, do_density = FALSE, h = 4, w = 8) {
    umap_use <- umap_use[, 1:2]
    colnames(umap_use) <- c('X1', 'X2')
    plt_df <- umap_use %>% data.frame() %>% 
        cbind(meta_data) %>% 
        dplyr::sample_frac(1L) 
    plt_df$given_name <- plt_df[[label_name]]
    
    if (!missing(nice_names)) {
        plt_df %<>%
            dplyr::inner_join(nice_names, by = "given_name") %>% 
            subset(nice_name != "" & !is.na(nice_name))

        plt_df[[label_name]] <- plt_df$nice_name        
    }
        
    plt <- plt_df %>% 
        ggplot(aes(X1, X2, colour = .data[[label_name]], fill = .data[[label_name]])) + 
        theme_tufte(base_size = base_size) + 
        theme(panel.background = element_rect(fill = NA, color = "black")) + 
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 4)), alpha = FALSE) +
        scale_color_manual(values = palette_use) + 
        scale_fill_manual(values = palette_use) +    
        theme(plot.title = element_text(hjust = .5)) + 
        labs(x = "UMAP 1", y = "UMAP 2") 
    
    if (do_points) 
        plt <- plt + geom_point(size = 0.2)
    if (do_density) 
        plt <- plt + geom_density_2d()    
        

    if (no_guides)
        plt <- plt + guides("none")
    
    if (do_labels) 
        plt <- plt + geom_label_repel(data = data.table(plt_df)[, .(X1 = mean(X1), X2 = mean(X2)), by = label_name], label.size = NA,
                                      aes(label = .data[[label_name]]), color = "white", size = pt_size, alpha = 1, segment.size = 0) + 
        guides(col = FALSE, fill = FALSE)
    return(plt)
}


## -----------------------------------------------------------------------------
data(cell_lines)
V <- cell_lines$scaled_pcs
V_cos <- cosine_normalize(V, 1)
meta_data <- cell_lines$meta_data

## ---- warning=FALSE, fig.width=5, fig.height=3, fig.align="center"------------
do_scatter(V, meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Colored by dataset', x = 'PC1', y = 'PC2') +
do_scatter(V, meta_data, 'cell_type', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Colored by cell type', x = 'PC1', y = 'PC2') +
NULL

## -----------------------------------------------------------------------------

set.seed(1)
harmonyObj <- harmony::RunHarmony(
    data_mat = V, ## PCA embedding matrix of cells
    meta_data = meta_data, ## dataframe with cell labels
    theta = 1, ## cluster diversity enforcement
    vars_use = 'dataset', ## variable to integrate out
    nclust = 5, ## number of clusters in Harmony model
    max_iter = 0, ## stop after initialization
    return_object = TRUE ## return the full Harmony model object
)



## ---- fig.width=5, fig.height=3, fig.align="center"---------------------------
do_scatter(t(harmonyObj$Z_orig), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Z_orig', subtitle = 'Euclidean distance', x = 'PC1', y = 'PC2') +
do_scatter(t(harmonyObj$Z_cos), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Z_cos', subtitle = 'Induced Cosine distance', x = 'PC1', y = 'PC2')


## ---- fig.width=8, fig.height=3, out.width="100%"-----------------------------

harmonyObj$Z_cos %>% t %>% data.frame() %>% 
    cbind(meta_data) %>% 
    tidyr::gather(key, val, X1:X20) %>% 
    ggplot(aes(reorder(gsub('X', 'PC', key), as.integer(gsub('X', '', key))), val)) + 
        geom_boxplot(aes(color = dataset)) + 
        scale_color_manual(values = colors_use) + 
        labs(x = 'PC number', y = 'PC embedding value', title = 'Z_cos (unit scaled PCA embeddings) for all 20 PCs') + 
        theme_tufte(base_size = 10) + geom_rangeframe() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ---- fig.width=4, fig.height=3, fig.align="center"---------------------------

cluster_centroids <- harmonyObj$Y

do_scatter(t(harmonyObj$Z_cos), meta_data, 'dataset', no_guides = FALSE, do_labels = FALSE) + 
    labs(title = 'Initial kmeans cluster centroids', subtitle = '', x = 'PC1', y = 'PC2') +
    geom_point(
        data = data.frame(t(cluster_centroids)), 
        color = 'black', fill = 'black', alpha = .8,
        shape = 21, size = 6
    ) +
NULL


## -----------------------------------------------------------------------------
cluster_assignment_matrix <- harmonyObj$R


## ---- fig.height=5, fig.width=5-----------------------------------------------
t(harmonyObj$Z_cos) %>% data.frame() %>%
    cbind(meta_data) %>% 
    tibble::rowid_to_column('id') %>% 
    dplyr::inner_join(
        cluster_assignment_matrix %>% t() %>% data.table() %>% 
            tibble::rowid_to_column('id') %>%
            tidyr::gather(cluster, r, -id) %>% 
            dplyr::mutate(cluster = gsub('V', 'Cluster ', cluster)), 
        by = 'id'
    ) %>% 
    dplyr::sample_frac(1L) %>% 
    ggplot(aes(X1, X2, color = r)) + 
        geom_point(size=0.2) + 
        theme_tufte(base_size = 10) + theme(panel.background = element_rect()) + 
        facet_grid(cluster ~ dataset) + 
        scale_color_gradient(low = 'lightgrey', breaks = seq(0, 1, .1)) + 
        labs(x = 'Scaled PC1', y = 'Scaled PC2', title = 'Initial probabilistic cluster assignments')

## -----------------------------------------------------------------------------
observed_counts <- with(harmonyObj, R %*% t(as.matrix(Phi)))
round(observed_counts)



## -----------------------------------------------------------------------------
## observed counts
round(harmonyObj$O)

## observed counts
round(harmonyObj$E)


## -----------------------------------------------------------------------------
phi_celltype <- onehot(meta_data$cell_type) 
observed_cell_counts <- harmonyObj$R %*% t(phi_celltype)
round(observed_cell_counts)


## -----------------------------------------------------------------------------
harmonyObj$max_iter_kmeans

## -----------------------------------------------------------------------------
## we can specify how many rounds of clustering to do
harmonyObj$max_iter_kmeans <- 10
harmonyObj$cluster_cpp()

## -----------------------------------------------------------------------------
round(harmonyObj$O)

## ---- fig.height=5, fig.width=5-----------------------------------------------
new_cluster_assignment_matrix <- harmonyObj$R

t(harmonyObj$Z_cos) %>% data.frame() %>%
    cbind(meta_data) %>% 
    tibble::rowid_to_column('id') %>% 
    dplyr::inner_join(
        new_cluster_assignment_matrix %>% t() %>% data.table() %>% 
            tibble::rowid_to_column('id') %>%
            tidyr::gather(cluster, r, -id) %>% 
            dplyr::mutate(cluster = gsub('V', 'Cluster ', cluster)), 
        by = 'id'
    ) %>% 
    dplyr::sample_frac(1L) %>% 
    ggplot(aes(X1, X2, color = r)) + 
        geom_point(shape = '.') + 
        theme_tufte(base_size = 10) + theme(panel.background = element_rect()) + 
        facet_grid(cluster ~ dataset) + 
        scale_color_gradient(low = 'lightgrey', breaks = seq(0, 1, .1)) + 
        labs(x = 'Scaled PC1', y = 'Scaled PC2', title = 'New probabilistic cluster assignments')

## -----------------------------------------------------------------------------
phi_celltype <- onehot(meta_data$cell_type)
observed_cell_counts <- harmonyObj$R %*% t(phi_celltype)
round(observed_cell_counts)

## -----------------------------------------------------------------------------
round(apply(prop.table(observed_cell_counts, 1), 1, min) * 100, 3)

## -----------------------------------------------------------------------------

with(harmonyObj, {
    distance_matrix <- 2 * (1 - t(Y) %*% Z_cos)
    distance_score <- exp(-distance_matrix / as.numeric(sigma))
    diversity_score <- sweep(E / O, 2, theta, '/') %*% as.matrix(Phi)
    ## new assignments are based on distance and diversity
    R_new <- distance_score * diversity_score  
    ## normalize R so each cell sums to 1
    R_new <- prop.table(R_new, 2)    
})


## -----------------------------------------------------------------------------
## with theta = 0
with(harmonyObj, {
    (E / O) ^ 0
})

## -----------------------------------------------------------------------------
## with theta = 1
with(harmonyObj, {
    round((E / O) ^ 1, 2)
})


## -----------------------------------------------------------------------------
## as theta approach infinity
with(harmonyObj, {
    round((E / O) ^ 1e6, 2)
})


## -----------------------------------------------------------------------------
Y_unscaled <- with(harmonyObj, Z_cos %*% t(R))

## -----------------------------------------------------------------------------
Y_new <- cosine_normalize(Y_unscaled, 2)

## -----------------------------------------------------------------------------
harmonyObj$moe_correct_ridge_cpp()

## ---- fig.width=5, fig.height=3, fig.align="center"---------------------------

do_scatter(cosine_normalize(t(harmonyObj$Z_orig), 1), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Z_cos before MoE', x = 'PC1', y = 'PC2') +
do_scatter(t(harmonyObj$Z_cos), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Z_cos after MoE', x = 'PC1', y = 'PC2')

## ---- fig.width=8, fig.height=3, fig.align="center", out.width="100%"---------

do_scatter(t(harmonyObj$Z_orig), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Z_orig', subtitle = 'Original PCA embeddings', x = 'PC1', y = 'PC2') +
do_scatter(t(harmonyObj$Z_corr), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Z_corr', subtitle = '= Z_orig - correction_factors', x = 'PC1', y = 'PC2') +
do_scatter(t(harmonyObj$Z_cos), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Z_cos', subtitle = '= Unit_scaled(Z_corr)', x = 'Scaled PC1', y = 'Scaled PC2') +
NULL

## ---- fig.width=5, fig.height=3, fig.align="center"---------------------------

plt <- data.table(PC1_After = harmonyObj$Z_corr[1, ], PC1_Before = harmonyObj$Z_orig[1, ]) %>% 
    cbind(meta_data) %>% 
    dplyr::sample_frac(1L) %>% 
    ggplot(aes(PC1_Before, PC1_After)) + 
        geom_abline(slope = 1, intercept = 0) + 
        theme_tufte(base_size = 10) + geom_rangeframe() + 
        scale_color_tableau() + 
        guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16, size = 4))) + 
        NULL

plt + geom_point(shape = '.', aes(color = dataset)) + 
        labs(x = 'PC1 before correction', y = 'PC1 after correction', 
             title = 'PC1 correction for each cell', subtitle = 'Colored by Dataset') + 
plt + geom_point(shape = '.', aes(color = cell_type)) + 
        labs(x = 'PC1 before correction', y = 'PC1 after correction', 
             title = 'PC1 correction for each cell', subtitle = 'Colored by Cell Type') + 
NULL


## ---- echo=TRUE---------------------------------------------------------------

W <- list()
## Convert sparse data structures to dense matrix
Phi.moe <- as.matrix(harmonyObj$Phi_moe)
lambda <- diag(c(harmonyObj$lambda))
## Get beta coeeficients for all the clusters
for (k in 1:harmonyObj$K) {
    W[[k]] <- solve(Phi.moe %*% diag(harmonyObj$R[k, ]) %*% t(Phi.moe) + lambda) %*% (Phi.moe %*% diag(harmonyObj$R[k, ])) %*% t(harmonyObj$Z_orig)
}



## ---- fig.width=5, fig.height=5-----------------------------------------------

cluster_assignment_matrix <- harmonyObj$R

t(harmonyObj$Z_orig) %>% data.frame() %>%
    cbind(meta_data) %>% 
    tibble::rowid_to_column('id') %>% 
    dplyr::inner_join(
        cluster_assignment_matrix %>% t() %>% data.table() %>% 
            tibble::rowid_to_column('id') %>%
            tidyr::gather(cluster, r, -id) %>% 
            dplyr::mutate(cluster = gsub('V', 'Cluster ', cluster)), 
        by = 'id'
    ) %>% 
    dplyr::sample_frac(1L) %>% 
    ggplot(aes(X1, X2, color = r)) + 
        geom_point(shape = 0.2) + 
        theme_tufte(base_size = 10) + theme(panel.background = element_rect()) + 
        facet_grid(cluster ~ dataset) + 
        scale_color_gradient(low = 'grey', breaks = seq(0, 1, .2)) + 
        labs(x = 'PC1', y = 'PC2', title = 'Cluster assigned in original PCA space (Z_orig)')


## -----------------------------------------------------------------------------
plt_list <- lapply(1:harmonyObj$K, function(k) {
    plt_df <- W[[k]] %>% data.frame() %>% 
        dplyr::select(X1, X2)
    ## Append n
    plt_df <- plt_df %>% 
        cbind(
            data.frame(t(matrix(unlist(c(c(0, 0), rep(plt_df[1, ], 3))), nrow = 2))) %>% 
                dplyr::rename(x0 = X1, y0 = X2) 
        ) %>%
        cbind(type = c('intercept', unique(meta_data$dataset)))
    plt <- plt_df %>% 
        ggplot() + 
            geom_point(aes(X1, X2),
                       data = t(harmonyObj$Z_orig) %>% data.frame(),
                       size = 0.5,
                       color = 'grey'
            ) + 
            geom_segment(aes(x = x0, y = y0, xend = X1 + x0, yend = X2 + y0, color = type), linewidth=1) + 
            scale_color_manual(values = c('intercept' = 'black', colors_use)) + 
            theme_tufte(base_size = 10) + theme(panel.background = element_rect()) + 
            labs(x = 'PC 1', y = 'PC 2', title = sprintf('Cluster %d', k))
    plt <- plt + guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16)))    
    # if (k == harmonyObj$K) {
    # } else {
    #     plt <- plt + guides(color = FALSE)
    # }
    plt
})



## ---- fig.height=6, fig.width=6-----------------------------------------------
Reduce(`+`, plt_list) + 
  patchwork::plot_annotation(title = 'Mixture of experts beta terms before correction (Z_orig)') + 
  plot_layout(ncol = 2)

## ---- fig.width=4, fig.height=3, fig.align="center"---------------------------

plt_list <- lapply(1:harmonyObj$K, function(k) {
    plt_df <- W[[k]] %>% data.frame() %>% 
        dplyr::select(X1, X2)

    plt_df <- plt_df %>% 
        cbind(
            data.frame(t(matrix(unlist(c(c(0, 0), rep(plt_df[1, ], 3))), nrow = 2))) %>% 
                dplyr::rename(x0 = X1, y0 = X2) 
        ) %>%
        cbind(type = c('intercept', unique(meta_data$dataset))) 

    plt <- plt_df %>% 
        ggplot() + 
            geom_point(aes(X1, X2),
                data = t(harmonyObj$Z_corr) %>% data.frame(),
                shape = '.', 
                color = 'grey'
            ) + 
            geom_segment(aes(x = x0, y = y0, xend = X1 + x0, yend = X2 + y0, color = type), linewidth=1) + 
            scale_color_manual(values = c('intercept' = 'black', colors_use)) + 
            theme_tufte(base_size = 10) + theme(panel.background = element_rect()) + 
            labs(x = 'PC 1', y = 'PC 2', title = sprintf('Cluster %d', k))
    plt <- plt + guides(color = guide_legend(override.aes = list(stroke = 1, alpha = 1, shape = 16)))
    plt
})



## ---- fig.height=6, fig.width=6-----------------------------------------------
Reduce(`+`, plt_list) + 
  patchwork::plot_annotation(title = 'Mixture of experts beta terms after correction (Z_corr)') + 
  plot_layout(ncol = 2)

## ---- echo=TRUE---------------------------------------------------------------

Z_i <- harmonyObj$Z_orig[, 5]
Z_i_pred <- Reduce(`+`, lapply(1:harmonyObj$K, function(k) {
    W[[k]] * harmonyObj$Phi_moe[, 5] * harmonyObj$R[k, 5]
})) %>% colSums



## ---- fig.width=4, fig.height=3, fig.align="center"---------------------------
data.table(obs = Z_i, pred = Z_i_pred) %>% 
    tibble::rowid_to_column('PC') %>% 
    ggplot(aes(obs, pred)) + 
        geom_point(shape = 21) + 
        geom_label_repel(aes(label = PC)) + 
        geom_abline(slope = 1, intercept = 0) + 
        theme_tufte() + geom_rangeframe() + 
        labs(x = 'Observed PC score', 'Predicted PC score', title = 'Observed and predicted values of PC scores\nfor cell 5') + 
        NULL        

## -----------------------------------------------------------------------------
delta <- Reduce(`+`, lapply(1:harmonyObj$K, function(k) {
    W[[k]][2:4, ] * harmonyObj$Phi[, 5] * harmonyObj$R[k, 5]
})) %>% colSums

Z_corrected <- harmonyObj$Z_orig[, 5] - delta


## ---- fig.width=3, fig.height=3, fig.align="center"---------------------------


harmonyObj$Z_orig %>% t %>% data.frame() %>% 
    ggplot(aes(X1, X2)) + 
        geom_point(shape = '.') + 
        geom_point(
            data = data.frame(t(harmonyObj$Z_orig[, 5, drop = FALSE])), 
            color = 'red'
        ) + 
        geom_segment(
            data = data.table(x0 = harmonyObj$Z_orig[1, 5], 
                              y0 = harmonyObj$Z_orig[2, 5], 
                              x1 = Z_corrected[1],
                              y1 = Z_corrected[2]), 
            aes(x = x0, y = y0, xend = x1, yend = y1),
            linewidth = 1,
            color = 'red', 
            arrow = arrow(length = unit(0.05, "npc"), type = 'closed')            
        ) + 
        theme_tufte(base_size = 10) + geom_rangeframe() + 
        labs(x = 'PC1', y = 'PC2', title = 'Correction of cell #5')


## -----------------------------------------------------------------------------

harmonyObj <- RunHarmony(
    data_mat = V, ## PCA embedding matrix of cells
    meta_data = meta_data, ## dataframe with cell labels
    theta = 1, ## cluster diversity enforcement
    vars_use = 'dataset', ## (list of) variable(s) we'd like to Harmonize out
    nclust = 50, ## number of clusters in Harmony model
    max_iter = 0, ## don't actually run Harmony, stop after initialization
    return_object = TRUE ## return the full Harmony model object, not just the corrected PCA matrix
)


## ---- message=FALSE, fig.width=5, fig.height=3, fig.align="center"------------

i <- 0

do_scatter(t(harmonyObj$Z_cos), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = sprintf('Round %d', i), subtitle = 'Colored by dataset', x = 'Scaled PC1', y = 'Scaled PC2') +
do_scatter(t(harmonyObj$Z_cos), meta_data, 'cell_type', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = sprintf('Round %d', i), subtitle = 'Colored by cell type', x = 'Scaled PC1', y = 'Scaled PC2') +
NULL

## ---- fig.width=5, fig.height=3, fig.align="center", message=FALSE------------

for (i in 1:2) {
    harmony:::harmonize(harmonyObj, 1)
    plt <- do_scatter(t(harmonyObj$Z_cos), meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
        labs(title = sprintf('Round %d', i), subtitle = 'Colored by dataset', x = 'Scaled PC1', y = 'Scaled PC2') +
    do_scatter(t(harmonyObj$Z_cos), meta_data, 'cell_type', no_guides = TRUE, do_labels = TRUE) + 
        labs(title = sprintf('Round %d', i), subtitle = 'Colored by cell type', x = 'Scaled PC1', y = 'Scaled PC2') +
    NULL
    plot(plt)
}
    

## -----------------------------------------------------------------------------
sessionInfo()


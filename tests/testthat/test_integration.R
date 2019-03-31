context('Test main Harmony integration function: HarmonyMatrix')

obj <- HarmonyMatrix(pbmc.small$pca_embeddings, pbmc.small$meta_data, 
                     theta = 1, nclust = 50, 
                     max.iter.cluster = 10, max.iter.harmony = 5,
                     'stim', do_pca = FALSE, return_object = TRUE, verbose = FALSE)    

test_that('dimensions match in Harmony object data structures', {
    expect_equal(dim(obj$Y), c(obj$d, obj$K))
    expect_equal(dim(obj$Z_corr), c(obj$d, obj$N))
    expect_equal(dim(obj$Z_cos), c(obj$d, obj$N))
    expect_equal(dim(obj$R), c(obj$K, obj$N))    
})

test_that('R defines proper probability distributions', {
    expect_gte(min(obj$R), 0)
    expect_lte(max(obj$R), 1)
    expect_equal(colSums(obj$R), rep(1, obj$N))    
})

test_that('there are no null values in the corrected embedding', {
    expect_true(all(!is.infinite(obj$Z_corr)))
    expect_true(all(!is.na(obj$Z_corr)))
    expect_true(all(!is.infinite(obj$Z_cos)))
    expect_true(all(!is.na(obj$Z_cos)))
})


test_that('increasing theta decreases chi-squared between Cluster and Batch assignment', {
    obj0 <- HarmonyMatrix(pbmc.small$pca_embeddings, pbmc.small$meta_data, 
                         theta = 0, nclust = 20, 
                         max.iter.cluster = 5, max.iter.harmony = 2,
                         'stim', do_pca = FALSE, return_object = TRUE, verbose = FALSE)    
    obj1 <- HarmonyMatrix(pbmc.small$pca_embeddings, pbmc.small$meta_data, 
                         theta = 1, nclust = 20, 
                         max.iter.cluster = 5, max.iter.harmony = 2,
                         'stim', do_pca = FALSE, return_object = TRUE, verbose = FALSE)    
    
    expect_gt(
        sum(((obj0$O - obj0$E) ^ 2) / obj0$E), 
        sum(((obj1$O - obj1$E) ^ 2) / obj1$E)
    )
})

test_that('error messages work', {
    expect_error(
        HarmonyMatrix(pbmc.small$pca_embeddings, pbmc.small$meta_data, do_pca = FALSE)
    )
    expect_error(
        HarmonyMatrix(pbmc.small$pca_embeddings, pbmc.small$meta_data, 'fake_variable', do_pca = FALSE)
    )
    expect_error(
        HarmonyMatrix(pbmc.small$pca_embeddings, head(pbmc.small$meta_data, -1), 'stim', do_pca = FALSE)
    )
})
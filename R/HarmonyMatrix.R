#' @title HarmonyMatrix
#'
#' @description 
#' 
#' @param pc_mat 
#' @param batch_labels 
#' @param batch_labels2 
#' @param theta 
#' @param theta2 
#' @param sigma 
#' @param alpha strength of dirichlet prior for DKL penalty
#' @param nclust 
#' @param tau desired cells/cluster
#' @param block.size 
#' @param max.iter.harmony 
#' @param max.iter.cluster 
#' @param epsilon.cluster 
#' @param epsilon.harmony 
#' @param burn.in.time 
#' @param plot_convergence window size for kmeans convergenc
#'
#' @return
#' @export
#'
#' @examples
HarmonyMatrix <- function(pc_mat, 
                          batch_labels, 
                          batch_labels2 = NULL, 
                          theta = 1, 
                          theta2 = 1, 
                          sigma = 0.1, 
                          alpha = .1, 
                          nclust = 100, 
                          tau = 0, 
                          block.size = 0.05, 
                          max.iter.harmony = 10, 
                          max.iter.cluster = 200, 
                          epsilon.cluster = 1e-5, 
                          epsilon.harmony = 1e-4, 
                          burn.in.time = 10, 
                          plot_convergence = FALSE) {
    
    cells_as_cols <- TRUE
    if (length(batch_labels) != ncol(pc_mat)) {
        if (length(batch_labels) == nrow(pc_mat)) {
            pc_mat <- t(pc_mat)
            cells_as_cols <- FALSE
        } else {
            stop("ERROR: Number of labels do not correspond to number of samples in PC matrix.")
        }
    }    

    ## RUN HARMONY
    harmonyObj <- new(harmony, 0) ## 0 is a dummy variable - will change later
    batch_mat <- t(onehot(batch_labels))
    harmonyObj$setup(
        pc_mat, ## Z
        batch_mat, ## Phi
        sigma, ## sigma
        theta, ## theta
        max.iter.cluster, ## max.iter
        epsilon.cluster, ## kmeans converge.thresh
        epsilon.harmony, ## harmony epsilon
        TRUE, ## correct Z_orig only
        alpha, ## EXPERIMENTAL: alpha, strength of dirichlet prior for DKL penalty
        nclust, ## K
        tau, ## tau (desired cells/cluster)
        block.size, ## model$block.size
        rep(1, length(batch_labels)), ## EXPERIMENTAL FEATURE: each cell gets its own weight
        FALSE, ## do linear correction on Z_cos?
        rep(1, nrow(batch_mat)), ## EXPERIMENTAL FEATURE: only correct certain batches
        burn.in.time ## window size for kmeans convergence
    )
  
    ## OPTIONAL: 2nd level of batch defined
    if (!is.null(batch_labels2)) {
      harmonyObj$setup_batch2(t(onehot(batch_labels2)), theta2, tau)
    }
  
    harmonyObj$harmonize(max.iter.harmony)
    if (plot_convergence) plot(HarmonyConvergencePlot(harmonyObj))
    
    res <- as.matrix(harmonyObj$Z_corr)
    row.names(res) <- row.names(pc_mat)
    colnames(res) <- colnames(pc_mat)
    if (!cells_as_cols) 
        res <- t(res)
    return(res)
}

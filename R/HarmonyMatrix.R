#' @title HarmonyMatrix
#'
#' @description
#'
#' @param pc_mat
#' @param meta_data
#' @param batch_labels
#' @param theta
#' @param lambda
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
#' @importFrom graphics plot
#' @importFrom methods new
#' @return
#' @export
#'
#' @examples
HarmonyMatrix <- function(pc_mat,
                          meta_data,
                          vars_use,
                          theta = NULL,
                          lambda = NULL,
                          sigma = 0.1,
                          alpha = 0.1,
                          nclust = 100,
                          tau = 0,
                          block.size = 0.05,
                          max.iter.harmony = 10,
                          max.iter.cluster = 200,
                          epsilon.cluster = 1e-5,
                          epsilon.harmony = 1e-4,
                          burn.in.time = 10,
                          plot_convergence = FALSE,
                          return_object = FALSE) {

  ## TODO: check for
  ##    partially observed batch variables (WARNING)
  ##    batch variables with only 1 level (WARNING)
  ##    if lambda given, check correct length
  ##    if theta given, check correct length
  ##    very small batch size and tau=0: suggest tau>0
  ##    is PCA correct?

  N <- nrow(meta_data)
  cells_as_cols <- TRUE
  if (ncol(pc_mat) != N) {
    if (nrow(pc_mat) == N) {
      pc_mat <- t(pc_mat)
      cells_as_cols <- FALSE
    } else {
      stop("ERROR: Number of labels do not correspond to number of samples in PC matrix.")
    }
  }

  if (is.null(theta)) {
    theta <- rep(2, length(vars_use))
  }
  if (is.null(lambda)) {
    lambda <- rep(1, length(vars_use))
  }

  ## Pre-compute some useful statistics
  phi <- Reduce(rbind, lapply(vars_use, function(var_use) {
    t(onehot(meta_data[[var_use]]))
  }))
  N_b <- rowSums(phi)
  Pr_b <- N_b / N
  B_vec <- Reduce(c, lapply(vars_use, function(var_use) {
    length(unique(meta_data[[var_use]]))
  }))
  theta <- Reduce(c, lapply(1:length(B_vec), function(b) rep(theta[b], B_vec[b])))
  theta <- theta * (1 - exp(-(N_b / (nclust * tau))^2))

  lambda <- Reduce(c, lapply(1:length(B_vec), function(b) rep(lambda[b], B_vec[b])))
  lambda_mat <- diag(c(0, lambda))

  ## RUN HARMONY
  harmonyObj <- new(harmony, 0) ## 0 is a dummy variable - will change later
  harmonyObj$setup(
    pc_mat, ## Z
    phi, ## Phi
    Pr_b,
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
    #        rep(1, N), ## EXPERIMENTAL FEATURE: each cell gets its own weight
    FALSE, ## do linear correction on Z_cos?
    burn.in.time, ## window size for kmeans convergence
    lambda_mat
  )


  harmonyObj$harmonize(max.iter.harmony)
  if (plot_convergence) plot(HarmonyConvergencePlot(harmonyObj))

  ## Return either the R6 Harmony object or the corrected PCA matrix (default)
  if (return_object) {
    return(harmonyObj)
  } else {
    res <- as.matrix(harmonyObj$Z_corr)
    row.names(res) <- row.names(pc_mat)
    colnames(res) <- colnames(pc_mat)
    if (!cells_as_cols) {
      res <- t(res)
    }
    return(res)
  }
}

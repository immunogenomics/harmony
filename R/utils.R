#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom dplyr %>%
#' @examples
#' x <- 5 %>% sum(10)
#' 
#' @usage lhs \%>\% rhs
#' @return return value of rhs function. 
NULL


onehot <- function(x) {
    data.frame(x) %>% 
        tibble::rowid_to_column("row_id") %>% 
        dplyr::mutate(dummy = 1) %>% 
        tidyr::spread(x, .data$dummy, fill = 0) %>% 
        dplyr::select(-.data$row_id) %>% as.matrix
}

scaleData <- function(A, margin = 1, thresh = 10) {
    if (!"dgCMatrix" %in% class(A))
        A <- methods::as(A, "dgCMatrix")
    
    if (margin != 1) A <- t(A)
    
    res <- scaleRows_dgc(A@x, A@p, A@i, ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}

moe_correct_ridge <- function(harmonyObj) {
    harmonyObj$moe_correct_ridge_cpp()
}


moe_ridge_get_betas <- function(harmonyObj) {
    harmonyObj$moe_ridge_get_betas_cpp()
}


cluster <- function(harmonyObj) {
    if (harmonyObj$ran_init == FALSE) {
        stop('before clustering, run init_cluster')
    }
    harmonyObj$cluster_cpp()
}

harmonize <- function(harmonyObj, iter_harmony, verbose=TRUE) {
    if (iter_harmony < 1) {
        return(0)
    }
    
    for (iter in seq_len(iter_harmony)) {
        if (verbose) {
            message(gettextf('Harmony %d/%d', iter, iter_harmony))        
        }
        
        # STEP 1: do clustering
        err_status <- cluster(harmonyObj)
        if (err_status == -1) {
            stop('terminated by user')
        } else if (err_status != 0) {
            stop(gettextf('Harmony exited with non-zero exit status: %d', 
                            err_status))
        }
        
        # STEP 2: regress out covariates
        moe_correct_ridge(harmonyObj)
        
        # STEP 3: check for convergence
        if (harmonyObj$check_convergence(1)) {
            if (verbose) {
                message(gettextf("Harmony converged after %d iterations", 
                        iter))    
            }
            return(0)
        }
    }
}


init_cluster <- function(harmonyObj, cluster_prior=NULL, use_weights=FALSE) {
    if (harmonyObj$ran_setup == FALSE) {
        stop('before initializing cluster, run setup')
    }
    if (!is.null(cluster_prior)) {
        if (ncol(cluster_prior) != harmonyObj$N) {
            stop('cluster_prior must be defined by N cells')
        }
        if (nrow(cluster_prior) > harmonyObj$K) {
            stop('cluster_prior cannot contain more than K clusters')
        }
        C <- nrow(cluster_prior)
        harmonyObj$Y <- matrix(0, harmonyObj$d, harmonyObj$K)
        harmonyObj$Y[, seq_len(C)] <- compute_Y(harmonyObj$Z_cos,
                                                           cluster_prior)
        harmonyObj$R <- matrix(0, harmonyObj$K, harmonyObj$N)
        harmonyObj$R[seq_len(nrow(cluster_prior)), ] <- cluster_prior
        

        ## if needed, initialize K-C clusters        
        if (C < harmonyObj$K) {
#             cluster_res <- flexclust::cclust(
#                 t(harmonyObj$Z_cos), 
#                 harmonyObj$K - C, 
#                 dist = "euclidean", 
#                 method = "hardcl",
#                 weights=harmonyObj$weights
#             )
#             Ynew <- t(cluster_res@centers)
            if (use_weights) {
                Ynew <- soft_kmeans_weighted(
                    harmonyObj$Z_cos, 
                    harmonyObj$K - C,
                    harmonyObj$weights)$Y
            } else {
                Ynew <- t(stats::kmeans(t(harmonyObj$Z_cos), 
                                        centers = harmonyObj$K - C,
                                        iter.max = 25, nstart = 10)$centers)
            }
            
            harmonyObj$Y[, seq(1+C, harmonyObj$K)] <- Ynew
        }
        
        harmonyObj$init_cluster_cpp(C)
    } else {
#         cluster_res <- flexclust::cclust(
#             t(harmonyObj$Z_cos), 
#             harmonyObj$K, 
#             dist = "euclidean", 
#             method = "hardcl",
#             weights=harmonyObj$weights
#         )
#         harmonyObj$Y <- t(cluster_res@centers)
        if (use_weights) {
            harmonyObj$Y <- soft_kmeans_weighted(
                harmonyObj$Z_cos, 
                harmonyObj$K,
                harmonyObj$weights)$Y
        } else {
            harmonyObj$Y <- t(stats::kmeans(t(harmonyObj$Z_cos), 
                                            centers = harmonyObj$K,
                                            iter.max = 25, nstart = 10)$centers)
        }
        harmonyObj$init_cluster_cpp(0)
    }
}


HarmonyConvergencePlot <- function(
        harmonyObj, round_start=1, round_end=Inf, do_wrap=FALSE
    ) {  
    ## ignore initial value
    ## break down kmeans objective into rounds
    obj_fxn <- data.frame(
        kmeans_idx = Reduce(c, lapply(harmonyObj$kmeans_rounds, 
                        function(rounds) {
            seq_len(rounds)
        })),
        harmony_idx = Reduce(c, lapply(
            seq_len(length(harmonyObj$kmeans_rounds)),
            function(i) {rep(i, harmonyObj$kmeans_rounds[i])})
        ),
        val = utils::tail(harmonyObj$objective_kmeans, -1)
    ) %>%
        dplyr::filter(.data$harmony_idx >= round_start) %>% 
        dplyr::filter(.data$harmony_idx <= round_end) %>% 
        tibble::rowid_to_column("idx") 
    
    
    plt <- obj_fxn %>% ggplot2::ggplot(ggplot2::aes(.data$idx, .data$val,
                                                    col = .data$harmony_idx)) + 
        ggplot2::geom_point(shape = 21) + 
        ggplot2::labs(y = "Objective Function", x = "Iteration Number")
    
    if (do_wrap) {
        plt <- plt + ggplot2::facet_grid(.~.data$harmony_idx, scales = 'free',
            space = 'free_x')
    } 
    return(plt)
}

soft_kmeans_weighted <- function(X, k, w, sigma=0.1, iter.max=25, nstart=10, tol=1e-3) {
    seeds <- sample(1e6, nstart)
#     if (!missing(w)) {
        score_curr <- Inf
        for (iter in seq_len(nstart)) {
            set.seed(seeds[iter])
            mediods <- sample(ncol(X), k, prob = w)
            res_new <- soft_kmeans_weighted_cpp(X[, mediods], X, w, iter.max, sigma, tol)
            score_new <- tail(res_new$scores, 1)
            if (score_new < score_curr) {
                res <- res_new
                score_curr <- score_new
            }
        }
#     } else {
#         Y <- stats::kmeans(t(Z), k, iter.max, nstart)$centers %>% t() ## D x K
#         res <- soft_kmeans_cpp(Y, Z, iter.max, sigma)        
#     }
    return(res)
}









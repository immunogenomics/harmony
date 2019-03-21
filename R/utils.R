onehot <- function(x) {
  data.frame(x) %>% 
    tibble::rowid_to_column("id") %>% 
    dplyr::mutate(dummy = 1) %>% 
    tidyr::spread(x, dummy, fill = 0) %>% 
    dplyr::select(-id) %>% as.matrix
}


moe_correct_ridge <- function(harmonyObj) {
  harmonyObj$moe_correct_ridge_cpp()
}


cluster <- function(harmonyObj) {
  if (harmonyObj$ran_init == FALSE) {
    stop('Before clustering, run init_cluster')
  }
  harmonyObj$cluster_cpp()
}

harmonize <- function(harmonyObj, iter_harmony, verbose=TRUE) {
  if (iter_harmony < 1) {
    return(0)
  }
  
  for (iter in 1:iter_harmony) {
    if (verbose) {
        message(sprintf('Harmony %d/%d', iter, iter_harmony))        
    }
    
    # STEP 1: do clustering
    err_status <- cluster(harmonyObj)
    if (err_status == -1) {
      stop('terminated by user')
    } else if (err_status != 0) {
      stop(sprintf('Harmony exited with non-zero exit status: %d', err_status))
    }
    
    # STEP 2: regress out covariates
    moe_correct_ridge(harmonyObj)
    
    # STEP 3: check for convergence
    if (harmonyObj$check_convergence(1)) {
      if (verbose) {
          message(sprintf("Harmony converged after %d iterations", iter))    
      }
      return(0)
    }
  }
}


init_cluster <- function(harmonyObj) {
  if (harmonyObj$ran_setup == FALSE) {
    stop('Before initializing cluster, run setup')
  }
  
  harmonyObj$Y <- t(kmeans(t(harmonyObj$Z_cos), centers = harmonyObj$K, iter.max = 25, nstart = 10)$centers)
  harmonyObj$init_cluster_cpp()
}


HarmonyConvergencePlot <- function(harmonyObj, round_start=1, round_end=Inf, do_wrap=FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("WARNING: Failed to load ggplot2 package. Cannot plot convergence without this package.")
    return(NA)
  }
  
  ## ignore initial value
  ## break down kmeans objective into rounds
  obj_fxn <- data.frame(
    kmeans_idx = Reduce(c, lapply(harmonyObj$kmeans_rounds, function(rounds) {1:rounds})),
    harmony_idx = Reduce(c, lapply(1:length(harmonyObj$kmeans_rounds), function(i) {rep(i, harmonyObj$kmeans_rounds[i])})),
    val = tail(harmonyObj$objective_kmeans, -1)
  ) %>%
    subset(harmony_idx >= round_start & harmony_idx <= round_end) %>% 
    tibble::rowid_to_column("idx") 
  
  
  plt <- obj_fxn %>% ggplot(ggplot2::aes(idx, val, col = harmony_idx)) + 
    geom_point(shape = 21) + 
    labs(y = "Objective Function", x = "Iteration Number")
  
  if (do_wrap) {
    plt <- plt + facet_grid(.~harmony_idx, scales = 'free', space = 'free_x')
  } 
  return(plt)
}














void cosine_normalize(MATTYPE& X, int margin, bool do_safe) {
  // to avoid Inf values, first divide by max 
  if (margin == 1) {
    for (int r = 0; r < X.n_rows; r++) {
      if (do_safe)
        X.row(r) = X.row(r) / X.row(r).max(); 
      X.row(r) = X.row(r) / norm(X.row(r), 2);
    }     
  } else {
    for (int c = 0; c < X.n_cols; c++) {
      if (do_safe)
        X.col(c) = X.col(c) / X.col(c).max(); 
      X.col(c) = X.col(c) / norm(X.col(c), 2);
    }    
  } 
}


MATTYPE safe_entropy(const MATTYPE& X) {
  MATTYPE A = X % log(X);
  A.elem(find_nonfinite(A)).zeros();
  return(A);
}


/*
arma::uvec std_setdiff(const arma::uvec x, const arma::uvec& y) {
  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  return arma::conv_to<arma::uvec>::from(out);
}
*/

// Overload pow to work on a MATTYPErix and vector
MATTYPE pow(MATTYPE A, const VECTYPE & T) {
  for (int c = 0; c < A.n_cols; c++) {
      A.col(c) = pow(A.col(c), as_scalar(T.row(c)));
  }
  return(A);
}

/*
MATTYPE fuzzy_kmeans(const MATTYPE & X, float sigma, const int K, int max_iter = 20, 
                 float converge_thresh = 1e-10) {
  
  // Assume that X is already cosine normalized
  // Should do some random restarts or kmeans++ initialization
  
  int d = X.n_rows;
  int N = X.n_cols;
  MATTYPE Y = zeros<MATTYPE>(d, K);
  MATTYPE R = zeros<MATTYPE>(K, N);  
  uvec rand_idx = conv_to<uvec>::from(randi(K, distr_param(0, K - 1)));
  Y = X.cols(rand_idx);
  cosine_normalize(Y, 0, false); // normalize columns
  vector<float> objective;
  
//  Y.print("Y: ");
//  R.print("R: ");
  
  for (int iter = 0; iter < max_iter; iter++) {    
    R = - (1 / sigma) * 2 * (1 - Y.t() * X);  
    R.each_row() -= max(R, 0);
    R = exp(R);    
    R = normalise(R, 1, 0);
    
    float kmeans_error = as_scalar(accu(R % (2 * (1 - (Y.t() * X)))));
    float _entropy = as_scalar(accu(safe_entropy(R)));  

    objective.push_back(kmeans_error + sigma * _entropy);    
    float obj_old = objective[objective.size() - 2];
    float obj_new = objective[objective.size() - 1];
    
    float obj_change = -(obj_new - obj_old) / abs(obj_old);
    if (obj_change < converge_thresh) break;
    Y = X * R.t();    
    cosine_normalize(Y, 2, true);  
  }
  
  return(R);
}

*/


MATTYPE merge_R(const MATTYPE & R, float thresh = 0.8) {
  MATTYPE cor_res = cor(R.t());
  int K = R.n_rows;
  
  // define equivalence classes  
  uvec equiv_classes = linspace<uvec>(0, K - 1, K);
  int new_idx;
  for (int i = 0; i < K - 1; i++) {
    for (int j = i + 1; j < K; j++) {
      if (cor_res(i, j) > thresh) {
          new_idx = min(equiv_classes(i), equiv_classes(j));
          equiv_classes(i) = new_idx;
          equiv_classes(j) = new_idx;
      }
    }
  }

  // sum over equivalence classes
  uvec uclasses = unique(equiv_classes);
  MATTYPE R_new = zeros<MATTYPE>(uclasses.n_elem, R.n_cols); 
  for (int i = 0; i < R_new.n_rows; i++) {
      uvec idx = find(equiv_classes == uclasses(i)); 
      R_new.row(i) = sum(R.rows(idx), 0);
  }
  
  return R_new;
  
}

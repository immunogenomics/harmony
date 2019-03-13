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


// Overload pow to work on a MATTYPErix and vector
MATTYPE pow(MATTYPE A, const VECTYPE & T) {
  for (int c = 0; c < A.n_cols; c++) {
    A.col(c) = pow(A.col(c), as_scalar(T.row(c)));
  }
  return(A);
}


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

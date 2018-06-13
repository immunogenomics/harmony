void cosine_normalize(mat& X, int margin, bool do_safe) {
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


mat safe_entropy(const mat& X) {
  mat A = X % log(X);
  A.elem(find_nonfinite(A)).zeros();
  return(A);
}


arma::uvec std_setdiff(const arma::uvec x, const arma::uvec& y) {
  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  return arma::conv_to<arma::uvec>::from(out);
}

// Overload pow to work on a matrix and vector
mat pow(mat A, const vec & T) {
  for (int c = 0; c < A.n_cols; c++) {
      A.col(c) = pow(A.col(c), as_scalar(T.row(c)));
  }
  return(A);
}


#include "utils.h"
#include "types.h"

//[[Rcpp::export]]
arma::mat kmeans_centers(const arma::mat& X, const int K){
  
  // Environment 
  Rcpp::Environment stats_env("package:stats");
  // Cast function as callable from C++
  Rcpp::Function kmeans = stats_env["kmeans"];
  // Call the function and receive its list output
  Rcpp::List res = kmeans(Rcpp::_["x"] = X.t(),
                          Rcpp::_["centers"] = K,
                          Rcpp::_["iter.max"] = 25,
                          Rcpp::_["nstart"] = 10
                          );
  return res["centers"];
}


MATTYPE safe_entropy(const MATTYPE& X) {
  MATTYPE A = X % log(X);
  A.elem(find_nonfinite(A)).zeros();
  return(A);
}

// Overload pow to work on a MATTYPErix and vector
MATTYPE harmony_pow(MATTYPE A, const VECTYPE& T) {

  for (unsigned c = 0; c < A.n_cols; c++) {
    A.unsafe_col(c) = pow(A.unsafe_col(c), as_scalar(T.row(c)));
  }
  return(A);
}

VECTYPE calculate_norm(const MATTYPE& M){
  VECTYPE x(M.n_cols);
  for(unsigned i = 0; i < M.n_cols; i++){
    x(i) = norm(M.col(i));
  }
  return x;
}


//https://stackoverflow.com/questions/8377412/ceil-function-how-can-we-implement-it-ourselves
int my_ceil(float num) {
    int inum = (int)num;
    if (num == (float)inum) {
        return inum;
    }
    return inum + 1;
}

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


// [[Rcpp::export]]
MATTYPE scaleRows_dgc(const VECTYPE& x, const VECTYPE& p, const VECTYPE& i, int ncol, int nrow, float thresh) {
  
    // (0) fill in non-zero elements
    MATTYPE res = arma::zeros<MATTYPE>(nrow, ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            res(i[j], c) = x(j);
        }
    }

    // (1) compute means
    VECTYPE mean_vec = arma::zeros<VECTYPE>(nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            mean_vec(i[j]) += x[j];
        }
    }
    mean_vec /= ncol;

    // (2) compute SDs
    VECTYPE sd_vec = arma::zeros<VECTYPE>(nrow);
    arma::uvec nz = arma::zeros<arma::uvec>(nrow);
    nz.fill(ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            sd_vec(i[j]) += (x[j] - mean_vec(i[j])) * (x[j] - mean_vec(i[j])); // (x - mu)^2
            nz(i[j])--;
        }
    }

    // count for the zeros
    for (int r = 0; r < nrow; r++) {
        sd_vec(r) += nz(r) * mean_vec(r) * mean_vec(r);
    }

    sd_vec = arma::sqrt(sd_vec / (ncol - 1));

    // (3) scale values
    res.each_col() -= mean_vec;
    res.each_col() /= sd_vec;
    res.elem(find(res > thresh)).fill(thresh);
    res.elem(find(res < -thresh)).fill(-thresh);
    return res;
}


// [[Rcpp::export]]
double find_one_lambda_cpp(const arma::vec& cluster_size, const arma::vec& range) {
    double batch_max = cluster_size.max();
    double batch_min = cluster_size.min();
    double lambda = pow(batch_max, 0.5) * pow(batch_min, 0.5);
    lambda = std::min(range(1), lambda);
    lambda = std::max(range(0), lambda);
    return lambda;
}

// [[Rcpp::export]]
arma::vec find_lambda_cpp(const arma::vec& cluster_size, const arma::vec& range,
                          const std::vector<int>& B_vec) {
  arma::vec lambda_dym_vec(cluster_size.n_rows + 1, arma::fill::zeros);
  int current_idx = 1; // base indx for lambda_dy_vec; 1 to omit the intercept
  // find lambda_dym for every batch variable
  for(unsigned int b = 0; b < B_vec.size(); b++){
    arma::vec sub_cluster_size = cluster_size.subvec(current_idx,
                                           current_idx + B_vec[b]-1);
    double lambda_dym = find_one_lambda_cpp(sub_cluster_size, range);
    lambda_dym_vec.subvec(current_idx, current_idx + B_vec[b]-1) = arma::vec(
        B_vec[b], arma::fill::value(lambda_dym));
    current_idx += B_vec[b];
  }
  return lambda_dym_vec;
}

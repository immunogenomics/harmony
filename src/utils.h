#pragma once
#include "types.h"
#include <RcppArmadillo.h>

arma::mat kmeans_centers(const arma::mat& X, const int K);

MATTYPE safe_entropy(const MATTYPE& X);

MATTYPE harmony_pow(MATTYPE A, const VECTYPE& T);

VECTYPE calculate_norm(const MATTYPE& M);


int my_ceil(float num);



// [[Rcpp::export]]
double find_one_lambda_cpp(arma::vec cluster_size, arma::vec range){
    double batch_max = cluster_size.max();
    double batch_min = cluster_size.min();
    double lambda = pow(batch_max, 0.5) * pow(batch_min, 0.5);
    lambda = std::min(range(1), lambda);
    lambda = std::max(range(0), lambda);
    return lambda;
}

// [[Rcpp::export]]
arma::vec find_lambda_cpp(arma::vec cluster_size, arma::vec range,
                       std::vector<int> B_vec){
  arma::vec lambda_dym_vec(cluster_size.n_rows + 1, arma::fill::zeros);
  int current_idx = 0;
  arma::vec sub_cluster_size;
  for(int b = 0; b < B_vec.size(); b++){
    sub_cluster_size = cluster_size.subvec(current_idx,
                                           current_idx + B_vec[b]-1);
    double lambda_dym = find_one_lambda_cpp(sub_cluster_size, range);
    lambda_dym_vec.subvec(current_idx + 1, current_idx + B_vec[b]) = arma::vec(
        B_vec[b], arma::fill::value(lambda_dym)
    );
    current_idx += B_vec[b];
  }
  return lambda_dym_vec;
}
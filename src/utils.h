#pragma once
#include "types.h"
#include <RcppArmadillo.h>

arma::mat kmeans_centers(const arma::mat& X, const int K);

MATTYPE safe_entropy(const MATTYPE& X);

MATTYPE harmony_pow(MATTYPE A, const VECTYPE& T);

VECTYPE calculate_norm(const MATTYPE& M);


int my_ceil(float num);


double find_one_lambda_cpp(const arma::vec& cluster_size, const arma::vec& range);

arma::vec find_lambda_cpp(const arma::vec& cluster_size, const arma::vec& range,
                          const std::vector<int>& B_vec);

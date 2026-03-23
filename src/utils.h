#pragma once
#include "types.h"
#include <RcppArmadillo.h>
#define DEBUG false


MATTYPE kmeans_centers(const MATTYPE& X, const unsigned int K, bool verbose);

MATTYPE safe_entropy(const MATTYPE& X);

MATTYPE harmony_pow(MATTYPE A, const VECTYPE& T);

VECTYPE calculate_norm(const MATTYPE& M);


int my_ceil(float num);
float my_accu(const MATTYPE& X);

VECTYPE find_lambda_cpp(const float alpha, const VECTYPE& cluster_E);

std::vector< std::pair<unsigned,unsigned> > find_contigs(std::vector<unsigned>& keep_vectors);

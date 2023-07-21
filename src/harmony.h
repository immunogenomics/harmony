#include "types.h"
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <progress.hpp>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;

class harmony;
RCPP_EXPOSED_CLASS(harmony)
  
  
#include "harmony_types.h"
  
class harmony { 
public:
  
  harmony();
  
    void setup(const MATTYPE& __Z, const arma::sp_mat& __Phi,
	       const VECTYPE __sigma, const VECTYPE __theta, const int __max_iter_kmeans, 
	       const float __epsilon_kmeans, const float __epsilon_harmony,
	       const int __K, const float __block_size,
	       const MATTYPE __lambda, const bool __verbose, arma::vec __lambda_range,
         std::vector<int> __B_vec);
  
  /* METHODS */
  void moe_correct_ridge_cpp();
  // TODO may merge this function with moe_correct_ridge_cpp
  void mid_cap_moe_correct_ridge_cpp();
  CUBETYPE moe_ridge_get_betas_cpp();
  int cluster_cpp();

  void init_cluster_cpp(unsigned);
  void allocate_buffers();
  void compute_objective(); 
  int update_R();
  bool check_convergence(int type);
  void setY(const MATTYPE& Z);

  /* FIELDS */
  MATTYPE R, Z_orig, Z_corr, Z_cos, Y, all_lambda_dym_mat;
  arma::sp_mat Phi, Phi_moe, Phi_moe_t, Phi_t, lambda, Rk;
  VECTYPE Pr_b, theta, N_b, sigma, lambda_range;
  
  vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, objective_kmeans_cross, objective_harmony;
  vector<int> kmeans_rounds, B_vec; // OLD: Kb
  
  float block_size, epsilon_kmeans, epsilon_harmony;
  unsigned int N, K, B, d, max_iter_kmeans, window_size;

  // buffers
  MATTYPE _scale_dist, dist_mat, O, E, dir_prior, lambda_dym_mat; // N_k, N_kb, N_b, numerator, denominator, C;
  uvec update_order, cells_update;
  MATTYPE W;
    
  // flags
  bool ran_setup, ran_init, verbose; // do_merge_R;
  
};

  

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
  
    void setup(MATTYPE& __Z, arma::sp_mat& __Phi,
	       VECTYPE __sigma, VECTYPE __theta, int __max_iter_kmeans, 
	       float __epsilon_kmeans, float __epsilon_harmony,
	       int __K, float __block_size,
	       MATTYPE __lambda, bool __verbose);
  
  /* METHODS */
  void moe_correct_ridge_cpp();
  int cluster_cpp();

  void init_cluster_cpp(unsigned);
  void allocate_buffers();
  void compute_objective(); 
  int update_R();
  bool check_convergence(int type);
  void setY(const MATTYPE& Z);

  /* FIELDS */
  MATTYPE R, Z_orig, Z_corr, Z_cos, Y;
  arma::sp_mat Phi, Phi_moe, Phi_moe_t, Phi_t, lambda, Rk;
  VECTYPE Pr_b, theta, N_b, sigma;
  
  vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, objective_kmeans_cross, objective_harmony;
  vector<int> kmeans_rounds; // OLD: Kb
  
  float block_size, epsilon_kmeans, epsilon_harmony;
  unsigned int N, K, B, d, max_iter_kmeans, window_size;

  // buffers
  MATTYPE _scale_dist, dist_mat, O, E, dir_prior; // N_k, N_kb, N_b, numerator, denominator, C;
  uvec update_order, cells_update;
  MATTYPE W;
    
  // flags
  bool ran_setup, ran_init, verbose, Yset; // do_merge_R;
  
};

  

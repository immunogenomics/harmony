#include <RcppArmadillo.h>
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
  /* CONSTRUCTORS etc */
  harmony(int __K);    
  void setup(MATTYPE& __Z, MATTYPE& __Phi, VECTYPE __Pr_b,
             VECTYPE __sigma, VECTYPE __theta, int __max_iter_kmeans, 
             float __epsilon_kmeans, float __epsilon_harmony, 
             int __K, float tau, float __block_size, 
             MATTYPE __lambda, bool __verbose);
  
  /* METHODS */
  void moe_correct_ridge_cpp();
  void init_cluster_cpp();
  int cluster_cpp();
  
  void allocate_buffers();
  void compute_objective(); 
  int update_R();
  bool check_convergence(int type);

  /* FIELDS */
  MATTYPE R, Z_orig, Z_corr, Z_cos, Y, Y_unnormed, Phi, Phi_moe; 
  VECTYPE Pr_b, theta, N_b, sigma, sigma_prior;
  MATTYPE lambda; // diagonal MATTYPErix of ridge regression penalties
  vector<float> objective_harmony;
  vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, objective_kmeans_cross;
  vector<int> kmeans_rounds; // OLD: Kb
  
  //    vector<uvec> phi_map;
  float block_size, epsilon_kmeans, epsilon_harmony, merge_thresh_global;
  int N, K, B, d, max_iter_kmeans, window_size; 

  // buffers
  MATTYPE _scale_dist, dist_mat, O, E, dir_prior, Phi_Rk; // N_k, N_kb, N_b, numerator, denominator, C;
  uvec update_order, cells_update;
  MATTYPE W;
  
  // flags
  bool ran_setup, ran_init, verbose; // do_merge_R;
  
};

  
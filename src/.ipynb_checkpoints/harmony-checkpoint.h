#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <progress.hpp>

//#include <progress_bar.hpp>

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
                float __epsilon_kmeans, float __epsilon_harmony, bool __correct_with_Zorig,
                float __alpha, int __K, float tau, float __block_size, 
                //ROWVECTYPE& w_new, 
                bool __correct_with_cosine,
                int __window_size, MATTYPE __lambda);
    
    /* METHODS */
    void harmonize(int harmony_iter);
//    void gmm_correct_armadillo();
//    void moe_correct_contrast();
//    void moe_correct_onehot();
    void init_clusters_random_balanced();
    void moe_correct_ridge();
    void init_cluster();
    int cluster();

    void allocate_buffers();
    void compute_objective(); 
    int compute_R(bool random_order);
    bool check_convergence(int type);
    
//    void set_R_merge_flag(float merge_thresh_new);
//    void update_R_merge();  

    /* FIELDS */
    MATTYPE R, Z_orig, Z_corr, Z_cos, Y, Phi, Phi_moe; 
    VECTYPE Pr_b, theta, N_b, sigma, sigma_prior;
//    ROWVECTYPE w;
    MATTYPE lambda; // diagonal MATTYPErix of ridge regression penalties
    vector<float> objective_harmony;
    vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, objective_kmeans_cross;
    vector<int> kmeans_rounds; // OLD: Kb

//    vector<uvec> phi_map;
    float block_size, alpha, epsilon_kmeans, epsilon_harmony, merge_thresh_global;
    int N, K, B, d, max_iter_kmeans, window_size; 
    bool correct_with_Zorig, correct_with_cosine;
  
    // buffers
    MATTYPE _scale_dist, dist_mat, O, E, dir_prior, Phi_Rk; // N_k, N_kb, N_b, numerator, denominator, C;
    uvec update_order, cells_update;
//    CUBETYPE W;
//    MATTYPE mu_k, mu_k_r, mu_bk_r,
//    CUBETYPE mu_bk;
    MATTYPE W;

    // flags
    bool ran_setup, ran_init; // do_merge_R;
  
    // FOR DEBUGGING ONLY - SHOULD ERASE THESE
//    vector<MATTYPE> R_list;
};


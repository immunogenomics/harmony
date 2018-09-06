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


class harmony { 
  public: 
    /* CONSTRUCTORS etc */
    harmony(int __K);    
    void setup(arma::fmat& Z, arma::fmat& Phi, 
               double __sigma, double __theta, int __max_iter_kmeans, 
                float __epsilon_kmeans, float __epsilon_harmony, bool __correct_with_Zorig,
                double __alpha, int __K, float tau, float __block_size, 
               rowvec& w_new, bool __correct_with_cosine, vector<bool> batch_mask_new, int __window_size);
    void dummy(int __K) {
      K = __K;
    }
  
    /* METHODS */
    void harmonize(int harmony_iter);
    void gmm_correct_armadillo();
    void init_cluster();
    int cluster();

    void allocate_buffers();
    void set_thetas(float theta_max, float tau);  
    fmat compute_C(uvec& cells_in, uvec& cells_out);
    void compute_objective(); 
    int compute_R();
    bool check_convergence(int type);
  /*
    void init_batch_clusters(uvec & batch, float merge_thresh,
                              float sigma_local, int K_local);
    void compute_phi_hat(const uvec & batches, float merge_thresh,
                              float sigma_local, int K_local);
                              */
    void set_R_merge_flag(float merge_thresh_new);
    void update_R_merge();  

    /* FIELDS */
    fmat R, Z_orig, Z_corr, Z_cos, Y, Phi, phi_hat; 
    vec N_b, Pr_b, Pr_Kb, theta, theta2, N_Kb;
    rowvec w;
    vector<float> objective_harmony;
    vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, objective_kmeans_cross;
    vector<int> Kb, kmeans_rounds;
    vector<bool> batch_mask;
    vector<uvec> phi_map;
    float sigma, block_size, alpha, epsilon_kmeans, epsilon_harmony, theta_max, merge_thresh_global;
    int N, K, B, d, max_iter_kmeans, window_size; 
    bool correct_with_Zorig, correct_with_cosine;
  
    // buffers
    fmat _scale_dist, mu_k, mu_k_r, mu_bk_r, O, E, O2, E2, dir_prior, dir_prior2; // N_k, N_kb, N_b, numerator, denominator, C;
    uvec update_order;
    cube mu_bk;
  
    // flags
    bool ran_setup, ran_init, do_conservation, do_merge_R;
  
    // FOR DEBUGGING ONLY - SHOULD ERASE THESE
//    vector<fmat> R_list;
    uvec cells_update;
  
};


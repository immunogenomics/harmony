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
    void setup(fmat& Z, fmat& Phi, 
               float __sigma, float __theta, int __max_iter_kmeans, 
                float __epsilon_kmeans, float __epsilon_harmony, bool __correct_with_Zorig,
                float __alpha, int __K, float tau, float __block_size, 
               frowvec& w_new, bool __correct_with_cosine, vector<bool> batch_mask_new, int __window_size);
    void dummy(int __K) {
      K = __K;
    }
  
    /* METHODS */
    void harmonize(int harmony_iter);
    void gmm_correct_armadillo();
    void init_cluster();
    int cluster();

    void allocate_buffers();
    fvec set_thetas(float theta_max, float tau, fvec& Nb_use);  
    fmat compute_C(uvec& cells_in, uvec& cells_out);
    void compute_objective(); 
    int compute_R();
    bool check_convergence(int type);
  
  
//    void init_batch_clusters(uvec & batch, float merge_thresh,
//                              float sigma_local, int K_local);
//    void compute_phi_hat(const uvec & batches, float merge_thresh,
//                              float sigma_local, int K_local);

    void setup_batch2(fmat& Phi2_new, float theta2_new, float tau);
  
    void set_R_merge_flag(float merge_thresh_new);
    void update_R_merge();  

    /* FIELDS */
    fmat R, Z_orig, Z_corr, Z_cos, Y, Phi, Phi2; 
    fvec Pr_b, Pr_b2, theta, theta2, N_b, N_b2;
    frowvec w;
    vector<float> objective_harmony;
    vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, objective_kmeans_cross;
    vector<int> kmeans_rounds; // OLD: Kb
    vector<bool> batch_mask;
    vector<uvec> phi_map;
    float sigma, block_size, alpha, epsilon_kmeans, epsilon_harmony, merge_thresh_global;
    int N, K, B, B2, d, max_iter_kmeans, window_size; 
    bool correct_with_Zorig, correct_with_cosine;
  
    // buffers
    fmat _scale_dist, mu_k, mu_k_r, mu_bk_r, O, E, O2, E2, dir_prior, dir_prior2; // N_k, N_kb, N_b, numerator, denominator, C;
    uvec update_order;
    fcube mu_bk;
  
    // flags
    bool ran_setup, ran_init, do_merge_R, do_theta2;
  
    // FOR DEBUGGING ONLY - SHOULD ERASE THESE
    vector<fmat> R_list;
    uvec cells_update;
  
};


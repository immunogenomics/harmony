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
    harmony(int __K);  
  
    void setup(arma::mat& Z, arma::mat& Phi, 
               double __sigma, double __theta, int __max_iter_kmeans, 
                float __converge_thresh, bool __correct_with_Zorig,
                double __alpha, int __K, float tau, float __block_size);
    void dummy(int __K) {
      K = __K;
    }
    void harmonize(int harmony_iter);
    void gmm_correct_armadillo();
    void init_cluster();
    int cluster();
  
    mat compute_C(uvec& cells_in, uvec& cells_out);
    void compute_objective(); 
    void compute_R();
    bool check_convergence(int type);

    mat R, Z_orig, Z_corr, Z_cos, Y, Phi; 
    vec Pr_b, theta;
    vector<float> objective_harmony;
    vector<float> objective_kmeans;
    float sigma, block_size, alpha, converge_thresh;
    int N, K, B, d, max_iter_kmeans; 
    bool correct_with_Zorig;
  
    // buffers
    mat _scale_dist, mu_k, mu_k_r, mu_bk_r, O, E; // N_k, N_kb, N_b, numerator, denominator, C;
    cube mu_bk;
  
    // flags
    bool ran_setup, ran_init;
  
    mat Rtest1(mat R);
    mat Rtest2(mat R);
  
  
};


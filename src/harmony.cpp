#include <algorithm>
#include <chrono>

#include "harmony.h"
#include "types.h"
#include "utils.h"




harmony::harmony() :
    window_size(3),
    ran_setup(false),
    ran_init(false),
    verbose(false)
{}



void harmony::setup(const MATTYPE& __Z, const arma::sp_mat& __Phi,
                    const VECTYPE __sigma, const VECTYPE __theta, const int __max_iter_kmeans,
                    const float __epsilon_kmeans, const float __epsilon_harmony,
                    const int __K, const float __block_size,
                    const MATTYPE __lambda, const bool __verbose) {
    
  // Algorithm constants
  N = __Z.n_cols;
  B = __Phi.n_rows;
  d = __Z.n_rows;

  Z_orig = __Z;
  Z_cos = arma::normalise(__Z, 2, 0);
  Z_corr = zeros(size(Z_orig));

  
  Phi = __Phi;
  Phi_t = Phi.t();
  
  Pr_b = sum(Phi, 1) / N;

  
  epsilon_kmeans = __epsilon_kmeans;
  epsilon_harmony = __epsilon_harmony;

  // Hyperparameters
  K = __K;
  lambda = __lambda;
  sigma = __sigma;

  block_size = __block_size;
  theta = __theta;
  max_iter_kmeans = __max_iter_kmeans;

  verbose = __verbose;
  
  allocate_buffers();
  ran_setup = true;
  
}


void harmony::allocate_buffers() {
  
  _scale_dist = zeros<MATTYPE>(K, N);
  dist_mat = zeros<MATTYPE>(K, N);
  O = E = zeros<MATTYPE>(K, B);
  
  // Hack: create matrix of ones by creating zeros and then add one!
  arma::sp_mat intcpt = zeros<arma::sp_mat>(1, N);
  intcpt = intcpt+1;
  
  Phi_moe = join_cols(intcpt, Phi);
  Phi_moe_t = Phi_moe.t();
  W = zeros<MATTYPE>(B + 1, d);
  
}


void harmony::init_cluster_cpp(unsigned C) {
  
  Rcerr << "Hard k-means centroids initialization"  <<std::endl;
  Y = kmeans_centers(Z_cos, K).t();
  
  // Cosine normalization of data centrods
  Y = arma::normalise(Y, 2, 0);

  // (2) ASSIGN CLUSTER PROBABILITIES
  // using a nice property of cosine distance,
  // compute squared distance directly with cross product
  dist_mat = 2 * (1 - Y.t() * Z_cos);
  
  R = -dist_mat;
  R.each_col() /= sigma;
  R = exp(R);
  R.each_row() /= sum(R, 0);
  
  
  // (3) BATCH DIVERSITY STATISTICS
  E = sum(R, 1) * Pr_b.t();
  O = R * Phi_t;
  
  compute_objective();
  objective_harmony.push_back(objective_kmeans.back());
  
  dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed

  ran_init = true;
  
}

void harmony::compute_objective() {
  float kmeans_error = as_scalar(accu(R % dist_mat));
  float _entropy = as_scalar(accu(safe_entropy(R).each_col() % sigma)); // NEW: vector sigma
  float _cross_entropy = as_scalar(
      accu((R.each_col() % sigma) % ((arma::repmat(theta.t(), K, 1) % log((O + 1) / (E + 1))) * Phi)));

  // Push back the data
  objective_kmeans.push_back(kmeans_error + _entropy + _cross_entropy);
  objective_kmeans_dist.push_back(kmeans_error);
  objective_kmeans_entropy.push_back(_entropy);
  objective_kmeans_cross.push_back(_cross_entropy);
}


bool harmony::check_convergence(int type) {
  float obj_new, obj_old;
  switch (type) {
    case 0: 
      // Clustering 
      // compute new window mean
      obj_old = 0;
      obj_new = 0;
      for (unsigned i = 0; i < window_size; i++) {
        obj_old += objective_kmeans[objective_kmeans.size() - 2 - i];
        obj_new += objective_kmeans[objective_kmeans.size() - 1 - i];
      }
      if ((obj_old - obj_new) / abs(obj_old) < epsilon_kmeans) {
        return(true); 
      } else {
        return(false);
      }
    case 1:
      // Harmony
      obj_old = objective_harmony[objective_harmony.size() - 2];
      obj_new = objective_harmony[objective_harmony.size() - 1];
      if ((obj_old - obj_new) / abs(obj_old) < epsilon_harmony) {
        return(true);              
      } else {
        return(false);              
      }
  }
  
  // gives warning if we don't give default return value
  return(true);
}


int harmony::cluster_cpp() {
  int err_status = 0;
  Progress p(max_iter_kmeans, verbose);
  unsigned iter;
  
  // Z_cos has changed
  // R has assumed to not change
  // so update Y to match new integrated data  
  for (iter = 0; iter < max_iter_kmeans; iter++) {
      
      p.increment();
      if (Progress::check_abort())
	  return(-1);
    
      // STEP 1: Update Y (cluster centroids)
      Y = arma::normalise(Z_cos * R.t(), 2, 0);

      dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed

        
      // STEP 3: Update R    
      err_status = update_R();
    
      if (err_status != 0) {
	  // Rcout << "Compute R failed. Exiting from clustering." << endl;
	  return err_status;
      }
    
      // STEP 4: Check for convergence
      compute_objective();
    
      if (iter > window_size) {
	  bool convergence_status = check_convergence(0);
	  if (convergence_status) {
	      iter++;
	      break;
	  }
      }
  }
  
  kmeans_rounds.push_back(iter);
  objective_harmony.push_back(objective_kmeans.back());
  return 0;
}






int harmony::update_R() {

  // Generate the 0,N-1 indices
  uvec indices = linspace<uvec>(0, N - 1, N);
  update_order = shuffle(indices);
  
  // Inverse index
  uvec reverse_index(N, arma::fill::zeros);
  reverse_index.rows(update_order) = indices;
  
  _scale_dist = -dist_mat; // K x N
  _scale_dist.each_col() /= sigma; // NEW: vector sigma
  _scale_dist = exp(_scale_dist);
  _scale_dist = arma::normalise(_scale_dist, 1, 0);

  // GENERAL CASE: online updates, in blocks of size (N * block_size)
  unsigned n_blocks = (int)(my_ceil(1.0 / block_size));
  unsigned cells_per_block = my_ceil(N / n_blocks);
  
  // Allocate new matrices
  MATTYPE R_randomized = R.cols(update_order);
  arma::sp_mat Phi_randomized(Phi.cols(update_order));
  arma::sp_mat Phi_t_randomized(Phi_randomized.t());
  MATTYPE _scale_dist_randomized = _scale_dist.cols(update_order);

  for (unsigned i = 0; i < n_blocks; i++) {
      unsigned idx_max = min((i+1) * cells_per_block, N-1);
      auto Rcells = R_randomized.submat(0, i*cells_per_block, R_randomized.n_rows - 1, idx_max);
      auto Phicells = Phi_randomized.submat(0, i*cells_per_block, Phi_randomized.n_rows - 1, idx_max);
      auto Phi_tcells = Phi_t_randomized.submat(i*cells_per_block, 0, idx_max, Phi_t_randomized.n_cols - 1);
      auto _scale_distcells = _scale_dist_randomized.submat(0, i*cells_per_block, _scale_dist_randomized.n_rows - 1, idx_max);

      // Step 1: remove cells
      E -= sum(Rcells, 1) * Pr_b.t();
      O -= Rcells * Phi_tcells;

      // Step 2: recompute R for removed cells
      Rcells = _scale_distcells;
      Rcells = Rcells % (harmony_pow((E + 1) / (O + 1), theta) * Phicells);
      Rcells = normalise(Rcells, 1, 0); // L1 norm columns

      // Step 3: put cells back 
      E += sum(Rcells, 1) * Pr_b.t();
      O += Rcells * Phi_tcells;
  }
  this->R = R_randomized.cols(reverse_index);
  return 0;
}


void harmony::moe_correct_ridge_cpp() {

  Progress p(K, false);
  arma::sp_mat _Rk(N, N);
  
  Z_corr = Z_orig;
      
  for (unsigned k = 0; k < K; k++) {
      
      if (Progress::check_abort())
	  return;
      
      _Rk.diag() = R.row(k);
      arma::sp_mat Phi_Rk = Phi_moe * _Rk;
      W = arma::inv(arma::mat(Phi_Rk * Phi_moe_t + lambda)) * Phi_Rk * Z_orig.t();
      W.row(0).zeros(); // do not remove the intercept 
      Z_corr -= W.t() * Phi_Rk;
      
  }
  Z_cos = arma::normalise(Z_corr, 2, 0);
}

CUBETYPE harmony::moe_ridge_get_betas_cpp() {
  CUBETYPE W_cube(W.n_rows, W.n_cols, K); // rows, cols, slices

  arma::sp_mat _Rk(N, N);
  for (unsigned k = 0; k < K; k++) {
      _Rk.diag() = R.row(k);
      arma::sp_mat Phi_Rk = Phi_moe * _Rk;
      W_cube.slice(k) = arma::inv(arma::mat(Phi_Rk * Phi_moe_t + lambda)) * Phi_Rk * Z_orig.t();
  }

  return W_cube;
}

RCPP_MODULE(harmony_module) {
  class_<harmony>("harmony")
      .constructor()
      .field("Z_corr", &harmony::Z_corr)
      .field("Z_cos", &harmony::Z_cos)
      .field("Z_orig", &harmony::Z_orig)
      .field("Phi", &harmony::Phi)
      .field("Phi_moe", &harmony::Phi_moe)
      .field("N", &harmony::N)
      .field("B", &harmony::B)
      .field("K", &harmony::K)
      .field("d", &harmony::d)
      .field("O", &harmony::O)
      .field("E", &harmony::E)
      .field("Y", &harmony::Y)
      .field("Pr_b", &harmony::Pr_b)
      .field("W", &harmony::W)
      .field("R", &harmony::R)
      .field("kmeans_rounds", &harmony::kmeans_rounds)
      .field("objective_kmeans", &harmony::objective_kmeans)
      .field("objective_kmeans_dist", &harmony::objective_kmeans_dist)
      .field("objective_kmeans_entropy", &harmony::objective_kmeans_entropy)
      .field("objective_kmeans_cross", &harmony::objective_kmeans_cross)    
      .field("objective_harmony", &harmony::objective_harmony)
      .method("check_convergence", &harmony::check_convergence)
      .method("setup", &harmony::setup)
      .method("compute_objective", &harmony::compute_objective)
      .method("init_cluster_cpp", &harmony::init_cluster_cpp)
      .method("cluster_cpp", &harmony::cluster_cpp)	  
      .method("moe_correct_ridge_cpp", &harmony::moe_correct_ridge_cpp)
      .method("moe_ridge_get_betas_cpp", &harmony::moe_ridge_get_betas_cpp)
      ;
}

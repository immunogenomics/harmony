#include "harmony.h"
#include "utils.h"


// NOTE: This is a dummy constructor, needed by Rcpp
harmony::harmony(int __K): K(__K) {}



void harmony::setup(MATTYPE& __Z, MATTYPE& __Phi, MATTYPE& __Phi_moe, VECTYPE __Pr_b,
                    VECTYPE __sigma, VECTYPE __theta, int __max_iter_kmeans, 
                    float __epsilon_kmeans, float __epsilon_harmony, 
                    int __K, float tau, float __block_size, 
                    MATTYPE __lambda, bool __verbose) {
  
  Z_corr = MATTYPE(__Z);
  Z_orig = MATTYPE(__Z);
  Z_cos = MATTYPE(Z_orig);
  cosine_normalize(Z_cos, 0, true); // normalize columns
  
  Phi = __Phi;
  Phi_moe = __Phi_moe;
  N = Z_corr.n_cols;
  Pr_b = __Pr_b;
  B = Phi.n_rows;
  d = Z_corr.n_rows; 
  window_size = 3;
  epsilon_kmeans = __epsilon_kmeans;
  epsilon_harmony = __epsilon_harmony;
  
  
  lambda = __lambda;
  sigma = __sigma;
  sigma_prior = __sigma;
  block_size = __block_size;
  K = __K;
  max_iter_kmeans = __max_iter_kmeans;
  verbose = __verbose;
  
  theta = __theta;
  allocate_buffers();
  ran_setup = true;
}


void harmony::allocate_buffers() {
  _scale_dist = zeros<MATTYPE>(K, N);    
  dist_mat = zeros<MATTYPE>(K, N);    
  O = zeros<MATTYPE>(K, B);
  E = zeros<MATTYPE>(K, B);  
  W = zeros<MATTYPE>(B + 1, d); 
  Phi_Rk = zeros<MATTYPE>(B + 1, N);
}




void harmony::init_cluster_cpp() {
  // kmeans is called outside, in the R function
  cosine_normalize(Y, 0, false); // normalize columns
  
  // (2) ASSIGN CLUSTER PROBABILITIES
  // using a nice property of cosine distance,
  // compute squared distance directly with cross product
  dist_mat = 2 * (1 - Y.t() * Z_cos); // initial estimate based on Y_0 and Z_0
  R = - dist_mat;
  R.each_col() /= sigma;
  R.each_row() -= max(R, 0);  
  R = exp(R);
  R.each_row() /= sum(R, 0);
  
  // (3) BATCH DIVERSITY STATISTICS
  E = sum(R, 1) * Pr_b.t();
  O = R * Phi.t();
  
  compute_objective();
  objective_harmony.push_back(objective_kmeans.back());
  ran_init = true;
  
}

void harmony::compute_objective() {
  float kmeans_error = as_scalar(accu(R % dist_mat)); 
  float _entropy = as_scalar(accu(safe_entropy(R).each_col() % sigma)); // NEW: vector sigma
  float _cross_entropy;
  _cross_entropy = as_scalar(accu((R.each_col() % sigma) % ((arma::repmat(theta.t(), K, 1) % log((O + 1) / (E + 1))) * Phi)));
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
    for (int i = 0; i < window_size; i++) {
      obj_old += objective_kmeans[objective_kmeans.size() - 2 - i];
      obj_new += objective_kmeans[objective_kmeans.size() - 1 - i];
    }
    if (-(obj_new - obj_old) / obj_old < epsilon_kmeans) {
      //        Rcout << "kmeans old: " << obj_old << ", new: " << obj_new << ", diff: " << -(obj_new - obj_old) / obj_old << endl;
      return(true); 
    } else {
      return(false);
    }
  case 1:
    // Harmony
    obj_old = objective_harmony[objective_harmony.size() - 2];
    obj_new = objective_harmony[objective_harmony.size() - 1];
    if (-(obj_new - obj_old) / obj_old < epsilon_harmony) {
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
  int iter; 
  Progress p(max_iter_kmeans, verbose);

  // Z_cos has changed
  // R has assumed to not change
  // so update Y to match new integrated data
  dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed
  for (iter = 0; iter < max_iter_kmeans; iter++) {
    p.increment();
    if (Progress::check_abort())
      return(-1);
    
    // STEP 1: Update Y
    Y = normalise(Z_cos * R.t(), 2, 0);
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
        //        Rcout << "... Breaking Clustering ..., status = " << convergence_status << endl;
        iter++;
        // Rcout << "Clustered for " << iter << " iterations" << endl;
        break;        
      }
    }
  }
  kmeans_rounds.push_back(iter);
  objective_harmony.push_back(objective_kmeans.back());
  return 0;
}




int harmony::update_R() {  
  update_order = shuffle(linspace<uvec>(0, N - 1, N));
  _scale_dist = -dist_mat;
  _scale_dist.each_col() /= sigma; // NEW: vector sigma          
  _scale_dist.each_row() -= max(_scale_dist, 0);
  _scale_dist = exp(_scale_dist);

  
  // GENERAL CASE: online updates, in blocks of size (N * block_size)
  for (int i = 0; i < ceil(1. / block_size); i++) {    
    // gather cell updates indices
    int idx_min = i * N * block_size;
    int idx_max = min((int)((i + 1) * N * block_size), N - 1);
    if (idx_min > idx_max) break; // TODO: fix the loop logic so that this never happens
    uvec idx_list = linspace<uvec>(idx_min, idx_max, idx_max - idx_min + 1);
    cells_update = update_order.rows(idx_list); 
    
    // Step 1: remove cells
    E -= sum(R.cols(cells_update), 1) * Pr_b.t();
    O -= R.cols(cells_update) * Phi.cols(cells_update).t();

    // Step 2: recompute R for removed cells
    R.cols(cells_update) = _scale_dist.cols(cells_update);    
    R.cols(cells_update) = R.cols(cells_update) % (pow((E + 1) / (O + 1), theta) * Phi.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns
    
    // Step 3: put cells back 
    E += sum(R.cols(cells_update), 1) * Pr_b.t();
    O += R.cols(cells_update) * Phi.cols(cells_update).t(); 
    
  }
  return 0;
}


void harmony::moe_correct_ridge_cpp() {
  Z_corr = Z_orig;
  for (int k = 0; k < K; k++) { 
    Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
    W = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();
    W.row(0).zeros(); // do not remove the intercept 
    Z_corr -= W.t() * Phi_Rk;
  }
  Z_cos = arma::normalise(Z_corr, 2, 0);
}


RCPP_MODULE(harmony_module) {
  class_<harmony>("harmony")
  .constructor<int>()
  
  .field("Z_corr", &harmony::Z_corr)  
  .field("Z_orig", &harmony::Z_orig)  
  .field("Z_cos", &harmony::Z_cos)  
  .field("R", &harmony::R)  
  .field("Y", &harmony::Y)  
  .field("Phi", &harmony::Phi)        
  .field("Phi_moe", &harmony::Phi_moe)
  .field("Pr_b", &harmony::Pr_b)    
  .field("objective_kmeans", &harmony::objective_kmeans)
  .field("objective_kmeans_dist", &harmony::objective_kmeans_dist)
  .field("objective_kmeans_entropy", &harmony::objective_kmeans_entropy)
  .field("objective_kmeans_cross", &harmony::objective_kmeans_cross)    
  .field("objective_harmony", &harmony::objective_harmony)
  .field("dist_mat", &harmony::dist_mat)
  .field("ran_setup", &harmony::ran_setup)
  .field("ran_init", &harmony::ran_init)
  
  
  .field("N", &harmony::N)
  .field("K", &harmony::K)
  .field("B", &harmony::B)
  .field("d", &harmony::d)
  .field("W", &harmony::W)
  .field("max_iter_kmeans", &harmony::max_iter_kmeans)
  
  .field("sigma", &harmony::sigma)
  .field("theta", &harmony::theta)
  .field("lambda", &harmony::lambda)
  .field("O", &harmony::O) 
  .field("E", &harmony::E)    
  .field("update_order", &harmony::update_order)    
  .field("cells_update", &harmony::cells_update)    
  .field("kmeans_rounds", &harmony::kmeans_rounds)    
  .field("epsilon_kmeans", &harmony::epsilon_kmeans)    
  .field("epsilon_harmony", &harmony::epsilon_harmony)
  
  // .method("init_cluster", &harmony::init_cluster)
  .method("check_convergence", &harmony::check_convergence)
  .method("setup", &harmony::setup)
  .method("compute_objective", &harmony::compute_objective)
  .method("update_R", &harmony::update_R)
  .method("init_cluster_cpp", &harmony::init_cluster_cpp)
  .method("cluster_cpp", &harmony::cluster_cpp)
  .method("moe_correct_ridge_cpp", &harmony::moe_correct_ridge_cpp)
  
  ;
}








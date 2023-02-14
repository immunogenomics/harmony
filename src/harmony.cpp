#include <algorithm>
#include <chrono>
#include <Rcpp.h>
#include "harmony.h"
#include "types.h"
#include "utils.h"





using namespace std::chrono;
typedef high_resolution_clock::time_point Timepoint;
class Timer{
public:
  Timepoint start;
  std::string task_name;
  double& timer;
  double initial;
  Timer(std::string _task_name, double& _timer): timer(_timer){
    this->task_name = _task_name;
    this->start = high_resolution_clock::now();
  }
  double getLapse(){
    return duration_cast<duration<double>>(high_resolution_clock::now() - this->start).count();
  }
  
  ~Timer(){
    auto time_span = duration_cast<duration<double>>(high_resolution_clock::now() - this->start);
    timer += time_span.count();
    // Rcout << "Task: " << this->task_name << " took " << time_span.count() << " seconds" <<std::endl;
  }
  
  
};




harmony::harmony() :
    window_size(3),
    ran_setup(false),
    ran_init(false),
    verbose(false)
{}



void harmony::setup(MATTYPE& __Z, MATTYPE& __Phi,
                    VECTYPE __sigma, VECTYPE __theta, int __max_iter_kmeans,
                    float __epsilon_kmeans, float __epsilon_harmony,
                    int __K, float tau, float __block_size,
                    MATTYPE __lambda, bool __verbose) {
  
  Z_orig = MATTYPE(__Z);
  Z_cos = arma::normalise(__Z, 2, 0);
  Z_corr = zeros(size(Z_orig));
  
  
  Phi = __Phi;
  Pr_b = sum(Phi, 1);

  // Algorithm constants
  N = Z_corr.n_cols;
  B = Phi.n_rows;
  d = Z_corr.n_rows;

  epsilon_kmeans = __epsilon_kmeans;
  epsilon_harmony = __epsilon_harmony;

  // Hyperparameters
  K = __K;
  lambda = __lambda;
  sigma = __sigma;
  sigma_prior = __sigma;
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
  
  MATTYPE intcpt = ones(1, N);
  Phi_moe = join_cols(intcpt, Phi);
  
  W = zeros<MATTYPE>(B + 1, d);
  Phi_Rk = zeros<MATTYPE>(B + 1, N);  
  
}

void harmony::setY(const MATTYPE& _Y){
  Y = _Y;
  Yset = true;
}


void harmony::init_cluster_cpp(unsigned C) {
  
  if(!Yset){
    Rcerr << "Hard k-means centroids initialization"  <<std::endl;
    Y = kmeans_centers(Z_cos, K).t();
  }else{
    Rcerr << "Ommiting computation hard k-means"  <<std::endl;
  }

  // Cosine normalization of columns
  Y = arma::normalise(Y, 2, 0); 
  
  dist_mat = 2 * (1 - Y.t() * Z_cos);
  R = dist_mat;
  R.each_col() /= sigma;
  R = exp(-dist_mat);
  R.each_row() /= sum(R, 0);
  
  // (2) ASSIGN CLUSTER PROBABILITIES
  // using a nice property of cosine distance,
  // compute squared distance directly with cross product
  dist_mat = 2 * (1 - Y.t() * Z_cos);
  
  // (3) BATCH DIVERSITY STATISTICS
  E = sum(R, 1) * Pr_b.t();
  O = R * Phi.t();
  
  compute_objective();
  objective_harmony.push_back(objective_kmeans.back());
  
  dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed

  ran_init = true;
  
}

void harmony::compute_objective() {
  float kmeans_error = as_scalar(accu(R % dist_mat));
  float _entropy = as_scalar(accu(safe_entropy(R).each_col() % sigma)); // NEW: vector sigma
  float _cross_entropy;
  _cross_entropy = as_scalar(accu((R.each_col() % sigma) % ((arma::repmat(theta.t(), K, 1) % log((O + 1) / (E + 1))) * Phi)));

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
    // p.increment();
    if (Progress::check_abort())
      return(-1);
    
    // STEP 1: Update Y (cluster centroids)
    Y = arma::normalise(Z_cos * R.t(), 2, 0);
    dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed
    Rcout << "Max: " << dist_mat.max() << " Min:" << dist_mat.min() << std::endl;

    auto norm = calculate_norm(Y);
    Rcout << "YMax: " << norm.max() << " YMin:" << norm.min() << ", Size: " << size(norm) << std::endl;
    
    auto norm2 = calculate_norm(Z_cos);
    Rcout << "Z_cosMax: " << norm2.max() << " Z_cosMin:" << norm2.min() << ", Size: " << size(norm2) << std::endl;
    
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

  Rcout << std::endl;
  kmeans_rounds.push_back(iter);
  objective_harmony.push_back(objective_kmeans.back());
  Rcout << "compute_Y(): " << compute_Y_timer<< std::endl;
  Rcout << "compute_objective(): " << compute_objective_timer<< std::endl;
  Rcout << "update_R(): " << update_R_timer<< std::endl;
  return 0;
}




int harmony::update_R() {
  update_order = shuffle(linspace<uvec>(0, N - 1, N));
  _scale_dist = -dist_mat; // K x N
  _scale_dist.each_col() /= sigma; // NEW: vector sigma
  _scale_dist.each_row() -= max(_scale_dist, 0);
  _scale_dist = exp(_scale_dist);

  // GENERAL CASE: online updates, in blocks of size (N * block_size)
  unsigned n_blocks = (int)(ceil(1.0 / block_size));
  unsigned cells_per_block = (N / n_blocks) + 1;
  for (unsigned i = 0; i < n_blocks; i++) {
    // gather cell updates indices
    unsigned idx_min = i * cells_per_block;
    unsigned idx_max = min(idx_min + cells_per_block, N);
    if (idx_min > idx_max) break;
    uvec idx_list = linspace<uvec>(idx_min, idx_max - 1, idx_max - idx_min);
    cells_update = update_order.rows(idx_list);

    // Step 1: remove cells
    E -= sum(R.cols(cells_update), 1) * Pr_b.t();
    O -= R.cols(cells_update) * Phi.cols(cells_update).t();

    // Step 2: recompute R for removed cells
    R.cols(cells_update) = _scale_dist.cols(cells_update);
    R.cols(cells_update) = R.cols(cells_update) % (harmony_pow((E + 1) / (O + 1), theta) * Phi.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns
    
    // Step 3: put cells back 
    E += sum(R.cols(cells_update), 1) * Pr_b.t();
    O += R.cols(cells_update) * Phi.cols(cells_update).t();
    
  }
  return 0;
}


void harmony::moe_correct_ridge_cpp() {

  Z_corr = Z_orig;
  for (unsigned k = 0; k < K; k++) { 
    Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
    W = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();
    W.row(0).zeros(); // do not remove the intercept 
    Z_corr -= W.t() * Phi_Rk;
  }

}

RCPP_MODULE(harmony_module) {
  class_<harmony>("harmony")
      .constructor()
      .method("check_convergence", &harmony::check_convergence)
      .method("setup", &harmony::setup)
      .method("compute_objective", &harmony::compute_objective)
      .method("init_cluster_cpp", &harmony::init_cluster_cpp)
      .method("cluster_cpp", &harmony::cluster_cpp)	  
      .method("moe_correct_ridge_cpp", &harmony::moe_correct_ridge_cpp)
      .method("setY", &harmony::setY)
  ;
}








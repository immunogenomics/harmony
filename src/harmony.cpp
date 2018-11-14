#include "harmony.h"
#include "utils.h"


// NOTE: This is a dummy constructor, needed by Rcpp
harmony::harmony(int __K): K(__K) {}


void harmony::setup(MATTYPE& __Z, MATTYPE& __Phi, VECTYPE __Pr_b,
                float __sigma, VECTYPE __theta, int __max_iter_kmeans, 
                float __epsilon_kmeans, float __epsilon_harmony, bool __correct_with_Zorig,
                float __alpha, int __K, float tau, float __block_size, 
//                ROWVECTYPE& __w, 
                bool __correct_with_cosine, vector<bool> __batch_mask, 
                int __window_size, MATTYPE __lambda, string __correct_mode) {
  
  correct_with_cosine = __correct_with_cosine;
  if (correct_with_cosine)
    cosine_normalize(__Z, 0, false); // normalize columns
  
  Z_corr = MATTYPE(__Z);
  Z_orig = MATTYPE(__Z);

//  Phi_moe = Phi_moe_new;
  Phi = __Phi;
  Phi_moe = ones<MATTYPE>(Phi.n_rows + 1, Phi.n_cols); // same as Phi plus an intercept term
  Phi_moe.rows(1, Phi_moe.n_rows - 1) = Phi;
  N = Z_corr.n_cols;
//  N_b = sum(Phi, 1);
//  Pr_b = N_b / N;
  Pr_b = __Pr_b;
  B = Phi.n_rows;
  d = Z_corr.n_rows; 
//  w = __w;
  window_size = __window_size;
  epsilon_kmeans = __epsilon_kmeans;
  epsilon_harmony = __epsilon_harmony;

  
  lambda = __lambda;
  
  batch_mask = __batch_mask; // For the special case where some batches won't be corrected
  
  correct_mode = __correct_mode;
    
  sigma = __sigma;
  block_size = __block_size;
  K = __K;
  alpha = __alpha;
  max_iter_kmeans = __max_iter_kmeans;
//  converge_thresh = __converge_thresh;
  correct_with_Zorig = __correct_with_Zorig;

  // map from 
//  for (int b = 0; b < B; b++) {
//    phi_map.push_back(find(Phi.row(b) > 0));
//  }
  
  theta = __theta;
//  theta = set_thetas(__theta, tau, N_b);  
  allocate_buffers();
  ran_setup = true;
//  do_merge_R = false; // (EXPERIMENTAL) try to merge redundant clusters?
  init_cluster();  
}


void harmony::allocate_buffers() {
//  mu_k = zeros<MATTYPE>(d, K); 
//  mu_bk = zeros<CUBETYPE>(d, K, B); // nrow, ncol, nslice
//  mu_bk_r = zeros<MATTYPE>(d, N);  
//  mu_k_r = zeros<MATTYPE>(d, N);
  _scale_dist = zeros<MATTYPE>(K, N);    
  __dist = zeros<MATTYPE>(K, N);    
  O = zeros<MATTYPE>(K, B);
  E = zeros<MATTYPE>(K, B);  
  W = zeros<MATTYPE>(B + 1, d); 
//  W = zeros<CUBETYPE>(B + 1, d, K); // (1+B) x d x K
  Phi_Rk = zeros<MATTYPE>(B + 1, N);
}

/*
VECTYPE harmony::set_thetas(float theta_max, float tau, VECTYPE& N_b) {
  VECTYPE res;
  res.set_size(N_b.n_rows);
  if (tau == 0) {
    res.fill(theta_max);
  } else {
    for (int b = 0; b < N_b.n_rows; b++) {
      res.row(b) = theta_max * (1 - exp(-pow(N_b.row(b) / (K * tau), 2)));
    }
  }
  return res;
//  theta.print("theta: ");
}
*/

/*
void harmony::set_R_merge_flag(float merge_thresh_new) {
  merge_thresh_global = merge_thresh_new;
  do_merge_R = true;
}

void harmony::update_R_merge() {
  R = merge_R(R, merge_thresh_global);
  K = R.n_rows;
  // Y will be updated in first round of clustering
  Rcout << "new K: " << K << endl;
}
*/
// BEGIN NUMERICAL METHODS 
void harmony::harmonize(int iter_harmony) {
  int err_status;
  for (int iter = 0; iter < iter_harmony; iter++) {
    std::ostringstream oss;
    oss << "Harmony " << iter + 1 << "/" << iter_harmony;
    Rcpp::Function msg("message"); 
    msg(std::string(oss.str()));
    
    // STEP 1: do clustering
    err_status = cluster();
    if (err_status == -1) {
      Rcout << "terminated by user" << endl;
      break;
    } else if (err_status != 0) {
      break;
    }
    
    // STEP 2: regress out covariates
    moe_correct_ridge();      
    
    // STEP 3: check for convergence
    if (check_convergence(1)) {
      Rcout << "Harmony converged after " << iter + 1 << " iterations\n" << endl;
      break;
    }    
  }
}




void harmony::init_cluster() {
  if (!ran_setup) {
    Rcout << "ERROR: before initializing cluster, run setup" << endl;
    return;
  }
  // TODO: kmeans++, then average for nearest neighboring centroids across batches
  // get K random points from each batch
  // average the points   
  Y.set_size(d, K);
  Y.fill(0);
  for (int b = 0; b < B; b++) {
    uvec q = find(Phi.row(b) > 0); // indices of cells belonging to batch (b)
//    Rcout << "batch " << b << ": " << q.n_elem << endl;
    uvec rand_idx = conv_to<uvec>::from(randi(K, distr_param(0, q.n_elem - 1)));
    Y += Z_corr.cols(q.elem(rand_idx));
  }
  Y /= B;  
  cosine_normalize(Y, 0, false); // normalize columns

  
//  if (correct_with_Zorig)
  Z_cos = MATTYPE(Z_orig);
//  else 
//    Z_cos = MATTYPE(Z_corr);
  cosine_normalize(Z_cos, 0, true); // normalize columns
  __dist = 2 * (1 - Y.t() * Z_cos); // initial estimate based on Y_0 and Z_0
  // using a nice property of cosine distance,
  // compute squared distance directly with cross product
  R = - (1 / sigma) * 2 * (1 - (Y.t() * Z_cos));  
  R.each_row() -= max(R, 0);  
  R = exp(R);
  R.each_row() /= sum(R, 0);

  E = sum(R, 1) * Pr_b.t();
  O = R * Phi.t();
  
  compute_objective();
  objective_harmony.push_back(objective_kmeans.back());
  ran_init = true;
  
}



// TODO: generalize to adaptive sigma values
// TODO: use cached distance computation from before
void harmony::compute_objective() {
  
//  float kmeans_error = as_scalar(accu((R.each_row() % w) % (2 * (1 - (Y.t() * Z_cos)))));
//  float _entropy = as_scalar(accu(safe_entropy(R).each_row() % w));  
  
//  float kmeans_error = as_scalar(accu(R % (2 * (1 - (Y.t() * Z_cos)))));
  float kmeans_error = as_scalar(accu(R % __dist));
  float _entropy = as_scalar(accu(safe_entropy(R)));
  
//  float _cross_entropy = as_scalar(accu(R % log((E / O) * Phi)));   
//  objective_kmeans.push_back(kmeans_error + sigma * _entropy +
//                      sigma * theta * _cross_entropy);

  float _cross_entropy;
  dir_prior = alpha * E; // here, alpha is in [0, Inf). Reflects strength of dirichlet prior. 
//  _cross_entropy = as_scalar(accu((R.each_row() % w) % ((arma::repmat(theta.t(), K, 1) % log((O + dir_prior) / (E + dir_prior))) * Phi)));
  _cross_entropy = as_scalar(accu(R % ((arma::repmat(theta.t(), K, 1) % log((O + dir_prior) / (E + dir_prior))) * Phi)));
  
//    _cross_entropy = as_scalar(accu((R.each_row() % w) % ((arma::repmat(theta.t(), K, 1) % log((O + alpha) / (E + alpha))) * Phi))); 
  
  
  objective_kmeans.push_back(kmeans_error + sigma * _entropy + sigma * _cross_entropy);
  objective_kmeans_dist.push_back(kmeans_error);
  objective_kmeans_entropy.push_back(sigma * _entropy); 
  objective_kmeans_cross.push_back(sigma * _cross_entropy);
  
//  Rcout << "OBJ: " << kmeans_error + sigma * _entropy + sigma * _cross_entropy << endl;
  
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
  
}



int harmony::cluster() {
  if (!ran_setup) {
    Rcout << "ERROR: before clustering, run init_cluster" << endl;
    return -1;
  }
  int err_status = 0;
  int iter; 
  Progress p(max_iter_kmeans, true);
  
  __dist = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed
  for (iter = 0; iter < max_iter_kmeans; iter++) {    
    p.increment();
    if (Progress::check_abort())
      return(-1);

    // STEP 1: Update Y
//    Y = Z_cos * (R.each_row() % w).t();
    Y = Z_cos * R.t();
    cosine_normalize(Y, 0, true);
    __dist = 2 * (1 - Y.t() * Z_cos); // Y was changed

    // STEP 2: Update R
    err_status = compute_R();
    if (err_status != 0) {
      Rcout << "Compute R failed. Exiting from clustering." << endl;
      return err_status;
    }

    // STEP 3: Check for convergence
    compute_objective();
    if (iter > window_size) {
      bool convergence_status = check_convergence(0); 
      if (convergence_status) {
//        Rcout << "... Breaking Clustering ..., status = " << convergence_status << endl;
        iter++;
        Rcout << "Clustered for " << iter << " iterations" << endl;
        break;        
      }
    }
  }
  kmeans_rounds.push_back(iter);
  objective_harmony.push_back(objective_kmeans.back());
  return 0;
}



void harmony::foo1() {
}

void harmony::foo2() {
}

void harmony::foo3() {
}

void harmony::foo4() {

}
void harmony::foo5() {

}
void harmony::foo6() {
}
void harmony::foo7() {
}
void harmony::foo8() {
}
void harmony::foo9() {

}
void harmony::foo10() {
  
}


int harmony::compute_R() {  
  update_order = shuffle(linspace<uvec>(0, N - 1, N));
// compute distance based scale of R
//  _scale_dist = - (1 / sigma) * 2 * (1 - Y.t() * Z_cos);  
  _scale_dist = - (1 / sigma) * __dist;
  _scale_dist.each_row() -= max(_scale_dist, 0);
  _scale_dist = exp(_scale_dist);

  // GENERAL CASE: online updates, in blocks of size (N * block_size)
  for (int i = 0; i < ceil(1. / block_size); i++) {
    // gather cell updates indices
    int idx_min = i * N * block_size;
    int idx_max = min((int)((i + 1) * N * block_size), N - 1);
    if (idx_min > idx_max) break; // TODO: fix the loop logic so that this never happens
    
    uvec idx_list = linspace<uvec>(0, idx_max - idx_min, idx_max - idx_min + 1);
    cells_update = update_order.rows(idx_list); // FOR DEBUGGING: using global cells_update vector    

    // Step 1: remove cells
    E -= sum(R.cols(cells_update), 1) * Pr_b.t();
    O -= R.cols(cells_update) * Phi.cols(cells_update).t();
    
    // Step 2: recompute R for cells
    R.cols(cells_update) = _scale_dist.cols(cells_update);    
    dir_prior = alpha * E;
    R.cols(cells_update) = R.cols(cells_update) % (pow((E + dir_prior) / (O + dir_prior), theta) * Phi.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns

    // Step 3: put cells back 
    E += sum(R.cols(cells_update), 1) * Pr_b.t();
    O += R.cols(cells_update) * Phi.cols(cells_update).t();   
    
  }
  return 0;
}


void harmony::moe_correct_ridge() {
  Z_corr = Z_orig;
  for (int k = 0; k < K; k++) { 
    Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
    W = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();    
    // do not remove the intercept 
    W.row(0).zeros(); 
    Z_corr -= W.t() * Phi_Rk;
  }
  
  Z_cos = MATTYPE(Z_corr);
  cosine_normalize(Z_cos, 0, true); // normalize columns
}



/*
void harmony::moe_correct_ridge() {
  // here, we model both an intercept and one-hot 
  // we overcome singularity by putting a ridge-regression penalty on the non-intercept terms
//  W = zeros<CUBETYPE>(Phi_moe.n_rows, Z_orig.n_rows, K); // (1+B) x d x K
  for (int k = 0; k < K; k++) { 
    W.slice(k) = inv(Phi_moe * diagmat(R.row(k)) * Phi_moe.t() + lambda) * Phi_moe * diagmat(R.row(k)) * Z_orig.t();
  }

  // (3) remove batch effects from each row, one at a time 
  for (int d = 0; d < Z_orig.n_rows; d++) {
    MATTYPE W_sub = W.subcube(1, d, 0, W.n_rows - 1, d, W.n_slices - 1); // can we just use .col(d)? 
    Z_corr.row(d) = Z_orig.row(d) - sum(R % (W_sub.t() * Phi_moe.rows(1, Phi_moe.n_rows - 1)), 0);
  }

  Z_cos = MATTYPE(Z_corr);
  cosine_normalize(Z_cos, 0, true); // normalize columns
}
*/

/*
void harmony::moe_correct_onehot() {
  // (1) constrain the means of every cell to the batch-agnostic centroids
//  mu_k = Z_orig * (R.each_row() % w).t(); // d x K   
//  mu_k *= diagmat(1 / sum((R.each_row() % w), 1)); // divide by the effective number of cells in each cluster
  mu_k = Z_orig * R.t(); // d x K   
  mu_k *= diagmat(1 / sum(R, 1)); // divide by the effective number of cells in each cluster
  mu_k_r = mu_k * R; // expected location of each cell (by cluster mixture) d * N    
  
  
  // (2) model the differences of y = Z - mu with batch variables and residuals
  MATTYPE y = Z_orig - mu_k_r;
  W = zeros<CUBETYPE>(Phi.n_rows, Z_orig.n_rows, K); // B x d x K
  for (int k = 0; k < K; k++) { 
    W.slice(k) = inv(Phi * diagmat(R.row(k)) * Phi.t()) * Phi * diagmat(R.row(k)) * y.t();
  }

  // (3) remove batch effects from each row, one at a time 
  for (int d = 0; d < Z_orig.n_rows; d++) {
    MATTYPE W_sub = W.subcube(0, d, 0, W.n_rows - 1, d, W.n_slices - 1); // can we just use .col(d)? 
    Z_corr.row(d) = Z_orig.row(d) - sum(R % (W_sub.t() * Phi), 0);
  }

  Z_cos = MATTYPE(Z_corr);
  cosine_normalize(Z_cos, 0, true); // normalize columns
}
*/


/*
void harmony::moe_correct_contrast() {  
  // (1) model the Z directly
  W = zeros<CUBETYPE>(Phi_moe.n_rows, Z_orig.n_rows, K); // B x d x K
  for (int k = 0; k < K; k++) { 
    W.slice(k) = inv(Phi_moe * diagmat(R.row(k)) * Phi_moe.t()) * Phi_moe * diagmat(R.row(k)) * Z_orig.t();
  }

  // (2) remove batch effects from each row, one at a time 
  // NOTE: first row of W is intercept, which we do not remove as batch effect
  for (int d = 0; d < Z_orig.n_rows; d++) {
    MATTYPE W_sub = W.subcube(1, d, 0, W.n_rows - 1, d, W.n_slices - 1);
    Z_corr.row(d) = Z_orig.row(d) - sum(R % (W_sub.t() * Phi_moe.rows(1, Phi_moe.n_rows - 1)), 0);
  }

  Z_cos = MATTYPE(Z_corr);
  cosine_normalize(Z_cos, 0, true); // normalize columns  
}
*/




/*    
// This method operates in Euclidean PCA space (no cosine normalization)
void harmony::gmm_correct_armadillo() {  
  if (correct_with_Zorig) {
//    mu_k = Z_orig * (R.each_row() % w).t(); // d x K    
    mu_k = Z_orig * R.t(); // d x K    
  } else {
//    mu_k = Z_corr * (R.each_row() % w).t(); // d x K    
    mu_k = Z_corr * R.t(); // d x K    
  }
  
//  mu_k *= diagmat(1 / sum((R.each_row() % w), 1)); // divide by the effective number of cells in each cluster
  mu_k *= diagmat(1 / sum(R, 1)); // divide by the effective number of cells in each cluster
  mu_k_r = mu_k * R; // expected location of each cell (by cluster mixture) d * N    

  // Assumes that all cells within a batch have equal weight
  for (int b = 0; b < Phi.n_rows; b++) {
    uvec q = find(Phi.row(b) == 1); 
    
    // Asymmetric Harmony: don't correct this batch
    if (batch_mask[b] == false) {
      mu_bk_r.cols(q).zeros();
      mu_k_r.cols(q).zeros();
      continue;
    }
        
    if (correct_with_Zorig)
      mu_bk.slice(b) = Z_orig.cols(q) * R.cols(q).t();
    else 
      mu_bk.slice(b) = Z_corr.cols(q) * R.cols(q).t(); 
    
    mu_bk.slice(b) = mu_bk.slice(b) * diagmat(1 / sum(R.cols(q), 1));    
    mu_bk_r.cols(q) = mu_bk.slice(b) * R.cols(q);
  }
  
  if (correct_with_Zorig) 
    Z_corr = Z_orig - mu_bk_r + mu_k_r;
  else 
    Z_corr = Z_corr - mu_bk_r + mu_k_r;
    
  Z_cos = MATTYPE(Z_corr);
  cosine_normalize(Z_cos, 0, true); // normalize columns  
}
*/


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
    
    
  .field("N", &harmony::N)
  .field("K", &harmony::K)
  .field("B", &harmony::B)
  .field("d", &harmony::d)
//  .field("w", &harmony::w)
  .field("W", &harmony::W)
    

  .field("max_iter_kmeans", &harmony::max_iter_kmeans)

  .field("sigma", &harmony::sigma)
  .field("theta", &harmony::theta)
  .field("alpha", &harmony::alpha)
  .field("lambda", &harmony::lambda)
  .field("O", &harmony::O) 
  .field("E", &harmony::E)    
  .field("update_order", &harmony::update_order)    
  .field("cells_update", &harmony::cells_update)    
  .field("kmeans_rounds", &harmony::kmeans_rounds)    
  .field("epsilon_kmeans", &harmony::epsilon_kmeans)    
  .field("epsilon_harmony", &harmony::epsilon_harmony)    


    
  .method("harmonize", &harmony::harmonize)
  .method("init_cluster", &harmony::init_cluster)
  .method("check_convergence", &harmony::check_convergence)
  .method("setup", &harmony::setup)
//  .method("set_thetas", &harmony::set_thetas)
//  .method("set_R_merge_flag", &harmony::set_R_merge_flag)
//  .method("update_R_merge", &harmony::update_R_merge)
  .method("cluster", &harmony::cluster)
//  .method("gmm_correct_armadillo", &harmony::gmm_correct_armadillo)   
//  .method("moe_correct_onehot", &harmony::moe_correct_onehot)   
  .method("moe_correct_ridge", &harmony::moe_correct_ridge)   
//  .method("moe_correct_ridge2", &harmony::moe_correct_ridge2)   
//  .method("moe_correct_ridge3", &harmony::moe_correct_ridge3)
  .method("compute_objective", &harmony::compute_objective)
  .method("compute_R", &harmony::compute_R)

  .method("foo1", &harmony::foo1)
  .method("foo2", &harmony::foo2)
  .method("foo3", &harmony::foo3)
  .method("foo4", &harmony::foo4)
  .method("foo5", &harmony::foo5)
  .method("foo6", &harmony::foo6)
  .method("foo7", &harmony::foo7)
  .method("foo8", &harmony::foo8)
  .method("foo9", &harmony::foo9)
  .method("foo10", &harmony::foo10)    
    
  ;
}








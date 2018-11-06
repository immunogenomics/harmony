#include "harmony.h"
#include "utils.h"


// NOTE: This is a dummy constructor, needed by Rcpp
harmony::harmony(int __K): K(__K) {}


void harmony::setup(mat& Z_new, mat& Phi_new, 
                        float __sigma, float __theta, int __max_iter_kmeans, 
                        float __epsilon_kmeans, float __epsilon_harmony, bool __correct_with_Zorig,
                        float __alpha, int __K, float tau, float __block_size, 
                        rowvec& w_new, bool __correct_with_cosine, vector<bool> batch_mask_new, int __window_size) {
  
  correct_with_cosine = __correct_with_cosine;
  if (correct_with_cosine)
    cosine_normalize(Z_new, 0, false); // normalize columns
  
  Z_corr = mat(Z_new);
  Z_orig = mat(Z_new);

  Phi = Phi_new;
  N = Z_corr.n_cols;
  N_b = sum(Phi, 1);
  Pr_b = N_b / N;
  B = Phi.n_rows;
  d = Z_corr.n_rows; 
  w = w_new;
  window_size = __window_size;
  epsilon_kmeans = __epsilon_kmeans;
  epsilon_harmony = __epsilon_harmony;
  
  batch_mask = batch_mask_new; // For the special case where some batches won't be corrected
  
  sigma = __sigma;
  block_size = __block_size;
  K = __K;
  alpha = __alpha;
  max_iter_kmeans = __max_iter_kmeans;
//  converge_thresh = __converge_thresh;
  correct_with_Zorig = __correct_with_Zorig;

  // map from 
  for (int b = 0; b < B; b++) {
    phi_map.push_back(find(Phi.row(b) > 0));
  }  
  
  theta = set_thetas(__theta, tau, N_b);  
  allocate_buffers();
  ran_setup = true;
  do_merge_R = false; // (EXPERIMENTAL) try to merge redundant clusters?
  do_theta2 = false; // (EXPERIMENTAL) initially this is false
  init_cluster();  
}

void harmony::setup_batch2(mat& Phi2_new, float theta2_new, float tau) {
  if (theta2_new == 0) return;
  
  Phi2 = Phi2_new;
  N_b2 = sum(Phi2, 1);
  Pr_b2 = N_b2 / N;
  B2 = Phi.n_rows;
  theta2 = set_thetas(theta2_new, tau, N_b2);    
  E2 = sum(R, 1) * Pr_b2.t();
  O2 = R * Phi2.t();
  do_theta2 = true;
  
}


void harmony::allocate_buffers() {
  mu_k = zeros<mat>(d, K); 
  mu_bk = zeros<cube>(d, K, B); // nrow, ncol, nslice
  mu_bk_r = zeros<mat>(d, N);  
  mu_k_r = zeros<mat>(d, N);
  _scale_dist = zeros<mat>(K, N);    
  O = zeros<mat>(K, B);
  E = zeros<mat>(K, B);  
}


vec harmony::set_thetas(float theta_max, float tau, vec& N_b) {
  vec res;
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

// BEGIN NUMERICAL METHODS 
void harmony::harmonize(int iter_harmony) {
  int err_status;
  for (int iter = 0; iter < iter_harmony; iter++) {
    std::ostringstream oss;
    oss << "Harmony " << iter + 1 << "/" << iter_harmony;
    Rcpp::Function msg("message"); 
    msg(std::string(oss.str()));
    
    err_status = cluster();
    if (err_status == -1) {
      Rcout << "terminated by user" << endl;
      break;
    } else if (err_status != 0) {
      break;
    }
    
    gmm_correct_armadillo();

    // NOTE: this does not work. For now, run all iterations. 
    if (check_convergence(1)) {
      Rcout << "Harmony converged after " << iter + 1 << " iterations\n" << endl;
//      Rprintf("Converged after %d iterations\n", iter);
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

  
  if (correct_with_Zorig)
    Z_cos = mat(Z_orig);
  else 
    Z_cos = mat(Z_corr);
  cosine_normalize(Z_cos, 0, true); // normalize columns
  
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


/*
// OPTIONAL: create batch specific covariates
//           to preserve structure inside batches
void harmony::init_batch_clusters(uvec & batches, float merge_thresh,
                              float sigma_local, int K_local) {  
  // TODO: check if batch is 1-indexed. 
  
  compute_phi_hat(batches, merge_thresh, sigma_local, K_local);
  N_Kb = sum(phi_hat, 1);
  Pr_Kb = N_Kb / N;
  E2 = sum(R, 1) * Pr_Kb.t();
  O2 = R * phi_hat.t();
  
  // TODO: figure out if theta2 needs to be scaled on cluster size
  theta2 = zeros<vec>(N_Kb.n_rows);
  int i = 0;
  for (int b = 0; b < B; b++) {
    for (int k = 0; k < Kb[b]; k++) {
      theta2(i) = theta[b] * (N_Kb(i) / N_b(b));
      i++;
    }
  }  
  do_conservation = true;

}
*/


/*

void harmony::compute_phi_hat(const uvec & batches, float merge_thresh,
                              float sigma_local, int K_local) {
  mat X = mat(Z_orig);
  cosine_normalize(X, 2, true);
  
  // (1) list of per-batch cluster matrices
//  vector<mat> R_list;
  unsigned Kb_total = 0;
  for (int b = 0; b < B; b++) {
    uvec idx = find(batches == b); 
    R_list.push_back(fuzzy_kmeans(X.cols(idx), sigma_local, K_local, max_iter_kmeans, epsilon_kmeans));
    R_list[b] = merge_R(R_list[b], merge_thresh);
    Kb_total += R_list[b].n_rows;
//    Rcout << R_list[b].n_rows << endl;  
    Kb.push_back(R_list[b].n_rows);
  }

  // (2) merge them into a single sparse matrix
  int N = X.n_cols;  
  phi_hat = zeros<mat>(Kb_total, N);
  unsigned offset = 0;
//  Rcout << phi_hat.n_rows << " " << phi_hat.n_cols << endl;
  for (int b = 0; b < B; b++) {
    uvec c_idx = find(batches == b); 
    uvec r_idx = linspace<uvec>(offset, offset + R_list[b].n_rows - 1, R_list[b].n_rows);
    phi_hat(r_idx, c_idx) = R_list[b];                
    offset += R_list[b].n_rows;
  } 
//  Rcout << phi_hat.n_rows << " " << phi_hat.n_cols << endl;    
}
*/

// TODO: generalize to adaptive sigma values
// TODO: use cached distance computation from before
void harmony::compute_objective() {
  
  float kmeans_error = as_scalar(accu((R.each_row() % w) % (2 * (1 - (Y.t() * Z_cos)))));
  float _entropy = as_scalar(accu(safe_entropy(R).each_row() % w));  
  
//  float _cross_entropy = as_scalar(accu(R % log((E / O) * Phi)));   
//  objective_kmeans.push_back(kmeans_error + sigma * _entropy +
//                      sigma * theta * _cross_entropy);

  float _cross_entropy;
  dir_prior = alpha * E; // here, alpha is in [0, Inf). Reflects strength of dirichlet prior. 
  _cross_entropy = as_scalar(accu((R.each_row() % w) % ((arma::repmat(theta.t(), K, 1) % log((O + dir_prior) / (E + dir_prior))) * Phi)));

  if (do_theta2) {
    dir_prior2 = alpha * E2; // here, alpha is in [0, Inf). Reflects strength of dirichlet prior. 
    _cross_entropy =+ as_scalar(accu((R.each_row() % w) % ((arma::repmat(theta2.t(), K, 1) % log((O2 + dir_prior2) / (E2 + dir_prior2))) * Phi2)));
    
  }
  
  
//    _cross_entropy = as_scalar(accu((R.each_row() % w) % ((arma::repmat(theta.t(), K, 1) % log((O + alpha) / (E + alpha))) * Phi)));    
  
  
  objective_kmeans.push_back(kmeans_error + sigma * _entropy +
                      sigma * _cross_entropy);

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
  Progress p(max_iter_kmeans, true);
  
  int iter; 
//  Rcout << "max iter kmeans: " << max_iter_kmeans << endl;
  for (iter = 0; iter < max_iter_kmeans; iter++) {    
    p.increment();
    if (Progress::check_abort())
      return(-1);
  
//    Y = Z_cos * R.t();
    Y = Z_cos * (R.each_row() % w).t();
    cosine_normalize(Y, 2, true);
    err_status = compute_R();
    if (err_status != 0) {
      Rcout << "Compute R failed. Exiting from clustering." << endl;
      return err_status;
    }
    
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

int harmony::compute_R() {  
  update_order = shuffle(linspace<uvec>(0, N - 1, N));

  // compute distance based scale of R
  _scale_dist = - (1 / sigma) * 2 * (1 - Y.t() * Z_cos);  
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
    
    
    E -= sum(R.cols(cells_update), 1) * Pr_b.t();
    O -= R.cols(cells_update) * Phi.cols(cells_update).t();
    
    if (do_theta2) {
      E2 -= sum(R.cols(cells_update), 1) * Pr_b2.t();
      O2 -= R.cols(cells_update) * Phi2.cols(cells_update).t();      
    }
    
    R.cols(cells_update) = _scale_dist.cols(cells_update);
    
    dir_prior = alpha * E;
    R.cols(cells_update) = R.cols(cells_update) % (pow((E + dir_prior) / (O + dir_prior), theta) * Phi.cols(cells_update));
    if (do_theta2) {
      dir_prior2 = alpha * E2;
      R.cols(cells_update) = R.cols(cells_update) % (pow((E2 + dir_prior2) / (O2 + dir_prior2), theta2) * Phi2.cols(cells_update));
    }
    
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns

    E += sum(R.cols(cells_update), 1) * Pr_b.t();
    O += R.cols(cells_update) * Phi.cols(cells_update).t();   

    if (do_theta2) {
      E2 += sum(R.cols(cells_update), 1) * Pr_b2.t();
      O2 += R.cols(cells_update) * Phi2.cols(cells_update).t();      
    }
    
  }
  return 0;
}



// This method operates in Euclidean PCA space (no cosine normalization)
void harmony::gmm_correct_armadillo() {  
  if (correct_with_Zorig) {
    mu_k = Z_orig * (R.each_row() % w).t(); // d x K    
  } else {
    mu_k = Z_corr * (R.each_row() % w).t(); // d x K    
  }
  
  mu_k *= diagmat(1 / sum((R.each_row() % w), 1)); // divide by the effective number of cells in each cluster
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
    
  Z_cos = mat(Z_corr);
  cosine_normalize(Z_cos, 0, true); // normalize columns  
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
  .field("Phi2", &harmony::Phi2)    
  .field("Pr_b", &harmony::Pr_b)    
  .field("Pr_b2", &harmony::Pr_b2)    
  .field("objective_kmeans", &harmony::objective_kmeans)
  .field("objective_kmeans_dist", &harmony::objective_kmeans_dist)
  .field("objective_kmeans_entropy", &harmony::objective_kmeans_entropy)
  .field("objective_kmeans_cross", &harmony::objective_kmeans_cross)    
  .field("objective_harmony", &harmony::objective_harmony)
    
    
  .field("N", &harmony::N)
  .field("K", &harmony::K)
  .field("B", &harmony::B)
  .field("B2", &harmony::B2)
  .field("d", &harmony::d)
  .field("w", &harmony::w)
    

  .field("max_iter_kmeans", &harmony::max_iter_kmeans)

  .field("sigma", &harmony::sigma)
  .field("theta", &harmony::theta)
  .field("theta2", &harmony::theta2)
  .field("alpha", &harmony::alpha)
  .field("O", &harmony::O) 
  .field("E", &harmony::E)    
  .field("O2", &harmony::O2)    
  .field("E2", &harmony::E2)    
  .field("R_list", &harmony::R_list)    
  .field("update_order", &harmony::update_order)    
  .field("cells_update", &harmony::cells_update)    
  .field("kmeans_rounds", &harmony::kmeans_rounds)    
  .field("epsilon_kmeans", &harmony::epsilon_kmeans)    
  .field("epsilon_harmony", &harmony::epsilon_harmony)    


    
  .method("harmonize", &harmony::harmonize)
  .method("init_cluster", &harmony::init_cluster)
  .method("setup_batch2", &harmony::setup_batch2)
  .method("check_convergence", &harmony::check_convergence)
  .method("setup", &harmony::setup)
  .method("set_thetas", &harmony::set_thetas)
  .method("set_R_merge_flag", &harmony::set_R_merge_flag)
  .method("update_R_merge", &harmony::update_R_merge)
  .method("cluster", &harmony::cluster)
  .method("gmm_correct_armadillo", &harmony::gmm_correct_armadillo)   
  .method("compute_objective", &harmony::compute_objective)
  .method("compute_R", &harmony::compute_R)

    
    
  ;
}








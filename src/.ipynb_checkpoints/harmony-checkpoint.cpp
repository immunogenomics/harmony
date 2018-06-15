#include "harmony.h"
#include "utils.h"


// NOTE: This is a dummy constructor, needed by Rcpp
harmony::harmony(int __K): K(__K) {}


void harmony::setup(arma::mat& Z_new, arma::mat& Phi_new, 
                        double __sigma, double __theta, int __max_iter_kmeans, 
                        float __converge_thresh, bool __correct_with_Zorig,
                        double __alpha, int __K, float tau, float __block_size) {
  
  // unfortunately, need up to 3 copies of Z
  Z_corr = mat(Z_new);
  if (__correct_with_Zorig)
    Z_orig = mat(Z_new);
   
  Phi = Phi_new;
  N = Z_corr.n_cols;
  N_b = sum(Phi, 1);
  Pr_b = N_b / N;
  B = Phi.n_rows;
  d = Z_corr.n_rows; 
  
  sigma = __sigma;
  block_size = __block_size;
  K = __K;
  alpha = __alpha;
  max_iter_kmeans = __max_iter_kmeans;
  converge_thresh = __converge_thresh;
  correct_with_Zorig = __correct_with_Zorig;

  set_thetas(__theta, tau);  
  allocate_buffers();
  ran_setup = true;
  init_cluster();  
}

void harmony::allocate_buffers() {
  mu_k = zeros(d, K); 
  mu_bk = zeros<cube>(d, K, B); // nrow, ncol, nslice
  mu_bk_r = zeros<mat>(d, N);  
  mu_k_r = zeros<mat>(d, N);
  _scale_dist = zeros<mat>(K, N);    
  O = zeros<mat>(K, B);
  E = zeros<mat>(K, B);  
}


void harmony::set_thetas(float theta_max, float tau) {
  theta.set_size(B);
  if (tau == 0) {
    theta.fill(theta_max);
  } else {
    for (int b = 0; b < B; b++) {
      theta.row(b) = theta_max * (1 - exp(-pow(N_b.row(b) / (K * tau), 2)));
    }
  }  
  theta.print("theta: ");
}

/* BEGIN NUMERICAL METHODS */
void harmony::harmonize(int iter_harmony) {
  int err_status;
  for (int iter = 0; iter < iter_harmony; iter++) {
    Rcout << iter + 1 << "/" << iter_harmony << endl;
    err_status = cluster();
    if (err_status == -1) {
      Rcout << "terminated by user" << endl;
      break;
    };
    gmm_correct_armadillo();

    /* NOTE: this does not work. For now, run all iterations. 
    if (check_convergence(1)) {
      Rprintf("Converged after %d iterations\n", iter);
      break;
    } */   
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
    Rcout << "batch " << b << ": " << q.n_elem << endl;
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




// TODO: generalize to adaptive sigma values
// TODO: use cached distance computation from before
void harmony::compute_objective() {
  
  float kmeans_error = as_scalar(accu(R % (2 * (1 - (Y.t() * Z_cos)))));
  float _entropy = as_scalar(accu(safe_entropy(R)));  
  
//  float _cross_entropy = as_scalar(accu(R % log((E / O) * Phi)));   
//  objective_kmeans.push_back(kmeans_error + sigma * _entropy +
//                      sigma * theta * _cross_entropy);
  
  float _cross_entropy = as_scalar(accu(R % ((arma::repmat(theta.t(), K, 1) % log(E / O)) * Phi)));   
  objective_kmeans.push_back(kmeans_error + sigma * _entropy +
                      sigma * _cross_entropy);
  
  
//  if (verbose > 0) {
//    Rprintf("kmeans %0.2f; sig_entropy %0.2f; theta cross %0.2f; cross %0.2f; obj %0.2f\n", 
//           kmeans_error, sigma * _entropy, sigma * theta * _cross_entropy, _cross_entropy,
//           *objective_kmeans.end());
//  }
}



bool harmony::check_convergence(int type) {
  float obj_old, obj_new;
  switch (type) {
    case 0: 
      obj_old = objective_kmeans[objective_kmeans.size() - 2];
      obj_new = objective_kmeans[objective_kmeans.size() - 1];
    case 1:
      obj_old = objective_harmony[objective_harmony.size() - 2];
      obj_new = objective_harmony[objective_harmony.size() - 1];
  }
  
  float obj_change = -(obj_new - obj_old) / abs(obj_old);
//  float obj_change = abs((obj_new - obj_old) / obj_old);
  if (obj_change < converge_thresh) {
    Rcout << "obj conv with " << obj_change << " < " << converge_thresh << endl;
    Rcout << "OLD: " << obj_old << ", NEW: " << obj_new << endl;
    return(true);    
  }
  return(false);
}



int harmony::cluster() {
  if (!ran_setup) {
    Rcout << "ERROR: before clustering, run init_cluster" << endl;
    return(-1);
  }

  Progress p(max_iter_kmeans, true);
  for (int iter = 0; iter < max_iter_kmeans; iter++) {
    p.increment();
    if (Progress::check_abort()) {
      return(-1);
    }
    Y = Z_cos * R.t();
    cosine_normalize(Y, 2, true);  
    compute_R();
    compute_objective();
    if (check_convergence(0)) {
//      printf("Converged after %d iterations\n", iter);
      break;
    }
  }
  objective_harmony.push_back(objective_kmeans.back());
  return(0);
}

void harmony::compute_R() {
  uvec update_order = shuffle(linspace<uvec>(0, N - 1, N));
    
  // compute distance based scale of R
  _scale_dist = - (1 / sigma) * 2 * (1 - Y.t() * Z_cos);  
  _scale_dist.each_row() -= max(_scale_dist, 0);
  _scale_dist = exp(_scale_dist);
    
  for (int i = 0; i < ceil(1. / block_size); i++) {
    // gather cell updates indices
    int idx_min = i * N * block_size;
    int idx_max = min((int)((i + 1) * N * block_size), N - 1);
    uvec idx_list = linspace<uvec>(0, idx_max - idx_min, idx_max - idx_min + 1);
    uvec cells_update = update_order.rows(idx_list);
    
    E -= sum(R.cols(cells_update), 1) * Pr_b.t();
    O -= R.cols(cells_update) * Phi.cols(cells_update).t();

//    R.cols(cells_update) = (pow((E / O) * Phi.cols(cells_update), theta)) % 
//                              _scale_dist.cols(cells_update);
    R.cols(cells_update) = (pow(E / O, theta) * Phi.cols(cells_update)) % 
                              _scale_dist.cols(cells_update);
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0);
    
    E += sum(R.cols(cells_update), 1) * Pr_b.t();
    O += R.cols(cells_update) * Phi.cols(cells_update).t();
    
  }
}



// This method operates in Euclidean PCA space (no cosine normalization)
void harmony::gmm_correct_armadillo() {  
  if (correct_with_Zorig) 
    mu_k = Z_orig * R.t(); // d x K
  else 
    mu_k = Z_corr * R.t(); // d x K
  mu_k *= diagmat(1 / sum(R, 1));
  mu_k_r = mu_k * R; // d * N  
  
  for (int b = 0; b < Phi.n_rows; b++) {
    uvec q = find(Phi.row(b) == 1); 
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
  .field("Pr_b", &harmony::Pr_b)    
  .field("objective_kmeans", &harmony::objective_kmeans)
  .field("objective_harmony", &harmony::objective_harmony)

  .field("N", &harmony::N)
  .field("K", &harmony::K)
  .field("B", &harmony::B)
  .field("d", &harmony::d)
  .field("max_iter_kmeans", &harmony::max_iter_kmeans)

  .field("sigma", &harmony::sigma)
  .field("theta", &harmony::theta)
  .field("alpha", &harmony::alpha)    

  .method("harmonize", &harmony::harmonize)
  .method("init_cluster", &harmony::init_cluster)
  .method("check_convergence", &harmony::check_convergence)
  .method("setup", &harmony::setup)
  .method("set_thetas", &harmony::setup)
  .method("cluster", &harmony::cluster)
  .method("gmm_correct_armadillo", &harmony::gmm_correct_armadillo)   
    
  ;
}








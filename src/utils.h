void cosine_normalize(MATTYPE& X, int margin, bool do_safe) {
  // to avoid Inf values, first divide by max 
  if (margin == 1) {
    for (unsigned r = 0; r < X.n_rows; r++) {
      if (do_safe)
        X.row(r) = X.row(r) / X.row(r).max(); 
      X.row(r) = X.row(r) / norm(X.row(r), 2);
    }     
  } else {
    for (unsigned c = 0; c < X.n_cols; c++) {
      if (do_safe)
        X.col(c) = X.col(c) / X.col(c).max(); 
      X.col(c) = X.col(c) / norm(X.col(c), 2);
    }    
  } 
}


MATTYPE safe_entropy(const MATTYPE& X) {
  MATTYPE A = X % log(X);
  A.elem(find_nonfinite(A)).zeros();
  return(A);
}


// Overload pow to work on a MATTYPErix and vector
MATTYPE pow(MATTYPE A, const VECTYPE & T) {
  for (unsigned c = 0; c < A.n_cols; c++) {
    A.col(c) = pow(A.col(c), as_scalar(T.row(c)));
  }
  return(A);
}


MATTYPE merge_R(const MATTYPE & R, float thresh = 0.8) {
  MATTYPE cor_res = cor(R.t());
  int K = R.n_rows;
  
  // define equivalence classes
  uvec equiv_classes = linspace<uvec>(0, K - 1, K);
  int new_idx;
  for (int i = 0; i < K - 1; i++) {
    for (int j = i + 1; j < K; j++) {
      if (cor_res(i, j) > thresh) {
        new_idx = min(equiv_classes(i), equiv_classes(j));
        equiv_classes(i) = new_idx;
        equiv_classes(j) = new_idx;
      }
    }
  }
  
  // sum over equivalence classes
  uvec uclasses = unique(equiv_classes);
  MATTYPE R_new = zeros<MATTYPE>(uclasses.n_elem, R.n_cols); 
  for (unsigned i = 0; i < R_new.n_rows; i++) {
    uvec idx = find(equiv_classes == uclasses(i)); 
    R_new.row(i) = sum(R.rows(idx), 0);
  }
  return R_new;  
}

// // [[Rcpp::export]]
// MATTYPE compute_Y(const MATTYPE& Z_cos, const MATTYPE& R, const VECTYPE& weights) {
// //     return arma::normalise(Z_cos * R.t(), 2, 0);
//     MATTYPE Rw = R * arma::diagmat(weights);
//     Rw = arma::diagmat(1 / arma::sum(Rw, 1)) * (Rw * Z_cos.t());
//     return arma::normalise(Rw.t(), 2, 0);
// }


// [[Rcpp::export]]
MATTYPE scaleRows_dgc(const VECTYPE& x, const VECTYPE& p, const VECTYPE& i, 
                        int ncol, int nrow, float thresh) {
    // (0) fill in non-zero elements
    MATTYPE res = arma::zeros<MATTYPE>(nrow, ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            res(i[j], c) = x(j);
        }
    }

    // (1) compute means
    VECTYPE mean_vec = arma::zeros<VECTYPE>(nrow);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            mean_vec(i[j]) += x[j];
        }
    }
    mean_vec /= ncol;
    
    // (2) compute SDs
    VECTYPE sd_vec = arma::zeros<VECTYPE>(nrow);
    arma::uvec nz = arma::zeros<arma::uvec>(nrow);
    nz.fill(ncol);
    for (int c = 0; c < ncol; c++) {
        for (int j = p[c]; j < p[c + 1]; j++) {
            sd_vec(i[j]) += (x[j] - mean_vec(i[j])) * (x[j] - mean_vec(i[j])); // (x - mu)^2
            nz(i[j])--;
        }
    }
        
    // count for the zeros
    for (int r = 0; r < nrow; r++) {
        sd_vec(r) += nz(r) * mean_vec(r) * mean_vec(r);
    }
    
    sd_vec = arma::sqrt(sd_vec / (ncol - 1));
    
    // (3) scale values
    res.each_col() -= mean_vec;
    res.each_col() /= sd_vec;
    res.elem(find(res > thresh)).fill(thresh);
    res.elem(find(res < -thresh)).fill(-thresh);
    return res;
}

// [[Rcpp::export]]
double soft_kmeans_score_cpp(const arma::mat& R, const arma::rowvec& w, const arma::mat& dist_mat, float sigma) {
    float score_dist = arma::as_scalar(arma::accu((R * arma::diagmat(w)) % dist_mat)); 
//     float score_dist = arma::as_scalar(arma::accu((R.each_row() % w) % dist_mat)); 
    float score_entropy = arma::as_scalar(arma::accu((safe_entropy(R)) * arma::diagmat(w)));
    return score_dist + sigma*score_entropy;
}

bool kmeans_converged(const vector<float>& scores, float tol) {
    float s0 = scores[scores.size() - 2];
    float s1 = scores[scores.size() - 1];
    return (s0 - s1) / s0 < tol;
}

// [[Rcpp::export]]
List soft_kmeans_weighted_cpp(arma::mat Y, arma::mat Z, const arma::rowvec& w, unsigned max_iter, float sigma, float tol) {
    Y = arma::normalise(Y, 2, 0); // L2 normalize the columns
    Z = arma::normalise(Z, 2, 0); // L2 normalize the columns
    arma::mat R, Rw, dist_mat; 
    std::vector<float> scores;
    float s0, s1;
    for (unsigned i = 0; i < max_iter; i++) {
        dist_mat = 2 * (1 - Y.t() * Z); 
        R = -dist_mat / sigma;
        R.each_row() -= arma::max(R, 0);  
        R = exp(R);
        R.each_row() /= arma::sum(R, 0);
        Rw = R * arma::diagmat(w); 
        scores.push_back(soft_kmeans_score_cpp(R, w, dist_mat, sigma));
        Y = arma::normalise(Z * Rw.t(), 2, 0); 
        if (i > 0 && kmeans_converged(scores, tol)) break;
    }

    List result = List::create(_("R") = R , _["Y"] = Y, _["scores"] = scores);
    return result;    
}




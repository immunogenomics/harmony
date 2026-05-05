#include "utils.h"
#include "types.h"
#include <Rcpp.h>
#include <progress.hpp>
#include <time.h>
#include <cstdlib> 



MATTYPE initialize_centroids(const MATTYPE& X, const unsigned int K, bool verbose) {
  // K-means++ centroid initialization
  VECTYPE random_seeds(K, arma::fill::randu);
  int N = (X.n_cols-0.000001); // subtract a small number so randu=1 won't give out of bounds
  arma::uvec indices = arma::conv_to<arma::uvec>::from(floor(random_seeds * N));
  
  MATTYPE Y(X.cols(indices));    
  if (verbose) {
    Rcpp::Rcout << "Initializing centroids" << std::endl;
  }
  
  Progress p(K, verbose);
  std::set<unsigned> sup;
  
  // k-means++
  for (unsigned int i = 0; i < K; i++) {
    p.increment();
    
    VECTYPE distances = arma::abs((2* (1 - Y.col(i).t() * X)).as_col());    
    VECTYPE random_numbers(size(distances), arma::fill::randu);
    
    // Weighted Random Sampling, sample from different expontential
    // distributions with distance as different rate parameters
    VECTYPE prob = -arma::log(random_numbers) / distances;
    
    auto index = prob.index_min();
    
    // Make sure we have not selected the same point for cluster centroid already
    // This can be happen particularly in small datasets
    while (sup.find(index) != sup.end()) {
      Rcpp::Rcerr << index << "exists, retrying for cluster " << i << " " << distances(index) <<std::endl;
      prob[index] = prob.max();
      index = prob.index_min();      
    }
    sup.insert(index);
    Y.col(i) = X.col(index);
  }
  
  return Y;
}


//[[Rcpp::export]]
MATTYPE kmeans_centers(const MATTYPE& X, const unsigned int K, bool verbose) {

  MATTYPE Y = initialize_centroids(X, K, verbose);
  unsigned iterations = 10;
  for(unsigned i = 0; i < iterations; i++) {
    if (!arma::kmeans(Y, X, K, arma::keep_existing, 1, DEBUG)) {
      Rcpp::stop("Clustering failed");
    }
  }
    
  return Y;
}


float my_accu(const MATTYPE& X) {
  auto* X_mem = X.memptr();
  float sum=0;
  long len =X.n_rows * X.n_cols;
  for(long i = 0; i < len; ++i){
    sum+=X_mem[i];
  }
  return sum;
}

MATTYPE safe_entropy(const MATTYPE& X) {
  return X % trunc_log(X);
  // A.elem(find_nonfinite(A)).zeros();
  // return(A);
}

// Overload pow to work on a MATTYPErix and vector
MATTYPE harmony_pow(MATTYPE A, const VECTYPE& T) {

  for (unsigned c = 0; c < A.n_cols; c++) {
    A.unsafe_col(c) = pow(A.unsafe_col(c), as_scalar(T.row(c)));
  }
  return(A);
}

VECTYPE calculate_norm(const MATTYPE& M) {
  VECTYPE x(M.n_cols);
  for(unsigned i = 0; i < M.n_cols; i++){
    x(i) = norm(M.col(i));
  }
  return x;
}


//https://stackoverflow.com/questions/8377412/ceil-function-how-can-we-implement-it-ourselves
int my_ceil(float num) {
    int inum = (int)num;
    if (num == (float)inum) {
        return inum;
    }
    return inum + 1;
}


// [[Rcpp::export]]
MATTYPE scaleRows_dgc(const VECTYPE& x, const VECTYPE& p, const VECTYPE& i, int ncol, int nrow, float thresh) {
  
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
VECTYPE find_lambda_cpp(const float alpha, const VECTYPE& cluster_E) {
  VECTYPE lambda_dym_vec(cluster_E.n_rows + 1, arma::fill::zeros);
  lambda_dym_vec.subvec(1, lambda_dym_vec.n_rows - 1) = cluster_E * alpha;
  return lambda_dym_vec;
}





std::vector< std::pair<unsigned,unsigned> > find_contigs(std::vector<unsigned>& keep_vectors) {
  unsigned kprev = keep_vectors[0], k;
  std::vector< std::pair<unsigned,unsigned> > ranges;
  unsigned i0 = 0, i;
  ranges.reserve(10000);
  for (i = 1; i < keep_vectors.size(); ++i) {
    k = keep_vectors[i];
    if (k - kprev  != 1) {
      // std::cout << "i0 = " << i0 << " i= "  << i << std::endl;
      ranges.push_back({i0, i - 1});
      
      i0 = i;
    }
    kprev = k;
  }
  ranges.push_back({i0, i-1});
  return ranges;
}

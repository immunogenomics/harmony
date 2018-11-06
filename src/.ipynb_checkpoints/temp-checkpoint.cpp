#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
mat moe_correct_armadillo(const mat& Z, const mat& phi, const mat& R) {
//    mat Z_corr = zeros<mat>(Z.n_rows, Z.n_cols); 
    int K = R.n_rows;
    mat W = zeros<mat>(K, phi.n_rows);
    
    for (int k = 0; k < K; k++) { 
        W.col(k) = inv(phi * diagmat(R.row(k)) * phi.t()) * phi * diagmat(R.row(k)) * Z.t();
    }
    
//    Rcout << "HELLO 1" << endl;
//    Rcout << "HELLO 2" << endl;
    
    
    return(Z - sum(R % (W.rows(1, W.n_rows - 1).t() * phi.rows(1, phi.n_rows - 1)), 0));   
//    return(Z - (R % (cross(W.rows(1, W.n_rows - 1), phi.rows(1, phi.n_rows - 1)))));    
}
  




//    Z_corr <- Z - colSums(R * crossprod(W[2:nrow(W), , drop = FALSE], phi[2:nrow(phi), , drop = FALSE]))
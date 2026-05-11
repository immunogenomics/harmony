#pragma once
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

#ifdef HARMONY_SCALAR_DOUBLE
typedef double SCALAR;
#else
typedef float SCALAR;
#endif

typedef arma::Mat<SCALAR> MATTYPE;
typedef arma::mat RMAT;

typedef arma::SpMat<SCALAR> SPMAT;
typedef arma::sp_mat RSPMAT;

typedef arma::Col<SCALAR> VECTYPE;
typedef arma::vec RVEC;


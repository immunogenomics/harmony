#include <RcppArmadillo.h>

#if USE_FLOAT_TYPES==1
typedef arma::fmat MATTYPE;
typedef arma::fvec VECTYPE;
typedef arma::frowvec ROWVECTYPE;
typedef arma::fcube CUBETYPE;
#endif

#if USE_FLOAT_TYPES==0
typedef arma::mat MATTYPE;
typedef arma::vec VECTYPE;
typedef arma::rowvec ROWVECTYPE;
typedef arma::cube CUBETYPE;
#endif


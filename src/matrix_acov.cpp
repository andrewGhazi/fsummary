// #include <Rcpp.h>
#include <RcppArmadillo/Lightest>
using namespace Rcpp;
#define ARMA_USE_FFTW3

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat fftm(NumericMatrix x) {

  int nr = x.nrow();
  int nc = x.ncol();
  int ne = nr*nc;
  
  arma::cx_mat Y (2*nr, nc);
  
  Y = arma::fft(as<arma::mat>(x), 2*nr);

  for (int i=0; i < ne; i++) {
    // abs(Y)^2 involves a redundant sqrt then square
    Y[i] = pow(real(Y[i]), 2) + pow(imag(Y[i]), 2);
  }
  
  Y.set_imag(arma::zeros(2*nr,nc));
  
  return (real(ifft(Y)));
}

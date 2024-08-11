// #include <Rcpp.h>
#include <RcppArmadillo/Lightest>
using namespace Rcpp;
// #define ARMA_USE_FFTW3
// on my machine, FFTW3 is a tiny bit slower, plus it's an additional dependency.

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat fftm(NumericMatrix x) {

  int nr = x.nrow();
  int nc = x.ncol();

  arma::cx_mat Y (2*nr, nc);

  Y = arma::fft(as<arma::mat>(x), 2*nr);

  Y.set_real(square(real(Y)) + square(imag(Y)));
  Y.set_imag(arma::zeros(2*nr, nc));

  return (real(ifft(Y)));
}

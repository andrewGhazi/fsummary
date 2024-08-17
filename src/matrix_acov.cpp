// #include <Rcpp.h>
#include <RcppArmadillo/Lightest>
using namespace Rcpp;

// #define ARMA_USE_FFTW3
// on my machine, FFTW3 is a tiny bit slower, plus it's an additional dependency.

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat fftm(arma::mat x, int k) {

  // int nr = x.n_rows;
  int nc = x.n_cols;

  arma::cx_mat Y (2*k, nc, arma::fill::none);

  Y = arma::fft(x, 2*k);

  arma::cx_mat::iterator it     = Y.begin();
  arma::cx_mat::iterator it_end = Y.end();

  for(; it != it_end; ++it)
  {
    (*it).real((pow(real(*it), 2) + pow(imag(*it), 2)));
    (*it).imag(0);
  }

  // arma::mat res (k, nc, arma::fill::none);

  x = real(ifft(Y)).eval().head_rows(k);

  return (x);
}

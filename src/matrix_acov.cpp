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

// This is actually about 32ms slower for a 12288x2652 ddf
// // [[Rcpp::export]]
// arma::mat fftm2(arma::mat x) {
//
//   int nr = x.n_rows;
//   int nc = x.n_cols;
//
//   arma::cx_mat Y (2*nr, nc);
//
//   Y = arma::fft(x, 2*nr);
//
//   Y.set_real(square(real(Y)) + square(imag(Y)));
//   Y.set_imag(arma::zeros(2*nr, nc));
//
//   return (real(ifft(Y)));
// }

// [[Rcpp::export]]
arma::mat fftm3(arma::mat x) {

  int nr = x.n_rows;
  int nc = x.n_cols;

  arma::cx_mat Y (2*nr, nc, arma::fill::none);

  Y = arma::fft(x, 2*nr);

  Y.set_real(square(real(Y)) + square(imag(Y)));
  Y.set_imag(arma::zeros(2*nr, nc));

  return (real(ifft(Y)));
}

// [[Rcpp::export]]
arma::mat fftm4(arma::mat x) {

  int nr = x.n_rows;
  int nc = x.n_cols;

  arma::cx_mat Y (2*nr, nc, arma::fill::none);

  Y = arma::fft(x, 2*nr);

  arma::cx_mat::iterator it     = Y.begin();
  arma::cx_mat::iterator it_end = Y.end();

  for(; it != it_end; ++it)
  {
    (*it).real(pow(real(*it), 2) + pow(imag(*it), 2));
    (*it).imag(0);
  }

  return (real(ifft(Y)));
}

// [[Rcpp::export]]
arma::mat fftm5(NumericMatrix x) {

  int nr = x.nrow();
  int nc = x.ncol();

  arma::cx_mat Y(2*nr, nc, arma::fill::none);

  Y = arma::fft(as<arma::mat>(x), 2*nr);

  arma::cx_mat::iterator it     = Y.begin();
  arma::cx_mat::iterator it_end = Y.end();

  for(; it != it_end; ++it)
  {
    (*it).real(pow(real(*it), 2) + pow(imag(*it), 2));
    (*it).imag(0);
  }

  return (real(ifft(Y)));
}

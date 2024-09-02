// #include <Rcpp.h>
#include <RcppArmadillo/Lightest>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::export]]
// arma::mat fftm(arma::mat x, int k) {
//
//   // int nr = x.n_rows;
//   int nc = x.n_cols;
//
//   arma::cx_mat Y (2*k, nc, arma::fill::none);
//
//   Y = arma::fft(x, 2*k);
//
//   arma::cx_mat::iterator it     = Y.begin();
//   arma::cx_mat::iterator it_end = Y.end();
//
//   for(; it != it_end; ++it)
//   {
//     (*it).real((pow(real(*it), 2) + pow(imag(*it), 2)));
//     (*it).imag(0);
//   }
//
//   // arma::mat res (k, nc, arma::fill::none);
//
//   x = real(ifft(Y)).eval().head_rows(k);
//
//   return (x);
// }

// [[Rcpp::export]]
arma::mat cov_head(arma::mat x, int n, int offset) {

  int nc = x.n_cols;
  int nr = x.n_rows;

  arma::mat res (n, nc);

  for (int j=0; j<nc; j++) {
    for (int t=0; t < n; t++) {
      int t_off = t+offset;
      res(t,j) += sum(x.col(j).head(nr-t_off) % x.col(j).tail(nr-t_off)) ;
    }
  }

  // res *= (nr-1)/nr;

  return (res);
}

// [[Rcpp::export]]
arma::vec myrank(arma::vec v) {

  int n = v.n_elem;
  uvec o = sort_index(v);
  arma::vec res = conv_to< arma::vec >::from(sort_index(o) + 1);
  arma::vec sv = v(o);

  int uniq_obs = 1;

  for (int i=1; i<n; i++) {
    // scan through the sorted vector
    int need_div = 0;

    if (sv(i) == sv(i-1)) {
      uniq_obs += 1;

      if ((i < (n-1) && sv(i) != sv(i+1)) || ((i+2) > n)) {
        need_div = 1;
      }
    }

    if (need_div) {
      // set the previous uniq_obs elements of res to their mean
      double mv = 0;

      for (int j=0; j<uniq_obs; j++) {
        mv += res(o(i-j));
      }

      mv /= uniq_obs;

      for (int j=0; j<uniq_obs; j++) {
        res(o(i-j)) = mv;
      }

      uniq_obs = 1;
      need_div = 0;
    }

  }

  return (res.as_col());
}

// [[Rcpp::export]]
List center_split_df(DataFrame df, IntegerVector c_id, int n_chain, int n_iter) {
  List L(n_chain);

  int nc = df.ncol();

  vec sum_vec(n_chain);

  for (int k=0; k<n_chain; k++) {
    int imin = k*n_iter;
    int imax = (k+1) * n_iter;

    mat Lk = arma::zeros(n_iter, nc);

    for (int j=0; j<nc; j++) {
      NumericVector v = df[j];
      double param_mean = 0;

      for (int i=imin; i<imax; i++) {
        param_mean += v(i);
      }

      param_mean /= n_iter;

      for (int i=0; i<n_iter; i++) {
        Lk(i,j) = v(i+k*n_iter) - param_mean;
      }
    }

    L[k] = Lk;
  }

  return (L);
}

// // [[Rcpp::export]]
// DataFrame dfranks(DataFrame df) {
//
//   int nr = df.nrow();
//   int nc = df.ncol();
//   Function f2("qnorm");
//
//   for (int i=0; i<(nc-3); i++) {
//     // skip the last three columns - these are chain information
//
//     // (myrank(df[i]) - 3.0/8.0) / (nr - (2*3.0/8.0) + 1)
//     df[i] = as<arma::vec>(f2((- 3.0/8.0 + myrank(df[i]) ) / (nr - (2*3.0/8.0) + 1)));
//
//   }
//
//   return(df);
// }
// ^ This doesn't help at all

// // [[Rcpp::export]]
// DataFrame df_zscale(DataFrame df) {
//
//   int nc = df.ncol();
//   int nr = df.ncol();
//   // Function rank("rank");
//   // Function qnorm("qnorm");
//
//   for (int i=0; i<(nc-3); i++) {
//
//     // z_scaled = ddff |> # 3.7s
//     // mtt(across(variables,
//     //            \(x) qnorm((rank_fun(x) - 3/8) / (nrow(ddff) - 2 * 3/8 + 1))))
//     // skip the last three columns - these are chain information
//     NumericVector v = df[i];
//     NumericVector rankv = wrap(myrank(v));
//     v = qnorm((-3.0/8.0 + rankv ) / (-1*(2.0 * 3.0/8.0) + 1.0 + nr));
//   }
//
//   return(df);
// }

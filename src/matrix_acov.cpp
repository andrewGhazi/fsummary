// #include <Rcpp.h>
#include <RcppArmadillo/Lightest>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

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

// [[Rcpp::export]]
arma::vec fqnorm3v(arma::vec p) {
  double A0 = 3.3871328727963666080e0;
  double A1 = 1.3314166789178437745e2;
  double A2 = 1.9715909503065514427e3;
  double A3 = 1.3731693765509461125e4;
  double A4 = 4.5921953931549871457e4;
  double A5 = 6.7265770927008700853e4;
  double A6 = 3.3430575583588128105e4;
  double A7 = 2.5090809287301226727e3;
  double B1 = 4.2313330701600911252e1;
  double B2 = 6.8718700749205790830e2;
  double B3 = 5.3941960214247511077e3;
  double B4 = 2.1213794301586595867e4;
  double B5 = 3.9307895800092710610e4;
  double B6 = 2.8729085735721942674e4;
  double B7 = 5.2264952788528545610e3;

  double C0 = 1.42343711074968357734e0;
  double C1 = 4.63033784615654529590e0;
  double C2 = 5.76949722146069140550e0;
  double C3 = 3.64784832476320460504e0;
  double C4 = 1.27045825245236838258e0;
  double C5 = 2.41780725177450611770e-1;
  double C6 = 2.27238449892691845833e-2;
  double C7 = 7.74545014278341407640e-4;
  double D1 = 2.05319162663775882187e0;
  double D2 = 1.67638483018380384940e0;
  double D3 = 6.89767334985100004550e-1;
  double D4 = 1.48103976427480074590e-1;
  double D5 = 1.51986665636164571966e-2;
  double D6 = 5.47593808499534494600e-4;
  double D7 = 1.05075007164441684324e-9;

  double E0 = 6.65790464350110377720e0;
  double E1 = 5.46378491116411436990e0;
  double E2 = 1.78482653991729133580e0;
  double E3 = 2.96560571828504891230e-1;
  double E4 = 2.65321895265761230930e-2;
  double E5 = 1.24266094738807843860e-3;
  double E6 = 2.71155556874348757815e-5;
  double E7 = 2.01033439929228813265e-7;
  double F1 = 5.99832206555887937690e-1;
  double F2 = 1.36929880922735805310e-1;
  double F3 = 1.48753612908506148525e-2;
  double F4 = 7.86869131145613259100e-4;
  double F5 = 1.84631831751005468180e-5;
  double F6 = 1.42151175831644588870e-7;
  double F7 = 2.04426310338993978564e-15;

  arma::vec q = p - .5;
  arma::vec res(p.n_elem, fill::none);
  arma::vec r(p.n_elem, fill::none);

  double ri = 1.0;
  double qi = 1.0;

  for (arma::uword i=0; i<p.n_elem; i++) {

    ri = r[i];
    qi = q[i];

    // Rcout << "The value of qi : " << (fabs(qi) < .425) << "\n"; arma::abs would return true here wrongly for some reason
    if (fabs(qi) < .425) {

      ri = .180625 - qi*qi;

      res[i] = qi * (((((((ri * A7 + A6) * ri + A5) * ri + A4) * ri + A3) * ri + A2) * ri + A1) * ri + A0) /
        (((((((ri * B7 + B6) * ri + B5) * ri + B4) * ri + B3) * ri + B2) * ri + B1) * ri + 1.0);
    } else {
      if (qi < 0) {
        ri = p[i];
      } else {
        ri = 1.0 - p[i];
      }

      ri = sqrt(-log(ri));

      if (ri <= 5.0) {
        ri = ri - 1.6;

        res[i] = (((((((ri * C7 + C6) * ri + C5) * ri + C4) * ri + C3) * ri + C2) * ri + C1) * ri + C0) /
          (((((((ri * D7 + D6) * ri + D5) * ri + D4) * ri + D3) * ri + D2) * ri + D1) * ri + 1.0);
      } else {
        ri = ri - 5.0;

        res[i] = (((((((ri * E7 + E6) * ri + E5) * ri + E4) * ri + E3) * ri + E2) * ri + E1) * ri + E0) /
          (((((((ri * F7 + F6) * ri + F5) * ri + F4) * ri + F3) * ri + F2) * ri + F1) * ri + 1.0);
      }

      if (q[i] < 0) res[i] = -res[i];
    }
  }

  return (res);
}

// [[Rcpp::export]]
arma::vec fqnorm(arma::vec p) {
  // Adapted from:
  // Wichura, M. J. (1988) Algorithm AS 241: The percentage points of the normal distribution. Applied Statistics, 37, 477–484; doi:10.2307/2347330.

  double A0 = 3.3871328727963666080e0;
  double A1 = 1.3314166789178437745e2;
  double A2 = 1.9715909503065514427e3;
  double A3 = 1.3731693765509461125e4;
  double A4 = 4.5921953931549871457e4;
  double A5 = 6.7265770927008700853e4;
  double A6 = 3.3430575583588128105e4;
  double A7 = 2.5090809287301226727e3;
  double B1 = 4.2313330701600911252e1;
  double B2 = 6.8718700749205790830e2;
  double B3 = 5.3941960214247511077e3;
  double B4 = 2.1213794301586595867e4;
  double B5 = 3.9307895800092710610e4;
  double B6 = 2.8729085735721942674e4;
  double B7 = 5.2264952788528545610e3;

  double C0 = 1.42343711074968357734e0;
  double C1 = 4.63033784615654529590e0;
  double C2 = 5.76949722146069140550e0;
  double C3 = 3.64784832476320460504e0;
  double C4 = 1.27045825245236838258e0;
  double C5 = 2.41780725177450611770e-1;
  double C6 = 2.27238449892691845833e-2;
  double C7 = 7.74545014278341407640e-4;
  double D1 = 2.05319162663775882187e0;
  double D2 = 1.67638483018380384940e0;
  double D3 = 6.89767334985100004550e-1;
  double D4 = 1.48103976427480074590e-1;
  double D5 = 1.51986665636164571966e-2;
  double D6 = 5.47593808499534494600e-4;
  double D7 = 1.05075007164441684324e-9;

  double E0 = 6.65790464350110377720e0;
  double E1 = 5.46378491116411436990e0;
  double E2 = 1.78482653991729133580e0;
  double E3 = 2.96560571828504891230e-1;
  double E4 = 2.65321895265761230930e-2;
  double E5 = 1.24266094738807843860e-3;
  double E6 = 2.71155556874348757815e-5;
  double E7 = 2.01033439929228813265e-7;
  double F1 = 5.99832206555887937690e-1;
  double F2 = 1.36929880922735805310e-1;
  double F3 = 1.48753612908506148525e-2;
  double F4 = 7.86869131145613259100e-4;
  double F5 = 1.84631831751005468180e-5;
  double F6 = 1.42151175831644588870e-7;
  double F7 = 2.04426310338993978564e-15;

  double ri = 1.0;
  double qi = 1.0;

  arma::mat::iterator it = p.begin();
  arma::mat::iterator it_end = p.end();

  for (; it != it_end; ++it) {

    qi = (*it) - 0.5;

    // Rcout << "The value of qi : " << (fabs(qi) < .425) << "\n"; arma::abs
    // would return true here wrongly for some reason. Guess it's only for vecs,
    // not doubles?
    if (fabs(qi) < .425) {

      ri = .180625 - qi*qi;

      (*it) = qi * (((((((ri * A7 + A6) * ri + A5) * ri + A4) * ri + A3) * ri + A2) * ri + A1) * ri + A0) /
        (((((((ri * B7 + B6) * ri + B5) * ri + B4) * ri + B3) * ri + B2) * ri + B1) * ri + 1.0);
    } else {

      if (qi < 0) {
        ri = qi+0.5;
      } else {
        ri = .5 - qi;
      }

      ri = sqrt(-log(ri));

      if (ri <= 5.0) {
        ri = ri - 1.6;

        (*it) = (((((((ri * C7 + C6) * ri + C5) * ri + C4) * ri + C3) * ri + C2) * ri + C1) * ri + C0) /
          (((((((ri * D7 + D6) * ri + D5) * ri + D4) * ri + D3) * ri + D2) * ri + D1) * ri + 1.0);
      } else {
        ri = ri - 5.0;

        (*it) = (((((((ri * E7 + E6) * ri + E5) * ri + E4) * ri + E3) * ri + E2) * ri + E1) * ri + E0) /
          (((((((ri * F7 + F6) * ri + F5) * ri + F4) * ri + F3) * ri + F2) * ri + F1) * ri + 1.0);
      }

      if (qi < 0) (*it) = -(*it);
    }
  }

  return(p);
}

// [[Rcpp::export]]
arma::vec fqnorm6v(arma::vec p) {
  double A0 = 3.3871328727963666080e0;
  double A1 = 1.3314166789178437745e2;
  double A2 = 1.9715909503065514427e3;
  double A3 = 1.3731693765509461125e4;
  double A4 = 4.5921953931549871457e4;
  double A5 = 6.7265770927008700853e4;
  double A6 = 3.3430575583588128105e4;
  double A7 = 2.5090809287301226727e3;
  double B1 = 4.2313330701600911252e1;
  double B2 = 6.8718700749205790830e2;
  double B3 = 5.3941960214247511077e3;
  double B4 = 2.1213794301586595867e4;
  double B5 = 3.9307895800092710610e4;
  double B6 = 2.8729085735721942674e4;
  double B7 = 5.2264952788528545610e3;

  double C0 = 1.42343711074968357734e0;
  double C1 = 4.63033784615654529590e0;
  double C2 = 5.76949722146069140550e0;
  double C3 = 3.64784832476320460504e0;
  double C4 = 1.27045825245236838258e0;
  double C5 = 2.41780725177450611770e-1;
  double C6 = 2.27238449892691845833e-2;
  double C7 = 7.74545014278341407640e-4;
  double D1 = 2.05319162663775882187e0;
  double D2 = 1.67638483018380384940e0;
  double D3 = 6.89767334985100004550e-1;
  double D4 = 1.48103976427480074590e-1;
  double D5 = 1.51986665636164571966e-2;
  double D6 = 5.47593808499534494600e-4;
  double D7 = 1.05075007164441684324e-9;

  double E0 = 6.65790464350110377720e0;
  double E1 = 5.46378491116411436990e0;
  double E2 = 1.78482653991729133580e0;
  double E3 = 2.96560571828504891230e-1;
  double E4 = 2.65321895265761230930e-2;
  double E5 = 1.24266094738807843860e-3;
  double E6 = 2.71155556874348757815e-5;
  double E7 = 2.01033439929228813265e-7;
  double F1 = 5.99832206555887937690e-1;
  double F2 = 1.36929880922735805310e-1;
  double F3 = 1.48753612908506148525e-2;
  double F4 = 7.86869131145613259100e-4;
  double F5 = 1.84631831751005468180e-5;
  double F6 = 1.42151175831644588870e-7;
  double F7 = 2.04426310338993978564e-15;


  arma::vec res = p;

  res.transform( [&A7, &A6, &A5, &A4, &A3, &A2, &A1, &A0, &B7, &B6, &B5, &B4, &B3, &B2, &B1, &C7, &C6, &C5, &C4, &C3, &C2, &C1, &C0, &D7, &D6, &D5, &D4, &D3, &D2, &D1, &E7, &E6, &E5, &E4, &E3, &E2, &E1, &E0, &F7, &F6, &F5, &F4, &F3, &F2, &F1](double val) {
    double qi = val - 0.5;
    double ri = 1.0;
    double resi = 0.0;

    // Rcout << "The value of qi : " << (fabs(qi) < .425) << "\n"; arma::abs would return true here wrongly for some reason
    if (fabs(qi) < .425) {

      ri = .180625 - qi*qi;

      resi = qi * (((((((ri * A7 + A6) * ri + A5) * ri + A4) * ri + A3) * ri + A2) * ri + A1) * ri + A0) /
        (((((((ri * B7 + B6) * ri + B5) * ri + B4) * ri + B3) * ri + B2) * ri + B1) * ri + 1.0);
    } else {
      if (qi < 0) {
        ri = val;
      } else {
        ri = 1.0 - val;
      }

      ri = sqrt(-log(ri));

      if (ri <= 5.0) {
        ri = ri - 1.6;

        resi = (((((((ri * C7 + C6) * ri + C5) * ri + C4) * ri + C3) * ri + C2) * ri + C1) * ri + C0) /
          (((((((ri * D7 + D6) * ri + D5) * ri + D4) * ri + D3) * ri + D2) * ri + D1) * ri + 1.0);
      } else {
        ri = ri - 5.0;

        resi = (((((((ri * E7 + E6) * ri + E5) * ri + E4) * ri + E3) * ri + E2) * ri + E1) * ri + E0) /
          (((((((ri * F7 + F6) * ri + F5) * ri + F4) * ri + F3) * ri + F2) * ri + F1) * ri + 1.0);
      }

      if (qi < 0) resi = -resi;
    }

    return(resi);
  });

  return(res);
}

// [[Rcpp::export]]
arma::vec fqnorm7v(arma::vec p) {
  arma::vec A = {3.3871328727963666080e0,
  1.3314166789178437745e2,
  1.9715909503065514427e3,
  1.3731693765509461125e4,
  4.5921953931549871457e4,
  6.7265770927008700853e4,
  3.3430575583588128105e4,
  2.5090809287301226727e3};
  arma::vec B = {4.2313330701600911252e1,
   6.8718700749205790830e2,
   5.3941960214247511077e3,
   2.1213794301586595867e4,
   3.9307895800092710610e4,
   2.8729085735721942674e4,
   5.2264952788528545610e3};

  arma::vec C = {1.42343711074968357734e0,
   4.63033784615654529590e0,
   5.76949722146069140550e0,
   3.64784832476320460504e0,
   1.27045825245236838258e0,
   2.41780725177450611770e-1,
   2.27238449892691845833e-2,
   7.74545014278341407640e-4};
  arma::vec D = {2.05319162663775882187e0,
   1.67638483018380384940e0,
   6.89767334985100004550e-1,
   1.48103976427480074590e-1,
   1.51986665636164571966e-2,
   5.47593808499534494600e-4,
   1.05075007164441684324e-9};

  arma::vec E = {6.65790464350110377720e0,
   5.46378491116411436990e0,
   1.78482653991729133580e0,
   2.96560571828504891230e-1,
   2.65321895265761230930e-2,
   1.24266094738807843860e-3,
   2.71155556874348757815e-5,
   2.01033439929228813265e-7};
  arma::vec F = {5.99832206555887937690e-1,
   1.36929880922735805310e-1,
   1.48753612908506148525e-2,
   7.86869131145613259100e-4,
   1.84631831751005468180e-5,
   1.42151175831644588870e-7,
   2.04426310338993978564e-15};

  double ri = 1.0;
  double qi = 1.0;
  double n = 1.0;
  double d = 1.0;
  double rp = 1.0;

  arma::mat::iterator it = p.begin();
  arma::mat::iterator it_end = p.end();

  for (; it != it_end; ++it){

    qi = (*it) - 0.5;

    // Rcout << "The value of qi : " << (fabs(qi) < .425) << "\n"; arma::abs would return true here wrongly for some reason
    if (fabs(qi) < .425) {

      ri = .180625 - qi*qi;
      n = A[0];
      d = 1.0;
      rp = 1.0;

      for (arma::uword j=1; j<8; j++) {
        rp *= ri;
        n += A[j]*rp;
        d += B[j-1]*rp;
      }

      (*it) = qi * n / d;
    } else {

      if (qi < 0) {
        ri = qi+0.5;
      } else {
        ri = .5 - qi;
      }

      ri = sqrt(-log(ri));

      if (ri <= 5.0) {
        ri = ri - 1.6;

        n = C[0];
        d = 1.0;
        rp = 1.0;

        for (arma::uword j=1; j<8; j++) {
          rp *= ri;
          n += C[j]*rp;
          d += D[j-1]*rp;
        }

        (*it) = n / d;

      } else {
        ri = ri - 5.0;

        n = E[0];
        d = 1.0;
        rp = 1.0;

        for (arma::uword j=1; j<8; j++) {
          rp *= ri;
          n += E[j]*rp;
          d += F[j-1]*rp;
        }

        (*it) = n / d;
      }

      if (qi < 0) (*it) = -(*it);
    }
  }

  return(p);
}

// // [[Rcpp::export]]
// arma::vec fqnorm5v(arma::vec p) {
//   arma::vec A = {3.3871328727963666080e0,
//   1.3314166789178437745e2,
//   1.9715909503065514427e3,
//   1.3731693765509461125e4,
//   4.5921953931549871457e4,
//   6.7265770927008700853e4,
//   3.3430575583588128105e4,
//   2.5090809287301226727e3};
//   arma::vec B = {4.2313330701600911252e1,
//    6.8718700749205790830e2,
//    5.3941960214247511077e3,
//    2.1213794301586595867e4,
//    3.9307895800092710610e4,
//    2.8729085735721942674e4,
//    5.2264952788528545610e3};
//
//   arma::vec C = {1.42343711074968357734e0,
//    4.63033784615654529590e0,
//    5.76949722146069140550e0,
//    3.64784832476320460504e0,
//    1.27045825245236838258e0,
//    2.41780725177450611770e-1,
//    2.27238449892691845833e-2,
//    7.74545014278341407640e-4};
//   arma::vec D = {2.05319162663775882187e0,
//    1.67638483018380384940e0,
//    6.89767334985100004550e-1,
//    1.48103976427480074590e-1,
//    1.51986665636164571966e-2,
//    5.47593808499534494600e-4,
//    1.05075007164441684324e-9};
//
//   arma::vec E = {6.65790464350110377720e0,
//    5.46378491116411436990e0,
//    1.78482653991729133580e0,
//    2.96560571828504891230e-1,
//    2.65321895265761230930e-2,
//    1.24266094738807843860e-3,
//    2.71155556874348757815e-5,
//    2.01033439929228813265e-7};
//   arma::vec F = {5.99832206555887937690e-1,
//    1.36929880922735805310e-1,
//    1.48753612908506148525e-2,
//    7.86869131145613259100e-4,
//    1.84631831751005468180e-5,
//    1.42151175831644588870e-7,
//    2.04426310338993978564e-15};
//
//   double ri = 1.0;
//   double qi = 1.0;
//
//   arma::mat::iterator it = p.begin();
//   arma::mat::iterator it_end = p.end();
//
//   for (; it != it_end; ++it){
//
//     qi = (*it) - 0.5;
//
//     // Rcout << "The value of qi : " << (fabs(qi) < .425) << "\n"; arma::abs would return true here wrongly for some reason
//     if (fabs(qi) < .425) {
//
//       ri = .180625 - qi*qi;
//       arma::vec riv = {ri, ri, ri, ri, ri, ri, ri, ri};
//       riv = pow(riv, linspace(0.0, 7.0, 8));
//
//       (*it) = qi * sum(A % riv) / (sum(B % riv.tail(7)) + 1.0) ;
//     } else {
//
//       if (qi < 0) {
//         ri = qi+0.5;
//       } else {
//         ri = .5 - qi;
//       }
//
//       ri = sqrt(-log(ri));
//
//       if (ri <= 5.0) {
//         ri = ri - 1.6;
//
//         arma::vec riv = {ri, ri, ri, ri, ri, ri, ri, ri};
//         riv = pow(riv, linspace(0.0, 7.0, 8));
//         (*it) = sum(C % riv) / (sum(D % riv.tail(7) + 1.0));
//       } else {
//         ri = ri - 5.0;
//
//         arma::vec riv = {ri, ri, ri, ri, ri, ri, ri, ri};
//         riv = pow(riv, linspace(0.0, 7.0, 8));
//         (*it) = sum(E % riv) / (sum(F % riv.tail(7) + 1.0)) ;
//       }
//
//       if (qi < 0) (*it) = -(*it);
//     }
//   }
//
//   return(p);
// }
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


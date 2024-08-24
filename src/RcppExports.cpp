// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cov_head
arma::mat cov_head(arma::mat x, int n, int offset);
RcppExport SEXP _fsummary_cov_head(SEXP xSEXP, SEXP nSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(cov_head(x, n, offset));
    return rcpp_result_gen;
END_RCPP
}
// myrank
arma::vec myrank(arma::vec v);
RcppExport SEXP _fsummary_myrank(SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(myrank(v));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fsummary_cov_head", (DL_FUNC) &_fsummary_cov_head, 3},
    {"_fsummary_myrank", (DL_FUNC) &_fsummary_myrank, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_fsummary(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

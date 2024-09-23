// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sample_m
double sample_m(double n);
RcppExport SEXP _NMLmulti_sample_m(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_m(n));
    return rcpp_result_gen;
END_RCPP
}
// reject
Rcpp::NumericVector reject(double n);
RcppExport SEXP _NMLmulti_reject(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(reject(n));
    return rcpp_result_gen;
END_RCPP
}
// gibbs
Rcpp::NumericVector gibbs(NumericVector n, NumericVector ks);
RcppExport SEXP _NMLmulti_gibbs(SEXP nSEXP, SEXP ksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ks(ksSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs(n, ks));
    return rcpp_result_gen;
END_RCPP
}
// burnin
Rcpp::NumericVector burnin(int burn, NumericVector init, NumericVector ks, NumericVector Ns);
RcppExport SEXP _NMLmulti_burnin(SEXP burnSEXP, SEXP initSEXP, SEXP ksSEXP, SEXP NsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type init(initSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ks(ksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ns(NsSEXP);
    rcpp_result_gen = Rcpp::wrap(burnin(burn, init, ks, Ns));
    return rcpp_result_gen;
END_RCPP
}
// gen_chain
Rcpp::NumericMatrix gen_chain(NumericVector ks, NumericVector Ns, int batchsize, int thin, NumericVector startvec);
RcppExport SEXP _NMLmulti_gen_chain(SEXP ksSEXP, SEXP NsSEXP, SEXP batchsizeSEXP, SEXP thinSEXP, SEXP startvecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type ks(ksSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< int >::type batchsize(batchsizeSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type startvec(startvecSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_chain(ks, Ns, batchsize, thin, startvec));
    return rcpp_result_gen;
END_RCPP
}
// SCR
double SCR(NumericVector Q0, NumericVector chd, NumericVector NN);
RcppExport SEXP _NMLmulti_SCR(SEXP Q0SEXP, SEXP chdSEXP, SEXP NNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q0(Q0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chd(chdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NN(NNSEXP);
    rcpp_result_gen = Rcpp::wrap(SCR(Q0, chd, NN));
    return rcpp_result_gen;
END_RCPP
}
// LT
double LT(NumericVector Q0, NumericVector chd, NumericVector NN);
RcppExport SEXP _NMLmulti_LT(SEXP Q0SEXP, SEXP chdSEXP, SEXP NNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q0(Q0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chd(chdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NN(NNSEXP);
    rcpp_result_gen = Rcpp::wrap(LT(Q0, chd, NN));
    return rcpp_result_gen;
END_RCPP
}
// SCR_V
double SCR_V(NumericVector Q0, NumericVector chd, NumericVector NN);
RcppExport SEXP _NMLmulti_SCR_V(SEXP Q0SEXP, SEXP chdSEXP, SEXP NNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q0(Q0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chd(chdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NN(NNSEXP);
    rcpp_result_gen = Rcpp::wrap(SCR_V(Q0, chd, NN));
    return rcpp_result_gen;
END_RCPP
}
// SCR_G
double SCR_G(NumericVector Q0, NumericVector chd, NumericVector NN);
RcppExport SEXP _NMLmulti_SCR_G(SEXP Q0SEXP, SEXP chdSEXP, SEXP NNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q0(Q0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chd(chdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NN(NNSEXP);
    rcpp_result_gen = Rcpp::wrap(SCR_G(Q0, chd, NN));
    return rcpp_result_gen;
END_RCPP
}
// ifmu
double ifmu(NumericVector Q0, NumericVector chd, NumericVector NN);
RcppExport SEXP _NMLmulti_ifmu(SEXP Q0SEXP, SEXP chdSEXP, SEXP NNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q0(Q0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chd(chdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NN(NNSEXP);
    rcpp_result_gen = Rcpp::wrap(ifmu(Q0, chd, NN));
    return rcpp_result_gen;
END_RCPP
}
// ufmu
double ufmu(NumericVector Q0, NumericVector chd, NumericVector NN);
RcppExport SEXP _NMLmulti_ufmu(SEXP Q0SEXP, SEXP chdSEXP, SEXP NNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q0(Q0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chd(chdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NN(NNSEXP);
    rcpp_result_gen = Rcpp::wrap(ufmu(Q0, chd, NN));
    return rcpp_result_gen;
END_RCPP
}
// qm
double qm(NumericVector Q0, NumericVector chd, NumericVector NN);
RcppExport SEXP _NMLmulti_qm(SEXP Q0SEXP, SEXP chdSEXP, SEXP NNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type Q0(Q0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chd(chdSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type NN(NNSEXP);
    rcpp_result_gen = Rcpp::wrap(qm(Q0, chd, NN));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NMLmulti_sample_m", (DL_FUNC) &_NMLmulti_sample_m, 1},
    {"_NMLmulti_reject", (DL_FUNC) &_NMLmulti_reject, 1},
    {"_NMLmulti_gibbs", (DL_FUNC) &_NMLmulti_gibbs, 2},
    {"_NMLmulti_burnin", (DL_FUNC) &_NMLmulti_burnin, 4},
    {"_NMLmulti_gen_chain", (DL_FUNC) &_NMLmulti_gen_chain, 5},
    {"_NMLmulti_SCR", (DL_FUNC) &_NMLmulti_SCR, 3},
    {"_NMLmulti_LT", (DL_FUNC) &_NMLmulti_LT, 3},
    {"_NMLmulti_SCR_V", (DL_FUNC) &_NMLmulti_SCR_V, 3},
    {"_NMLmulti_SCR_G", (DL_FUNC) &_NMLmulti_SCR_G, 3},
    {"_NMLmulti_ifmu", (DL_FUNC) &_NMLmulti_ifmu, 3},
    {"_NMLmulti_ufmu", (DL_FUNC) &_NMLmulti_ufmu, 3},
    {"_NMLmulti_qm", (DL_FUNC) &_NMLmulti_qm, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_NMLmulti(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "longBet_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// predict
Rcpp::List predict(arma::mat X, arma::mat t, Rcpp::XPtr<std::vector<std::vector<tree>>> tree_pnt);
RcppExport SEXP _longBet_predict(SEXP XSEXP, SEXP tSEXP, SEXP tree_pntSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<std::vector<std::vector<tree>>> >::type tree_pnt(tree_pntSEXP);
    rcpp_result_gen = Rcpp::wrap(predict(X, t, tree_pnt));
    return rcpp_result_gen;
END_RCPP
}
// predict_beta
Rcpp::List predict_beta(arma::mat t_test, arma::mat t_train, arma::mat res, arma::mat A_diag, arma::mat Sig_diag, double sig_knl, double lambda_knl);
RcppExport SEXP _longBet_predict_beta(SEXP t_testSEXP, SEXP t_trainSEXP, SEXP resSEXP, SEXP A_diagSEXP, SEXP Sig_diagSEXP, SEXP sig_knlSEXP, SEXP lambda_knlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type t_test(t_testSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_train(t_trainSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type res(resSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A_diag(A_diagSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sig_diag(Sig_diagSEXP);
    Rcpp::traits::input_parameter< double >::type sig_knl(sig_knlSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_knl(lambda_knlSEXP);
    rcpp_result_gen = Rcpp::wrap(predict_beta(t_test, t_train, res, A_diag, Sig_diag, sig_knl, lambda_knl));
    return rcpp_result_gen;
END_RCPP
}
// r_to_json
Rcpp::StringVector r_to_json(double y_mean, Rcpp::XPtr<std::vector<std::vector<tree>>> tree_pnt);
RcppExport SEXP _longBet_r_to_json(SEXP y_meanSEXP, SEXP tree_pntSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type y_mean(y_meanSEXP);
    Rcpp::traits::input_parameter< Rcpp::XPtr<std::vector<std::vector<tree>>> >::type tree_pnt(tree_pntSEXP);
    rcpp_result_gen = Rcpp::wrap(r_to_json(y_mean, tree_pnt));
    return rcpp_result_gen;
END_RCPP
}
// json_to_r
Rcpp::List json_to_r(Rcpp::StringVector json_string_r);
RcppExport SEXP _longBet_json_to_r(SEXP json_string_rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type json_string_r(json_string_rSEXP);
    rcpp_result_gen = Rcpp::wrap(json_to_r(json_string_r));
    return rcpp_result_gen;
END_RCPP
}
// sample_int_crank
IntegerVector sample_int_crank(int n, int size, NumericVector prob);
RcppExport SEXP _longBet_sample_int_crank(SEXP nSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int_crank(n, size, prob));
    return rcpp_result_gen;
END_RCPP
}
// sample_int_ccrank
SEXP sample_int_ccrank(int n, int size, NumericVector prob);
RcppExport SEXP _longBet_sample_int_ccrank(SEXP nSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int_ccrank(n, size, prob));
    return rcpp_result_gen;
END_RCPP
}
// sample_int_expj
IntegerVector sample_int_expj(int n, int size, NumericVector prob);
RcppExport SEXP _longBet_sample_int_expj(SEXP nSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int_expj(n, size, prob));
    return rcpp_result_gen;
END_RCPP
}
// sample_int_expjs
IntegerVector sample_int_expjs(int n, int size, NumericVector prob);
RcppExport SEXP _longBet_sample_int_expjs(SEXP nSEXP, SEXP sizeSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type size(sizeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_int_expjs(n, size, prob));
    return rcpp_result_gen;
END_RCPP
}
// longBet_cpp
Rcpp::List longBet_cpp(arma::mat y, arma::mat X, arma::mat X_tau, arma::mat z, arma::mat t_con, arma::mat t_mod, arma::mat post_t, arma::mat trt_time, size_t num_sweeps, size_t burnin, size_t max_depth, size_t n_min, size_t num_cutpoints, double no_split_penality, size_t mtry_pr, size_t mtry_trt, size_t p_categorical_pr, size_t p_categorical_trt, size_t num_trees_pr, double alpha_pr, double beta_pr, double tau_pr, double kap_pr, double s_pr, bool pr_scale, size_t num_trees_trt, double alpha_trt, double beta_trt, double tau_trt, double kap_trt, double s_trt, bool trt_scale, bool verbose, bool parallel, bool set_random_seed, size_t random_seed, bool sample_weights_flag, bool a_scaling, bool b_scaling, bool split_time_ps, bool split_time_trt, double sig_knl, double lambda_knl);
RcppExport SEXP _longBet_longBet_cpp(SEXP ySEXP, SEXP XSEXP, SEXP X_tauSEXP, SEXP zSEXP, SEXP t_conSEXP, SEXP t_modSEXP, SEXP post_tSEXP, SEXP trt_timeSEXP, SEXP num_sweepsSEXP, SEXP burninSEXP, SEXP max_depthSEXP, SEXP n_minSEXP, SEXP num_cutpointsSEXP, SEXP no_split_penalitySEXP, SEXP mtry_prSEXP, SEXP mtry_trtSEXP, SEXP p_categorical_prSEXP, SEXP p_categorical_trtSEXP, SEXP num_trees_prSEXP, SEXP alpha_prSEXP, SEXP beta_prSEXP, SEXP tau_prSEXP, SEXP kap_prSEXP, SEXP s_prSEXP, SEXP pr_scaleSEXP, SEXP num_trees_trtSEXP, SEXP alpha_trtSEXP, SEXP beta_trtSEXP, SEXP tau_trtSEXP, SEXP kap_trtSEXP, SEXP s_trtSEXP, SEXP trt_scaleSEXP, SEXP verboseSEXP, SEXP parallelSEXP, SEXP set_random_seedSEXP, SEXP random_seedSEXP, SEXP sample_weights_flagSEXP, SEXP a_scalingSEXP, SEXP b_scalingSEXP, SEXP split_time_psSEXP, SEXP split_time_trtSEXP, SEXP sig_knlSEXP, SEXP lambda_knlSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_tau(X_tauSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_con(t_conSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_mod(t_modSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type post_t(post_tSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type trt_time(trt_timeSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_sweeps(num_sweepsSEXP);
    Rcpp::traits::input_parameter< size_t >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< size_t >::type max_depth(max_depthSEXP);
    Rcpp::traits::input_parameter< size_t >::type n_min(n_minSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_cutpoints(num_cutpointsSEXP);
    Rcpp::traits::input_parameter< double >::type no_split_penality(no_split_penalitySEXP);
    Rcpp::traits::input_parameter< size_t >::type mtry_pr(mtry_prSEXP);
    Rcpp::traits::input_parameter< size_t >::type mtry_trt(mtry_trtSEXP);
    Rcpp::traits::input_parameter< size_t >::type p_categorical_pr(p_categorical_prSEXP);
    Rcpp::traits::input_parameter< size_t >::type p_categorical_trt(p_categorical_trtSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_trees_pr(num_trees_prSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_pr(alpha_prSEXP);
    Rcpp::traits::input_parameter< double >::type beta_pr(beta_prSEXP);
    Rcpp::traits::input_parameter< double >::type tau_pr(tau_prSEXP);
    Rcpp::traits::input_parameter< double >::type kap_pr(kap_prSEXP);
    Rcpp::traits::input_parameter< double >::type s_pr(s_prSEXP);
    Rcpp::traits::input_parameter< bool >::type pr_scale(pr_scaleSEXP);
    Rcpp::traits::input_parameter< size_t >::type num_trees_trt(num_trees_trtSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_trt(alpha_trtSEXP);
    Rcpp::traits::input_parameter< double >::type beta_trt(beta_trtSEXP);
    Rcpp::traits::input_parameter< double >::type tau_trt(tau_trtSEXP);
    Rcpp::traits::input_parameter< double >::type kap_trt(kap_trtSEXP);
    Rcpp::traits::input_parameter< double >::type s_trt(s_trtSEXP);
    Rcpp::traits::input_parameter< bool >::type trt_scale(trt_scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type parallel(parallelSEXP);
    Rcpp::traits::input_parameter< bool >::type set_random_seed(set_random_seedSEXP);
    Rcpp::traits::input_parameter< size_t >::type random_seed(random_seedSEXP);
    Rcpp::traits::input_parameter< bool >::type sample_weights_flag(sample_weights_flagSEXP);
    Rcpp::traits::input_parameter< bool >::type a_scaling(a_scalingSEXP);
    Rcpp::traits::input_parameter< bool >::type b_scaling(b_scalingSEXP);
    Rcpp::traits::input_parameter< bool >::type split_time_ps(split_time_psSEXP);
    Rcpp::traits::input_parameter< bool >::type split_time_trt(split_time_trtSEXP);
    Rcpp::traits::input_parameter< double >::type sig_knl(sig_knlSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_knl(lambda_knlSEXP);
    rcpp_result_gen = Rcpp::wrap(longBet_cpp(y, X, X_tau, z, t_con, t_mod, post_t, trt_time, num_sweeps, burnin, max_depth, n_min, num_cutpoints, no_split_penality, mtry_pr, mtry_trt, p_categorical_pr, p_categorical_trt, num_trees_pr, alpha_pr, beta_pr, tau_pr, kap_pr, s_pr, pr_scale, num_trees_trt, alpha_trt, beta_trt, tau_trt, kap_trt, s_trt, trt_scale, verbose, parallel, set_random_seed, random_seed, sample_weights_flag, a_scaling, b_scaling, split_time_ps, split_time_trt, sig_knl, lambda_knl));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_longBet_predict", (DL_FUNC) &_longBet_predict, 3},
    {"_longBet_predict_beta", (DL_FUNC) &_longBet_predict_beta, 7},
    {"_longBet_r_to_json", (DL_FUNC) &_longBet_r_to_json, 2},
    {"_longBet_json_to_r", (DL_FUNC) &_longBet_json_to_r, 1},
    {"_longBet_sample_int_crank", (DL_FUNC) &_longBet_sample_int_crank, 3},
    {"_longBet_sample_int_ccrank", (DL_FUNC) &_longBet_sample_int_ccrank, 3},
    {"_longBet_sample_int_expj", (DL_FUNC) &_longBet_sample_int_expj, 3},
    {"_longBet_sample_int_expjs", (DL_FUNC) &_longBet_sample_int_expjs, 3},
    {"_longBet_longBet_cpp", (DL_FUNC) &_longBet_longBet_cpp, 43},
    {NULL, NULL, 0}
};

RcppExport void R_init_longBet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

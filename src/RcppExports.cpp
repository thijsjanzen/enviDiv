// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// create_tree_cpp
List create_tree_cpp(NumericVector parameters, NumericVector waterlevel_changes, int seed, float crown_age, int max_lin);
RcppExport SEXP _enviDiv_create_tree_cpp(SEXP parametersSEXP, SEXP waterlevel_changesSEXP, SEXP seedSEXP, SEXP crown_ageSEXP, SEXP max_linSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type waterlevel_changes(waterlevel_changesSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< float >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    rcpp_result_gen = Rcpp::wrap(create_tree_cpp(parameters, waterlevel_changes, seed, crown_age, max_lin));
    return rcpp_result_gen;
END_RCPP
}
// create_ref_table_cpp
List create_ref_table_cpp(int model, int num_repl, float crown_age, int min_lin, int max_lin, int num_threads);
RcppExport SEXP _enviDiv_create_ref_table_cpp(SEXP modelSEXP, SEXP num_replSEXP, SEXP crown_ageSEXP, SEXP min_linSEXP, SEXP max_linSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< int >::type num_repl(num_replSEXP);
    Rcpp::traits::input_parameter< float >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type min_lin(min_linSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(create_ref_table_cpp(model, num_repl, crown_age, min_lin, max_lin, num_threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_enviDiv_create_tree_cpp", (DL_FUNC) &_enviDiv_create_tree_cpp, 5},
    {"_enviDiv_create_ref_table_cpp", (DL_FUNC) &_enviDiv_create_ref_table_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_enviDiv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// create_tree_cpp
List create_tree_cpp(std::vector<float> parameters, std::vector<float> waterlevel_changes, int seed, float crown_age, int max_lin);
RcppExport SEXP _enviDiv_create_tree_cpp(SEXP parametersSEXP, SEXP waterlevel_changesSEXP, SEXP seedSEXP, SEXP crown_ageSEXP, SEXP max_linSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<float> >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< std::vector<float> >::type waterlevel_changes(waterlevel_changesSEXP);
    Rcpp::traits::input_parameter< int >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< float >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    rcpp_result_gen = Rcpp::wrap(create_tree_cpp(parameters, waterlevel_changes, seed, crown_age, max_lin));
    return rcpp_result_gen;
END_RCPP
}
// test_envidiv_tbb
Rcpp::List test_envidiv_tbb(int model, float crown_age);
RcppExport SEXP _enviDiv_test_envidiv_tbb(SEXP modelSEXP, SEXP crown_ageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< float >::type crown_age(crown_ageSEXP);
    rcpp_result_gen = Rcpp::wrap(test_envidiv_tbb(model, crown_age));
    return rcpp_result_gen;
END_RCPP
}
// create_ref_table_tbb_serial
List create_ref_table_tbb_serial(int model, int num_repl, float crown_age, int min_lin, int max_lin, int num_threads);
RcppExport SEXP _enviDiv_create_ref_table_tbb_serial(SEXP modelSEXP, SEXP num_replSEXP, SEXP crown_ageSEXP, SEXP min_linSEXP, SEXP max_linSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< int >::type num_repl(num_replSEXP);
    Rcpp::traits::input_parameter< float >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type min_lin(min_linSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(create_ref_table_tbb_serial(model, num_repl, crown_age, min_lin, max_lin, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// create_ref_table_tbb_par
List create_ref_table_tbb_par(int model, int num_repl, float crown_age, int min_lin, int max_lin, int num_threads);
RcppExport SEXP _enviDiv_create_ref_table_tbb_par(SEXP modelSEXP, SEXP num_replSEXP, SEXP crown_ageSEXP, SEXP min_linSEXP, SEXP max_linSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< int >::type num_repl(num_replSEXP);
    Rcpp::traits::input_parameter< float >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type min_lin(min_linSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(create_ref_table_tbb_par(model, num_repl, crown_age, min_lin, max_lin, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// sq_numbers_cpp_tbb
NumericVector sq_numbers_cpp_tbb(int n, int num_threads);
RcppExport SEXP _enviDiv_sq_numbers_cpp_tbb(SEXP nSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(sq_numbers_cpp_tbb(n, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// param_from_prior_cpp
std::vector<float> param_from_prior_cpp();
RcppExport SEXP _enviDiv_param_from_prior_cpp() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(param_from_prior_cpp());
    return rcpp_result_gen;
END_RCPP
}
// param_from_prior_exp_cpp
std::vector<float> param_from_prior_exp_cpp();
RcppExport SEXP _enviDiv_param_from_prior_exp_cpp() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(param_from_prior_exp_cpp());
    return rcpp_result_gen;
END_RCPP
}
// get_waterlevel_cpp
std::vector<float> get_waterlevel_cpp(int water_model, float maximum_time);
RcppExport SEXP _enviDiv_get_waterlevel_cpp(SEXP water_modelSEXP, SEXP maximum_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type water_model(water_modelSEXP);
    Rcpp::traits::input_parameter< float >::type maximum_time(maximum_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(get_waterlevel_cpp(water_model, maximum_time));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_enviDiv_create_tree_cpp", (DL_FUNC) &_enviDiv_create_tree_cpp, 5},
    {"_enviDiv_test_envidiv_tbb", (DL_FUNC) &_enviDiv_test_envidiv_tbb, 2},
    {"_enviDiv_create_ref_table_tbb_serial", (DL_FUNC) &_enviDiv_create_ref_table_tbb_serial, 6},
    {"_enviDiv_create_ref_table_tbb_par", (DL_FUNC) &_enviDiv_create_ref_table_tbb_par, 6},
    {"_enviDiv_sq_numbers_cpp_tbb", (DL_FUNC) &_enviDiv_sq_numbers_cpp_tbb, 2},
    {"_enviDiv_param_from_prior_cpp", (DL_FUNC) &_enviDiv_param_from_prior_cpp, 0},
    {"_enviDiv_param_from_prior_exp_cpp", (DL_FUNC) &_enviDiv_param_from_prior_exp_cpp, 0},
    {"_enviDiv_get_waterlevel_cpp", (DL_FUNC) &_enviDiv_get_waterlevel_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_enviDiv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

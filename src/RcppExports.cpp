// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// create_tree_cpp
List create_tree_cpp(std::vector<double> parameters, std::vector<double> waterlevel_changes, double crown_age, int max_lin);
RcppExport SEXP _enviDiv_create_tree_cpp(SEXP parametersSEXP, SEXP waterlevel_changesSEXP, SEXP crown_ageSEXP, SEXP max_linSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type waterlevel_changes(waterlevel_changesSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    rcpp_result_gen = Rcpp::wrap(create_tree_cpp(parameters, waterlevel_changes, crown_age, max_lin));
    return rcpp_result_gen;
END_RCPP
}
// test_envidiv_tbb
Rcpp::List test_envidiv_tbb(int model, double crown_age);
RcppExport SEXP _enviDiv_test_envidiv_tbb(SEXP modelSEXP, SEXP crown_ageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    rcpp_result_gen = Rcpp::wrap(test_envidiv_tbb(model, crown_age));
    return rcpp_result_gen;
END_RCPP
}
// create_ref_table_tbb_serial
List create_ref_table_tbb_serial(int model, int num_repl, double crown_age, int min_lin, int max_lin, int num_threads);
RcppExport SEXP _enviDiv_create_ref_table_tbb_serial(SEXP modelSEXP, SEXP num_replSEXP, SEXP crown_ageSEXP, SEXP min_linSEXP, SEXP max_linSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< int >::type num_repl(num_replSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type min_lin(min_linSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    Rcpp::traits::input_parameter< int >::type num_threads(num_threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(create_ref_table_tbb_serial(model, num_repl, crown_age, min_lin, max_lin, num_threads));
    return rcpp_result_gen;
END_RCPP
}
// create_ref_table_tbb_par
List create_ref_table_tbb_par(int model, int num_repl, double crown_age, int min_lin, int max_lin, int num_threads);
RcppExport SEXP _enviDiv_create_ref_table_tbb_par(SEXP modelSEXP, SEXP num_replSEXP, SEXP crown_ageSEXP, SEXP min_linSEXP, SEXP max_linSEXP, SEXP num_threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< int >::type num_repl(num_replSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
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
// sim_envidiv_cpp
Rcpp::List sim_envidiv_cpp(int model, std::vector<double> parameters, double crown_age, int max_lin);
RcppExport SEXP _enviDiv_sim_envidiv_cpp(SEXP modelSEXP, SEXP parametersSEXP, SEXP crown_ageSEXP, SEXP max_linSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    rcpp_result_gen = Rcpp::wrap(sim_envidiv_cpp(model, parameters, crown_age, max_lin));
    return rcpp_result_gen;
END_RCPP
}
// initial_draw_from_prior
Rcpp::NumericMatrix initial_draw_from_prior(int num_particles, double crown_age, int min_lin, int max_lin, bool verbose);
RcppExport SEXP _enviDiv_initial_draw_from_prior(SEXP num_particlesSEXP, SEXP crown_ageSEXP, SEXP min_linSEXP, SEXP max_linSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type num_particles(num_particlesSEXP);
    Rcpp::traits::input_parameter< double >::type crown_age(crown_ageSEXP);
    Rcpp::traits::input_parameter< int >::type min_lin(min_linSEXP);
    Rcpp::traits::input_parameter< int >::type max_lin(max_linSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(initial_draw_from_prior(num_particles, crown_age, min_lin, max_lin, verbose));
    return rcpp_result_gen;
END_RCPP
}
// param_from_prior_cpp
std::vector<double> param_from_prior_cpp(int model);
RcppExport SEXP _enviDiv_param_from_prior_cpp(SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(param_from_prior_cpp(model));
    return rcpp_result_gen;
END_RCPP
}
// param_from_prior_exp_cpp
std::vector<double> param_from_prior_exp_cpp(int model);
RcppExport SEXP _enviDiv_param_from_prior_exp_cpp(SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(param_from_prior_exp_cpp(model));
    return rcpp_result_gen;
END_RCPP
}
// get_waterlevel_cpp
std::vector<double> get_waterlevel_cpp(int water_model, double maximum_time, double rate);
RcppExport SEXP _enviDiv_get_waterlevel_cpp(SEXP water_modelSEXP, SEXP maximum_timeSEXP, SEXP rateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type water_model(water_modelSEXP);
    Rcpp::traits::input_parameter< double >::type maximum_time(maximum_timeSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    rcpp_result_gen = Rcpp::wrap(get_waterlevel_cpp(water_model, maximum_time, rate));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_enviDiv_create_tree_cpp", (DL_FUNC) &_enviDiv_create_tree_cpp, 4},
    {"_enviDiv_test_envidiv_tbb", (DL_FUNC) &_enviDiv_test_envidiv_tbb, 2},
    {"_enviDiv_create_ref_table_tbb_serial", (DL_FUNC) &_enviDiv_create_ref_table_tbb_serial, 6},
    {"_enviDiv_create_ref_table_tbb_par", (DL_FUNC) &_enviDiv_create_ref_table_tbb_par, 6},
    {"_enviDiv_sq_numbers_cpp_tbb", (DL_FUNC) &_enviDiv_sq_numbers_cpp_tbb, 2},
    {"_enviDiv_sim_envidiv_cpp", (DL_FUNC) &_enviDiv_sim_envidiv_cpp, 4},
    {"_enviDiv_initial_draw_from_prior", (DL_FUNC) &_enviDiv_initial_draw_from_prior, 5},
    {"_enviDiv_param_from_prior_cpp", (DL_FUNC) &_enviDiv_param_from_prior_cpp, 1},
    {"_enviDiv_param_from_prior_exp_cpp", (DL_FUNC) &_enviDiv_param_from_prior_exp_cpp, 1},
    {"_enviDiv_get_waterlevel_cpp", (DL_FUNC) &_enviDiv_get_waterlevel_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_enviDiv(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include "Rcpp.h"
#include "random_thijs.h"

Rcpp::NumericVector param_from_prior_cpp();
Rcpp::NumericVector param_from_prior_exp_cpp();
std::vector<float> get_waterlevel_cpp(int water_model,
                                      float maximum_time);

#endif

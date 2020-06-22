#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include "Rcpp.h"
#include "random_thijs.h"

Rcpp::NumericVector param_from_prior_cpp();
Rcpp::NumericVector param_from_prior_exp_cpp();

#endif

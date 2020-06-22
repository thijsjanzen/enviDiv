#include <cmath>
#include "random_thijs.h"
#include <string>

#include <Rcpp.h>
using namespace Rcpp;

//' draw parameter combinations from the prior
//' @return vector with 6 entries: extinction, sympatric speciation at high
//' water, sympatric speciation at low water, allopatric speciation, amount of
//' perturbation, the chosen water model
//' @export
// [[Rcpp::export]]
NumericVector param_from_prior_cpp() {

  rnd_t rndgen;

  NumericVector output(6);

  output[0] = powf(10, (-3 + 5 * rndgen.uniform()));  // extinction
  output[1] = powf(10, (-3 + 5 * rndgen.uniform()));  // symp spec high
  output[2] = powf(10, (-3 + 5 * rndgen.uniform()));  // symp spec low
  output[3] = powf(10, (-3 + 5 * rndgen.uniform()));  // allo spec
  output[4] = powf(10, (-3 + 3 * rndgen.uniform()));  // jiggle
  output[5] = 1 + rndgen.random_number(3); // model

  return(output);
}

//' draw parameter combinations from the prior
//' @return vector with 6  entries: extinction, sympatric speciation at high
//' water, sympatric speciation at low water, allopatric speciation, amount of
//' perturbation, the chosen water model
//' @export
// [[Rcpp::export]]
NumericVector param_from_prior_exp_cpp() {
  rnd_t rndgen;

  NumericVector output(6);
  output[0] = rndgen.Expon(0.05); // extinction
  output[1] = rndgen.Expon(0.1);  // symp spec high
  output[2] = rndgen.Expon( 0.1);  // symp spec low
  output[3] = rndgen.Expon(0.1);  // allo spec
  output[4] = rndgen.Expon( 0.05);  // jiggle
  output[5] = 1 + rndgen.random_number(3);

  return(output);
}


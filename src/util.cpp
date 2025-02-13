#include <cmath>
#include "random_thijs.h"
#include <string>
#include "Gillespie.h"
#include "util.h"

#include <Rcpp.h>
using namespace Rcpp;



std::vector<double> parameters_from_prior(rnd_t& rndgen_) {
  std::vector<double> output(7, 0.f);

  output[param_type::extinction_rate] = powf(10, (-3 + 5 * rndgen_.uniform()));  // extinction
  output[param_type::sym_high_rate]   = powf(10, (-3 + 5 * rndgen_.uniform()));  // symp spec high
  output[param_type::sym_low_rate]    = powf(10, (-3 + 5 * rndgen_.uniform()));  // symp spec low
  output[param_type::allo_rate]       = powf(10, (-3 + 5 * rndgen_.uniform()));  // allo spec
  output[param_type::wobble_rate]     = powf(10, (-3 + 3 * rndgen_.uniform()));  // jiggle
  output[param_type::water_rate]      = powf(10, ( 0 + 3 * rndgen_.uniform()));  // water rate
  output[param_type::model]           = 1 + rndgen_.random_number(4); // model

  return(output);
}

std::vector<double> parameters_from_prior(rnd_t& rndgen_,
                                          int model) {
  std::vector<double> output = parameters_from_prior(rndgen_);
  output[param_type::model] = model;

  return(output);
}

std::vector<double> get_waterlevel_changes(int water_model,
                                          double maximum_time,
                                          rnd_t& rndgen_,
                                          double rate) {
  if (water_model == 1) {
    std::vector<double> output(2);
    output[0] = 0.0;
    output[1] = maximum_time * 2;
    return output;
  }

  if (water_model == 4) { // fake testing model
    std::vector< double > output = {0.0};
    int water_level = 1;
    double time = 0;


    while (time < maximum_time) {
      time += rndgen_.Expon(rate);
      if (time > maximum_time) break;

      output.push_back(time);
      water_level = 1 - water_level;
    }

    return output;
  }

  if (water_model == 40) { // fake testing model
    std::vector<double> output;
    output.push_back(0.0);
    output.push_back(maximum_time * 0.5);
    output.push_back(maximum_time * 0.55);
    output.push_back(maximum_time * 2);
    return output;
  }

  std::vector< double > output;

  if (water_model == 3) {
    int water_level = 1;
    double time = 0;
    double temp_time = time;

    while (temp_time < (maximum_time - 1.1)) {
      temp_time = temp_time + rndgen_.Expon(rate);
      if (temp_time > (maximum_time - 1.1)) break;

      time = temp_time;
      output.push_back(time);
      water_level = 1 - water_level;
    }
    if (water_level == 0) {
      output.pop_back(); // water level has to be high
    }
  }

  output.push_back(maximum_time - 1.1)	;		// 0
  output.push_back(maximum_time - 0.55)	;	// 1
  output.push_back(maximum_time - 0.393);		// 0
  output.push_back(maximum_time - 0.363);		// 1
  output.push_back(maximum_time - 0.295);		// 0
  output.push_back(maximum_time - 0.262);		// 1
  output.push_back(maximum_time - 0.193);		// 0
  output.push_back(maximum_time - 0.169);		// 1
  output.push_back(maximum_time - 0.04)	;	// 0
  output.push_back(maximum_time - 0.035);	// 1
  output.push_back(maximum_time)				; // 1

  return(output);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////// some cpp functions to retain functionality in R  //////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//' draw parameter combinations from the prior
//' @return vector with 6 entries: extinction, sympatric speciation at high
//' water, sympatric speciation at low water, allopatric speciation, amount of
//' perturbation, the chosen water model
//' @export
// [[Rcpp::export]]
std::vector<double> param_from_prior_cpp() {

  rnd_t rndgen;

  std::vector<double> output = parameters_from_prior(rndgen);

  return(output);
}

//' draw parameter combinations from the prior
//' @return vector with 6  entries: extinction, sympatric speciation at high
//' water, sympatric speciation at low water, allopatric speciation, amount of
//' perturbation, the chosen water model
//' @export
// [[Rcpp::export]]
std::vector<double> param_from_prior_exp_cpp() {
  rnd_t rndgen;

  std::vector<double> output(6, 0.f);
  output[param_type::extinction_rate] = rndgen.Expon(0.05); // extinction
  output[param_type::sym_high_rate] = rndgen.Expon(0.1);  // symp spec high
  output[param_type::sym_low_rate] = rndgen.Expon( 0.1);  // symp spec low
  output[param_type::allo_rate] = rndgen.Expon(0.1);  // allo spec
  output[param_type::wobble_rate] = rndgen.Expon( 0.05);  // jiggle
  output[param_type::water_rate] = rndgen.Expon(10); // water rate
  output[param_type::model] = 1 + rndgen.random_number(4);

  return(output);
}

//' draw water level change2
//' @param water_model water model
//' @param maximum_time crown age
//' @param water level change rate per MY, default is 10
//' @return water level changes
//' @export
// [[Rcpp::export]]
std::vector<double> get_waterlevel_cpp(int water_model,
                                      double maximum_time,
                                      double rate = 10) {

  rnd_t rndgen;

  std::vector<double> output = get_waterlevel_changes(water_model,
                                                      maximum_time,
                                                      rndgen,
                                                      rate);

  return(output);
}

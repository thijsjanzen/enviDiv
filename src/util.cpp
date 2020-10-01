#include <cmath>
#include "random_thijs.h"
#include <string>
#include "Gillespie.h"
#include "util.h"

#include <Rcpp.h>
using namespace Rcpp;


std::vector<float> parameters_from_prior(rnd_t& rndgen_) {
  std::vector<float> output(6, 0.f);
/*
  output[0] = powf(10, (-5 + 6 * rndgen_.uniform()));  // extinction
  output[1] = powf(10, (-5 + 6 * rndgen_.uniform()));  // symp spec high
  output[2] = powf(10, (-5 + 6 * rndgen_.uniform()));  // symp spec low
  output[3] = powf(10, (-5 + 6 * rndgen_.uniform()));  // allo spec
  output[4] = powf(10, (-5 + 6 * rndgen_.uniform()));  // jiggle
  output[5] = 1 + rndgen_.random_number(3); // model
*/
  output[0] = rndgen_.Expon(10);
  output[1] = rndgen_.Expon(10);
  output[2] = rndgen_.Expon(10);
  output[3] = rndgen_.Expon(10);
  output[4] = rndgen_.Expon(10);
  output[5] = 1 + rndgen_.random_number(3); // model

  return(output);
}

std::vector<float> parameters_from_prior(rnd_t& rndgen_,
                                         int model) {
  std::vector<float> output = parameters_from_prior(rndgen_);
  output[5] = model; // model

  return(output);
}


std::vector<float> get_waterlevel_changes(int water_model,
                                          float maximum_time,
                                          rnd_t& rndgen_) {
  if (water_model == 1) {
    std::vector<float> output(2);
    output[0] = 0.0f;
    output[1] = maximum_time * 2;
    return output;
  }

  std::vector< float > output;

  if (water_model == 3) {
    int water_level = 1;
    float time = 0;
    float temp_time = time;

    rnd_t rndgen;

    while (temp_time < (maximum_time - 1.1)) {
      temp_time = temp_time + rndgen.Expon(10);
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
std::vector<float> param_from_prior_cpp() {

  rnd_t rndgen;

  std::vector<float> output = parameters_from_prior(rndgen);

  return(output);
}

//' draw parameter combinations from the prior
//' @return vector with 6  entries: extinction, sympatric speciation at high
//' water, sympatric speciation at low water, allopatric speciation, amount of
//' perturbation, the chosen water model
//' @export
// [[Rcpp::export]]
std::vector<float> param_from_prior_exp_cpp() {
  rnd_t rndgen;

  std::vector<float> output(6, 0.f);
  output[0] = rndgen.Expon(20); // extinction
  output[1] = rndgen.Expon(10);  // symp spec high
  output[2] = rndgen.Expon(10);  // symp spec low
  output[3] = rndgen.Expon(10);  // allo spec
  output[4] = rndgen.Expon(20);  // jiggle
  output[5] = 1 + rndgen.random_number(3);

  return(output);
}

//' draw water level change2
//' @param water_model water model
//' @param maximum_time crown age
//' @return water level changes
//' @export
// [[Rcpp::export]]
std::vector<float> get_waterlevel_cpp(int water_model,
                                      float maximum_time) {

  rnd_t rndgen;

  std::vector<float> output = get_waterlevel_changes(water_model,
                                                     maximum_time,
                                                     rndgen);

  return(output);
}

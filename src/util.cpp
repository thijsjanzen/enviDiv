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


// [[Rcpp::export]]
std::vector<float> get_waterlevel_cpp(int water_model,
                                 float maximum_time) {
  if (water_model == 1) {
    std::vector<float> output(2);
    output[0] = 0.0;
    output[1] = maximum_time * 2;
    return output;
  }

  if (water_model == 2) {
      std::vector<float> output(11);
      output[0] = (maximum_time - 1.1)	;		// 0
      output[1] = (maximum_time - 0.55)	;	// 1
      output[2] = (maximum_time - 0.393);		// 0
      output[3] = (maximum_time - 0.363);		// 1
      output[4] = (maximum_time - 0.295);		// 0
      output[5] = (maximum_time - 0.262);		// 1
      output[6] = (maximum_time - 0.193);		// 0
      output[7] = (maximum_time - 0.169);		// 1
      output[8] = (maximum_time - 0.04)	;	// 0
      output[9] = (maximum_time - 0.035);	// 1
      output[10] = (maximum_time)				; // 1
      return(output);
  }

  if (water_model == 3) {
    std::vector<float> output2;
    int water_level = 1;
    float time = 0;
    float temp_time = time;

    rnd_t rndgen;

    while (temp_time < (maximum_time - 1.1)) {
      temp_time = temp_time + rndgen.Expon(10);
      if (temp_time > (maximum_time - 1.1)) break;

      time = temp_time;
      output2.push_back(time);
      water_level = 1 - water_level;
    }
    if (water_level == 0) {
      output2.pop_back(); // water level has to be high
    }

    output2.push_back(maximum_time - 1.1)		;  // 0
    output2.push_back(maximum_time - 0.55)  ;	 // 1
    output2.push_back(maximum_time - 0.393)	;  // 0
    output2.push_back(maximum_time - 0.363)	;  // 1
    output2.push_back(maximum_time - 0.295)	;  // 0
    output2.push_back(maximum_time - 0.262)	;  // 1
    output2.push_back(maximum_time - 0.193)	;  // 0
    output2.push_back(maximum_time - 0.169)	;  // 1
    output2.push_back(maximum_time - 0.04)	;	 // 0
    output2.push_back(maximum_time - 0.035)	;  // 1
    output2.push_back(maximum_time)		;		     // 1

    return(output2);
  }

}





#include <cmath>
#include "random_thijs.h"
#include <string>
#include "Gillespie.h"
#include "util.h"

#include <Rcpp.h>
using namespace Rcpp;

void make_sleep(size_t ms) {
  std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}


std::vector<float> parameters_from_prior(rnd_t& rndgen_,
                                         int model = -1) {
  std::vector<float> output(6, 0.f);

  output[5] = model;
  if (model == -1) output[5] = 1 + rndgen_.random_number(3); // model

  if (output[5] == 1) {                                  // birth death model
    output[0] = powf(10, (-4 + 6 * rndgen_.uniform()));  // extinction
    output[1] = powf(10, (-4 + 6 * rndgen_.uniform()));  // symp spec high
    output[2] = powf(10, (-4 + 6 * rndgen_.uniform()));  // symp spec low
    output[3] = powf(10, (-4 + 6 * rndgen_.uniform()));  // allo spec
    output[4] = powf(10, (-4 + 6 * rndgen_.uniform()));  // jiggle
  }

  if (output[5] != 1) {                                  // water models
    output[0] = powf(10, (-4 + 6 * rndgen_.uniform()));  // extinction
    output[1] = powf(10, (-4 + 6 * rndgen_.uniform()));  // symp spec high
    output[2] = powf(10, (-2 + 4 * rndgen_.uniform()));  // symp spec low
    output[3] = powf(10, (-2 + 4 * rndgen_.uniform()));  // allo spec
    output[4] = powf(10, (-4 + 6 * rndgen_.uniform()));  // jiggle
  }

  return(output);
}

/*
std::vector<float> parameters_from_prior(rnd_t& rndgen_,
                                         int model) {
  std::vector<float> output = parameters_from_prior(rndgen_);
  output[5] = model; // model

  return(output);
                                         }*/


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


std::vector< std::vector< float >> create_l_table_float(
    const std::vector< species > & s1,
    const std::vector< species > & s2) {

//  Rcout << s1.size() << " " << s2.size() << "\n";
//  make_sleep(1);
//  force_output("");

  std::vector< std::vector< float > > output;
  int num_rows = s1.size() + s2.size();
  int num_cols = 4;
  for(int i = 0; i < num_rows; ++i) {
    std::vector< float > row(num_cols);
    output.push_back(row);
  }

  int i = 0;
  for (auto it = s1.begin(); it != s1.end(); ++it, ++i) {
    output[i][0] = (*it).birth_time;
    output[i][1] = (*it).get_parent();
    output[i][2] = (*it).get_ID();
    output[i][3] = (*it).death_time;
  }

  for (auto it = s2.begin(); it != s2.end(); ++it, ++i) {
    output[i][0]  = (*it).birth_time;
    output[i][1]  = -1 * ((*it).get_parent());
    output[i][2]  = -1 * ((*it).get_ID());
    output[i][3]  = (*it).death_time;
  }

  return output;
}

void force_output(std::string s) {
  Rcout << s << "\n";
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}


float calc_nltt(const std::vector< float >& emp_brts,
                const std::vector< std::vector< float > >& ltab) {

  std::vector< float > b1(emp_brts.size());
  std::vector< float > b2(ltab.size());
  for(int i = 0; i < emp_brts.size(); ++i) {
    b1[i] = -1 * emp_brts[i];
  }
  for(int i = 0; i < ltab.size(); ++i) {
    b2[i] = -1 * ltab[i][0];
  }
  b1.push_back(0.f);
  b2.push_back(0.f);
  std::sort(b1.begin(), b1.end(), std::greater<float>());
  std::sort(b2.begin(), b2.end(), std::greater<float>());
  return 0.f;


  /*
   * b1 <- -1 * c(rev(sort(unique(ltab1[, 1]))), 0)
   b2 <- -1 * c(rev(sort(unique(ltab2[, 1]))), 0)

   num_lin1 <- 2:length(b1)
   num_lin2 <- 2:length(b2)
   b1 <- 1 - b1 / min(b1)
   b2 <- 1 - b2 / min(b2)
   num_lin1 <- num_lin1 / max(num_lin1)
   num_lin2 <- num_lin2 / max(num_lin2)

   all_b_times <- unique(sort(c(b1, b2)))
   diff <- 0
   for (k in 2:length(all_b_times)) {
   tim <- all_b_times[k]
   index1 <- max(which(b1 < tim))
   index2 <- max(which(b2 < tim))
   lins1 <- num_lin1[max(index1, 1)]
   lins2 <- num_lin2[max(index2, 1)]
   dt <- all_b_times[k] - all_b_times[k - 1]
   diff <- diff + dt * abs(lins1 - lins2)
   }
   return(diff)
   *
   *
   *
   */



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
  int model = 1 + rndgen.random_number(3);
  std::vector<float> output = parameters_from_prior(rndgen, model);

  return(output);
}

//' draw parameter combinations from the prior
//' @param model chosen model
//' @return vector with 6 entries: extinction, sympatric speciation at high
//' water, sympatric speciation at low water, allopatric speciation, amount of
//' perturbation, the chosen water model
//' @export
// [[Rcpp::export]]
std::vector<float> param_from_prior_model_cpp(int model) {

  rnd_t rndgen;
  std::vector<float> output = parameters_from_prior(rndgen, model);

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
  std::vector<float> output = parameters_from_prior(rndgen);

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

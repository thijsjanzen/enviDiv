#include "sim2.h"


 //' simulate a tree using environmental diversification
 //' @param model chosen model
 //' @param parameters a vector of parameters in order: [extinction,
 //' sym_spec_high, sym_spec_low, allo_spec, perturbance, water_rate, model]
 //' @param crown_age crown age of the tree to be simulated
 //' @param max_lin maximum number of lineages
 //' @return RcppList
 //' @export
 // [[Rcpp::export]]
 Rcpp::List sim_envidiv2_cpp(std::vector<double> parameters,
                             int model,
                             double crown_age,
                             int max_lin,
                             int seed) {


   new_sim::simulation sim(parameters, model, crown_age, max_lin, seed);
   sim.run();


  return Rcpp::List::create( Rcpp::Named("code") = sim.run_info,
                                Rcpp::Named("Ltable") = sim.get_ltable(),
                                Rcpp::Named("num_spec") = sim.crowns[0] + sim.crowns[1],
                                Rcpp::Named("water") = sim.waterlevels);
 }



//' simulate a tree using environmental diversification
 //' @param model chosen model
 //' @param parameters a vector of parameters in order: [extinction,
 //' sym_spec_high, sym_spec_low, allo_spec, perturbance, water_rate, model]
 //' @param crown_age crown age of the tree to be simulated
 //' @param max_lin maximum number of lineages
 //' @return RcppList
 //' @export
 // [[Rcpp::export]]
 Rcpp::List sim_new_cond_cpp(int model,
                             double crown_age,
                             int min_lin,
                             int max_lin,
                             int num_tries) {

    new_sim::simulation sim(model, crown_age, max_lin);
    bool finished = false;
    for (int i = 0; i < num_tries; ++i) {
      sim.params_ = parameters_from_prior(sim.rnd, model);
      sim.run();
      if (sim.run_info == "done") {
         auto num_lin = sim.get_num_lin();
         if (num_lin >= min_lin && num_lin <= max_lin) {
            finished = true;
            break;
         }
      }
    }
    if (!finished) {
      sim.run_info = "failure_to_simulate";
    }

    return Rcpp::List::create( Rcpp::Named("code") = sim.run_info,
                               Rcpp::Named("Ltable") = sim.get_ltable(),
                               Rcpp::Named("num_spec") = sim.get_num_lin(),
                               Rcpp::Named("water") = sim.waterlevels,
                               Rcpp::Named("parameters") = sim.params_);
 }

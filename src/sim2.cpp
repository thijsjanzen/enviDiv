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

   std::array<double, 6> params = {parameters[0], parameters[1], parameters[2], parameters[3],
                                   parameters[4], parameters[5]};

   new_sim::simulation sim(params, model, crown_age, max_lin, seed);
   sim.run();


  return Rcpp::List::create( Rcpp::Named("code") = sim.run_info,
                                Rcpp::Named("Ltable") = sim.get_ltable(),
                                Rcpp::Named("num_spec") = sim.crowns[0] + sim.crowns[1],
                                Rcpp::Named("water") = sim.waterlevels);
 }

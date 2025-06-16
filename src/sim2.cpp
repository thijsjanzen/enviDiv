#include "sim2.h"
#include "util.h"

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
                             double crown_age,
                             int max_lin) {

   rnd_t reng;
   std::vector<double> waterlevel_changes = get_waterlevel_changes(model,
                                                                   crown_age,
                                                                   reng,
                                                                   parameters[ param_type::water_rate]);

   std::array<double, 4> params = {parameters[0], parameters[1], parameters[2], parameters[3]};

   new_sim::simulation sim(params, crown_age, waterlevel_changes, max_lin);
   sim.run();

   
  return Rcpp::List::create( Rcpp::Named("code") = sim.run_info,
                                Rcpp::Named("Ltable") = sim.get_ltable(),
                                Rcpp::Named("num_spec") = sim.crowns[0] + sim.crowns[1],
                                Rcpp::Named("water") = waterlevel_changes);
 }

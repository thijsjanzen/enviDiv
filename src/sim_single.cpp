#include <cmath>
#include <string>

#include "random_thijs.h"
#include "util.h"
#include "Gillespie.h"



//' simulate a tree using environmental diversification
//' @param model chosen model
//' @param parameters a vector of parameters in order: [extinction,
//' sym_spec_high, sym_spec_low, allo_spec, perturbance, water_rate, model]
//' @param crown_age crown age of the tree to be simulated
//' @param min_lin minimum number of lineages
//' @param max_lin maximum number of lineages
//' @return newick string
//' @export
// [[Rcpp::export]]
Rcpp::List sim_envidiv_cpp(int model,
                           std::vector<double> parameters,
                           double crown_age,
                           int max_lin) {

    rnd_t reng;
    std::vector<double> waterlevel_changes = get_waterlevel_changes(model,
                                                                   crown_age,
                                                                   reng,
                                                                   parameters[ param_type::water_rate]);

    ltab l_table;

    std::string code = do_run_r(parameters,
                                waterlevel_changes,
                                crown_age,
                                max_lin,
                                l_table,
                                reng);

    if (code == "success") {
       auto newick = ltable_to_newick(l_table, crown_age);
       return Rcpp::List::create( Rcpp::Named("code") = newick,
                                  Rcpp::Named("water") = waterlevel_changes);
    }

    return Rcpp::List::create( Rcpp::Named("code") = "failure",
                               Rcpp::Named("water") = waterlevel_changes);
}

Rcpp::NumericMatrix convert_to_matrix(
        const std::vector< std::vector<double>>& v) {

    int nrow = v.size();
    int ncol = v[0].size();

    Rcpp::NumericMatrix out(nrow, ncol);
    for (int i = 0; i < v.size(); ++i) {
        for (int j = 0; j < v[0].size(); ++j) {
            out(i, j) = v[i][j];
        }
    }
    return out;
}

//' generate particles from the prior
//' @param num_particles number of particles
//' @param crown_age crown age
//' @param min_lin minimum number of lineages
//' @param max_lin maximum number of lineages
//' @param verbose verbose output
//' @return numeric matrix
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix initial_draw_from_prior(int num_particles,
                                            double crown_age,
                                            int min_lin,
                                            int max_lin,
                                            bool verbose) {

    rnd_t reng;
    ltab l_table;
    std::vector< std::vector< double >> accepted_particles;

    int updateFreq = num_particles / 20;
    if(updateFreq < 1) updateFreq = 1;

    if(verbose) {
        Rcpp::Rcout << "0--------25--------50--------75--------100\n";
        Rcpp::Rcout << "*";
    }

    size_t prev_print = accepted_particles.size();

    const int num_models = 4;
    const int num_particles_per_model = num_particles / num_models;

    for (size_t model = 1; model <= num_models; ++model) {
        int num_accepted = 0;

        while(num_accepted < num_particles_per_model) {

            auto parameters = param_from_prior_cpp(model);
            // parameters[param_type::model] = model;

            auto rate = parameters[ param_type::water_rate];

            std::vector<double> waterlevel_changes = get_waterlevel_changes(model,
                                                                            crown_age,
                                                                            reng,
                                                                            rate);
            std::string code = do_run_r(parameters,
                                        waterlevel_changes,
                                        crown_age,
                                        max_lin,
                                        l_table,
                                        reng);
            if (code == "success") {
                int num_species = 0;

                for (size_t i = 0; i < l_table.size(); ++i) {
                    if (l_table[i][3] == -1) num_species++;
                }
                if (num_species >= min_lin) {
                    accepted_particles.push_back(parameters);
                    num_accepted++;
                }
            }

            if (accepted_particles.size() % updateFreq == 0 &&
                verbose &&
                accepted_particles.size() != prev_print) {
                Rcpp::Rcout << "**";
                prev_print = accepted_particles.size();
            }
            Rcpp::checkUserInterrupt();
        }
    }

    Rcpp::NumericMatrix out = convert_to_matrix(accepted_particles);

    return out;
}

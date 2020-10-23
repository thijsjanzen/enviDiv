#include "random_thijs.h"
#include "util.h"
#include "Gillespie.h"
#include "statistics.h"
#include "particle.h"

#include <Rcpp.h>
using namespace Rcpp;

// Function get_accepted_from_R = Environment::global_env()["accept_from_R"];

void update_weights(std::vector< std::vector< particle >>& p,
                    std::vector< float >& max_weights,
                    std::vector< float >& marginal_weights) {
  double s = 0.0;
  for(const auto& i : p) {
    for(const auto& j : i){
      s += j.weight;
    }
  }


  double factor = 1.0 / s;
  max_weights = std::vector<float>(p.size(), -1.0);
  marginal_weights = std::vector<float>(p.size(), 0.0);

  for(int i = 0; i < p.size(); ++i) {
    for(auto& j : p[i]) {
      j.weight *= factor;
      marginal_weights[i] += j.weight;
      if (j.weight > max_weights[i])
        max_weights[i] = j.weight;
    }
  }

  return;
}

std::string do_run(particle& p,
                   float maximum_time,
                   int max_lin,
                   rnd_t& rndgen);

std::vector<particle> convert_matrix(const NumericMatrix& m);
NumericMatrix convert_to_matrix(const std::vector<particle>& v);



//' simulate a tree using environmental diversification
//' @param m1 a matrix with particles of model 1
//' @param m2 a matrix with particles of model 2
//' @param m3 a matrix with particles of model 3
//' @param m_weights a vector with weight of each model
//' @param max_w a vector with the maximum weights of particles in each model
//' @param batch_size number of particles to generate
//' @param crown_age crown age
//' @param num_lin number of lineages to condition on
//' @param num_threads number of threads to use
//' @param sd_p standard deviation of parameter perturbation
//' @param self_prob_m probability of staying in the same model
//' @return list with generated particles
//' @export
// [[Rcpp::export]]
List smc_abc_batch(const NumericMatrix& m1,
                   const NumericMatrix& m2,
                   const NumericMatrix& m3,
                   const NumericVector& m_weights,
                   const NumericVector& max_w,
                   int batch_size,
                   float crown_age,
                   int num_lin,
                   int num_threads,
                   double sd_p,
                   double self_prob_m) {

  rnd_t reng;
  std::vector< particle > temp;
  std::vector< std::vector< particle > > pop(3, temp);

  force_output("conversion of input matrices to particles");
  pop[0] = convert_matrix(m1); force_output("m1 done");
  pop[1] = convert_matrix(m2); force_output("m2 done");
  pop[2] = convert_matrix(m3); force_output("m3 done");

  std::vector< std::vector < particle > > new_pop(3, temp);

  Rcout << "calculate max weights\n";
  std::vector<float> model_weights(m_weights.begin(), m_weights.end());

  Rcout << "model weights are:\n";
  for(int i = 0; i < 3; ++i) {
    Rcout << model_weights[i] << " ";
  }
  Rcout << "\n";

  std::vector< float > max_weights(max_w.begin(), max_w.end());
  Rcout << "max weights are:\n";
  for(int i = 0; i < 3; ++i) {
    Rcout << max_weights[i] << " ";
  }
  Rcout << "\n";

  Rcout << "first particle for each model is:\n";
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < pop[i][0].parameters.size(); ++j) {
      Rcout << pop[i][0].parameters[j] << " ";
    }
    Rcout << "\n";
  }



  int num_accepted = 0;
  Rcout << "starting generation of new particles\n";
  while(num_accepted < batch_size) {

    int model = sample_model(model_weights, reng);
    Rcout << "model drawn: " << model << "\n";
    int param_index = sample_param(pop[model], max_weights[model], reng);
    Rcout << "param_index drawn: " << param_index << "\n";

    particle new_particle = pop[model][param_index];
    new_particle.perturb(reng);

    for(int i = 0; i < new_particle.parameters.size(); ++i) {
      Rcout << new_particle.parameters[i] << " ";
    }
    Rcout << new_particle.model << "\n";

    std::string code = do_run(new_particle,
                              crown_age,
                              num_lin,
                              reng);

    new_particle.newick_string = ltable_to_newick(new_particle.l_table,
                                                  crown_age);

    new_pop[new_particle.model].push_back(new_particle);
    Rcout << num_accepted << "\n";
  }

  Rcout << "done, converting back to R format\n";
  NumericMatrix o1, o2, o3;
  o1 = convert_to_matrix(new_pop[0]);
  o2 = convert_to_matrix(new_pop[1]);
  o3 = convert_to_matrix(new_pop[2]);

  return List::create( Named("m1") = o1,
                       Named("m2") = o2,
                       Named("m3") = o3);
}



// [[Rcpp::export]]
int smc_abc_par(int num_particles,
                float crown_age,
                int num_lin,
                int num_threads,
                int batch_size,
                double sd_p,
                double self_prob_m,
                const std::vector<float>& emp_stats,
                const std::vector<float>& thresholds,
                const std::vector<float>& emp_brts) {





  // we first do serial, to be safe

  rnd_t reng;

  std::vector< particle > temp;
  std::vector< std::vector < particle > > pop(3, temp);
  std::vector< std::vector < particle > > new_pop(3, temp);


  for(int i = 0; i < num_particles; ++i) {
    particle new_particle(self_prob_m, sd_p);
    new_particle.random_parameters(reng, crown_age);
    pop[ new_particle.model ].push_back(new_particle);
  }
  std::vector< float > model_weights(3);
  std::vector< float > max_weights(3, 1.0);

  for(int i = 0; i < pop.size(); ++i) {
    model_weights[i] = pop[i].size() * 1.0 / num_particles;
  }

  for(int iter = 1; iter < thresholds.size(); ++iter) {
    // sample particle
    int num_accepted = 0;
    while(num_accepted < num_particles) {


      std::vector< particle > new_batch(batch_size);
      std::vector< std::string > newick_batch(batch_size);

      for(int r = 0; r < batch_size; ++r) {

        int model = sample_model(model_weights, reng);
        int param_index = sample_param(pop[model],
                                       max_weights[model],
                                       reng);

        particle new_particle = pop[model][param_index];


        std::string code = do_run(new_particle,
                                  crown_age,
                                  num_lin,
                                  reng);

        new_particle.newick_string = ltable_to_newick(new_particle.l_table,
                                                      crown_age);

        //   new_particle.update_sum_stats(emp_brts);

        new_batch[r] = new_particle;
        newick_batch[r] = new_particle.newick_string;

      }

      StringVector pass_newick;
      for(int i = 0; i < newick_batch.size(); ++i) {
        pass_newick[i] = newick_batch[i];
      }
      NumericVector accepted; // = get_accepted_from_R(emp_stats,
      //                     pass_newick,
      //                     thresholds[iter],
      //                     emp_brts);

      for(int i = 0; i < accepted.size(); ++i) {
        if (accepted[i]) {
          new_batch[i].calc_weight(pop[new_batch[i].model], model_weights);
          new_pop[new_batch[i].model].push_back(new_batch[i]);
          num_accepted++;
        }
      }
    }

    // normalize and calculate marginal weights
    update_weights(pop, max_weights, model_weights);

    std::swap(pop, new_pop);
  }

  return 0;
}

int count_extant(const std::vector<species>& s) {
  int num_extant = 0;
  for(const auto i : s) {
    if (i.death_time == -1) num_extant++;
  }
  return num_extant;
}

std::string do_run(particle& p,
                   float maximum_time,
                   int num_lin,
                   rnd_t& rndgen)
{
  std::vector< species > s1;

  float jiggle_amount = p.parameters[4];
  rndgen.set_normal_trunc(0.0f, jiggle_amount);

  //int idCount = 0;
  std::vector < std::vector < float > > l_table1;
  int error_code = run(p.parameters, p.waterlevel_changes,
                       s1, maximum_time, num_lin,
                       rndgen);

  if (error_code == 0) {
    return "extinction";
  }

  if (error_code == 1e6) {
    return "overflow";
  }

  int extant_species = count_extant(s1);
  num_lin -= extant_species;

  std::vector<species> s2;

  int error_code2 = run(p.parameters,
                        p.waterlevel_changes,
                        s2, maximum_time, num_lin,
                        rndgen);

  if (error_code2 == 0) {
    return "extinction";
  }
  if (error_code2 == 1e6) {
    return "overflow";
  }

  jiggle(s1, s2, maximum_time, jiggle_amount, rndgen);

  p.l_table = create_l_table_float(s1, s2);

  std::string output = "success";
  return output;
}


std::vector<particle> convert_matrix(const NumericMatrix& m) {
  std::vector< particle > output(m.nrow());
  force_output("start conversion");
  Rcout << m.nrow() << "\t" << m.ncol() << "\n";
  for(int i = 0; i < m.nrow(); ++i) {
    particle new_particle;
    NumericVector row_entry = m(i, _);

    for(int j = 0; j < row_entry.size(); ++j) {
      new_particle.parameters.push_back(row_entry[j]);
    }
    new_particle.weight = new_particle.parameters.back();
    new_particle.parameters.pop_back();
    //Rcout << m.nrow() << " " << m.ncol() << " " << i << " " << new_particle << "\t"; force_output("");
    output[i] = new_particle;
  }
  force_output("matrix conversion done\n");
  return(output);
}

NumericMatrix convert_to_matrix(const std::vector<particle>& v) {
  NumericMatrix output_matrix(v.size(), v[0].parameters.size() + 1);
  for(int i = 0; i < v.size(); ++i) {
    for(int j = 0; j < v[i].parameters.size(); ++j) {
      output_matrix(i, j) = v[i].parameters[j];
    }
    output_matrix(i, v[i].parameters.size() + 1) = v[i].weight;
  }
  return output_matrix;
}

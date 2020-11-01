#include <unistd.h>

#include "random_thijs.h"
#include "util.h"
#include "Gillespie.h"
#include "statistics.h"
#include "particle.h"
#include "parallel.h"

#include <tbb/tbb.h>

#include <Rcpp.h>
using namespace Rcpp;

// Function get_accepted_from_R = Environment::global_env()["accept_from_R"];

std::string do_run(const std::vector<float>& parameters,
                   const std::vector<float>& waterlevel_changes,
                   float maximum_time,
                   int max_lin,
                   std::vector< std::vector< float > >& l_table,
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
//' @param min_lin minimum number of lineages to condition on
//' @param max_lin maximum number of lineages to condition on
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
                   int min_lin,
                   int max_lin,
                   int num_threads,
                   double sd_p,
                   double self_prob_m) {

  std::vector< particle > temp;
  std::vector< std::vector< particle > > pop(3, temp);

  pop[0] = convert_matrix(m1);
  pop[1] = convert_matrix(m2);
  pop[2] = convert_matrix(m3);

  std::vector< std::vector < particle > > new_pop(3, temp);

  std::vector<float> model_weights(m_weights.begin(), m_weights.end());

  std::vector< float > max_weights(max_w.begin(), max_w.end());

  int num_accepted = 0;
  int num_tried = 0;
  float accept_rate = 1.0;

  while(num_accepted < batch_size) {

    int loop_size = batch_size - num_accepted;
    if (loop_size < 10) loop_size = 10;
    if (accept_rate > 0) loop_size *= 1.0 / accept_rate;
    if (loop_size > 1e3) loop_size = 1e3;

    std::vector< particle > add(loop_size);
    std::vector< bool > add_flag(loop_size, false);

    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, loop_size - 1),
      [&](const tbb::blocked_range<unsigned>& r) {

        rnd_t reng;
        for (unsigned i = r.begin(); i < r.end(); ++i) {
          int model = sample_model(model_weights, reng);

          int param_index = sample_param(pop[model], max_weights[model], reng);

          particle new_particle = pop[model][param_index];
          new_particle.perturb(reng);
          new_particle.update_waterlevel_changes(reng, crown_age);

          std::vector< std::vector< float > > temp_ltable;

          std::string code = do_run(new_particle.parameters,
                                    new_particle.waterlevel_changes,
                                    crown_age,
                                    max_lin,
                                    temp_ltable,
                                    reng);

          if (code == "success") {
            new_particle.l_table = temp_ltable;
            int num_lin = new_particle.count_num_lin();
            if (num_lin >= min_lin && num_lin <= max_lin) {
              new_particle.newick_string = ltable_to_newick(new_particle.l_table,
                                                            crown_age);
              add[i] = new_particle;
              add_flag[i] = true;
            }
          }
        }
      });

    for (int k = 0; k < add_flag.size(); ++k) {
      if (add_flag[k] == true) {
        new_pop[ add[k].model ].push_back(add[k]);
        num_accepted++;
      }
    }

    num_tried    += loop_size;
    accept_rate = 1.0f * num_accepted / num_tried;
    Rcout << num_accepted << " " << loop_size << " " <<
      num_tried << " " << accept_rate << "\n";
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
        new_particle.update_waterlevel_changes(reng, crown_age);


        std::string code = do_run(new_particle.parameters,
                                  new_particle.waterlevel_changes,
                                  crown_age,
                                  num_lin,
                                  new_particle.l_table,
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

std::string do_run(const std::vector<float>& parameters,
                   const std::vector<float>& waterlevel_changes,
                   float maximum_time,
                   int max_lin,
                   std::vector< std::vector< float > >& l_table,
                   rnd_t& rndgen)
{
  std::vector< species > s1;

  float jiggle_amount = parameters[4];
  rndgen.set_normal_trunc(0.0f, jiggle_amount);

  // force_output("run_left");
  int error_code = run(parameters, waterlevel_changes,
                       s1, maximum_time, max_lin,
                       rndgen);

  if (error_code == 0) {
    return "extinction";
  }

  if (error_code == 1e6) {
    return "overflow";
  }

  std::vector<species> s2;
  // force_output("run_right");
  int error_code2 = run(parameters,
                        waterlevel_changes,
                        s2, maximum_time, max_lin,
                        rndgen);

  if (error_code2 == 0) {
    return "extinction";
  }
  if (error_code2 == 1e6) {
    return "overflow";
  }

  jiggle(s1, s2, maximum_time, jiggle_amount, rndgen);

  //Rcout << s1.size() << " " << s2.size() << "\n";
  //force_output("creating ltable");
  //make_sleep(1);
  l_table = create_l_table_float(s1, s2);
  //force_output("ltable done");
  //make_sleep(1);

  std::string output = "success";
  return output;
}



std::vector<particle> convert_matrix(const NumericMatrix& m) {
  std::vector< particle > output(m.nrow());
  //  force_output("start conversion");
  // Rcout << m.nrow() << "\t" << m.ncol() << "\n";
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


std::vector<std::vector< float > > convert_matrix2(const NumericMatrix& m) {
  std::vector< std::vector< float > > output(m.nrow());
  for(int i = 0; i < m.nrow(); ++i) {

    NumericVector row_entry = m(i, _);
    std::vector< float > add(row_entry.begin(), row_entry.end());
    output[i] = add;
  }
  // force_output("matrix conversion done\n");
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

int sample_param_float(const std::vector< std::vector< float >>& p,
                       double max_w,
                       rnd_t& reng) {
  while(true) {
    int index =  reng.random_number(static_cast<int>(p.size()));
    double prob = p[index][5] * 1.0 / max_w;
    if (reng.uniform() < prob) {
      return index;
    }
  }
}

template <typename PK>
void perturb_param(std::vector< float >& parameters,
                   rnd_t& reng,
                   PK param_change) {

  for(int i = 0; i < 5; ++i) {
    parameters[i] = param_change.perturb(parameters[i], reng);
  }
  return;
}




//' simulate a tree using environmental diversification
//' @param model a vector of the four paramaters of the model
//' @param num_repl a vector that indicates the time points of water level changes
//' @param crown_age crown age of the tree to be simulated
//' @param min_lin minimum number of lineages
//' @param max_lin maximum number of lineages
//' @param num_threads number of threads
//' @return newick string
//' @export
// [[Rcpp::export]]
List abc_smc_2(const NumericMatrix& m1,
               const NumericMatrix& m2,
               const NumericMatrix& m3,
               const NumericVector& m_weights,
               const NumericVector& max_w,
               int batch_size,
               float crown_age,
               int min_lin,
               int max_lin,
               int num_threads,
               double sd_p,
               double self_prob_m) {

  statistics::model_dist model_change(self_prob_m);
  statistics::norm_dist_log param_change(sd_p);

  std::vector< std::vector< std::vector< float > > > pop(3);

  pop[0] = convert_matrix2(m1);
  pop[1] = convert_matrix2(m2);
  pop[2] = convert_matrix2(m3);

  std::vector< std::vector< std::vector< float > > > new_pop(3);

  std::vector<float> model_weights(m_weights.begin(), m_weights.end());

  std::vector< float > max_weights(max_w.begin(), max_w.end());

  try {

    std::vector< std::string > trees;
    std::vector< std::vector < float > > parameter_list;

    auto T0 = std::chrono::high_resolution_clock::now();
    int loop_size = batch_size - trees.size();

    float accept_rate = 1.0;
    int num_tried = 0;
    int num_accepted = 0;

    while(trees.size() < batch_size) {
      // loop size can be optimized further, depending on the average success rate
      // e.g. loop_size = loop_size * 1.0f / success_rate
      // this is especially interesting once only a few are left.
      loop_size = batch_size - trees.size();
      Rcout << "trees remaining: " << loop_size   <<
        " accept rate: "    << accept_rate << "\n";

      if (loop_size < 10) loop_size = 10;
      if (accept_rate > 0) loop_size *= 1.0 / accept_rate;
      if (loop_size > 1e4) loop_size = 1e4;


      std::vector< std::string > add(loop_size);
      std::vector< float > temp_filler(5);
      std::vector< std::vector< float > > add_params(loop_size, temp_filler);
      std::vector< bool > add_flag(loop_size, false);

      tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);

      tbb::parallel_for(
        tbb::blocked_range<unsigned>(0, loop_size - 1),
        [&](const tbb::blocked_range<unsigned>& r) {
          rnd_t reng;

          for (unsigned i = r.begin(); i < r.end(); ++i) {

            int model = sample_model_perturb(model_weights, reng, model_change);
            int param_index = sample_param_float(pop[model],
                                                 max_weights[model],
                                                 reng);

            std::vector< float > parameters = pop[model][param_index];
            parameters[5] = model + 1; // model indication is in R notation...
            perturb_param<statistics::norm_dist_log>(parameters, reng, param_change);

            std::vector<float> waterlevel_changes = get_waterlevel_changes(parameters[5],
                                                                           crown_age,
                                                                           reng);

            std::vector< std::vector< float > > l_table;

            std::string code = do_run_tbb(parameters,
                                          waterlevel_changes,
                                          crown_age,
                                          max_lin,
                                          l_table,
                                          reng);

            int num_lin = 0;
            for(int k = 0; k < l_table.size(); ++k) {
              if(l_table[k][3] == -1) num_lin++;
            }

            if(num_lin >= min_lin && num_lin <= max_lin) {
              add[i]        = ltable_to_newick(l_table, crown_age);
              add_params[i] = parameters;
              add_flag[i]   = true;
            }
          }
        });

      for(int j = 0; j < add_flag.size(); ++j) {
        if(add_flag[j]) {
          trees.push_back(add[j]);
          parameter_list.push_back(add_params[j]);
        }
      }

      num_accepted  =  trees.size();
      num_tried    += loop_size;
      accept_rate = 1.0f * num_accepted / num_tried;
    }

    Rcpp::List output(trees.size());
    for(int k = 0; k < trees.size(); ++k) {
      output[k] = trees[k];
    }

    Rcpp::NumericMatrix parameter_matrix(parameter_list.size(), 6);
    for(int k = 0; k < parameter_list.size(); ++k) {
      for(int j = 0; j < parameter_list[0].size(); ++j) {
        parameter_matrix(k, j) = parameter_list[k][j];
      }
    }


    auto T1 = std::chrono::high_resolution_clock::now();
    auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(T1 - T0).count());
    Rcout << "trees simulated in: " << elapsed << "seconds\n";
    return List::create(Named("trees") = trees,
                        Named("parameters") = parameter_matrix);
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    ::Rf_error("c++ exception (unknown reason)");
  }
  return NA_REAL;
}



double prior_pdf(const std::vector< float >& param) {
  double p = 1.0;
  static double range = pow(10, 6);
  for(int i = 0; i < 4; ++i) {
    if (log10(param[i]) < -4) {
      return 0.0;
    }
    p *= 1.0 / range;
  }
  p *= 1.0 / 3.0; // each model has 1/3 chance in the prior.
  return p;
}

template <typename PK>
double calc_prob_other(const std::vector< float > other_particle,
                       const std::vector< float > focal_particle,
                       PK param_change) {

  float p = other_particle[5];
  for(int i = 0; i < 4; ++i) {
    p *= param_change.pdf(other_particle[i], focal_particle[i]);
  }
  return(p);
}


template <typename PK>
double calc_prob_param(const std::vector< std::vector< float > >& other_particles,
                       const std::vector< float >& param,
                       PK param_change) {

  double p = 0.0;
  for(const auto& i : other_particles) {
    p += calc_prob_other<PK>(param, i, param_change);
  }
  return p;
}




template <typename MK, typename PK>
double calc_weight(const std::vector< std::vector< float > >& other_particles,
                   const std::vector< float >& param,
                   float s1,
                   MK model_change,
                   PK param_change) {

  double s2 = calc_prob_param<PK>(other_particles, param, param_change); // s2 = probability of drawing parameter from parameters

  double S = s1 * s2;
  double prior_prob = prior_pdf(param); // prior probability of both model and parameter value.

  return 1.0 * prior_prob / S;
}


//' calculate weights
//' @param m1 matrix with particles from model 1
//' @param m2 matrix with particles from model 2
//' @param m3 matrix with particles from model 3
//' @param p  matrix with particles for which to calculate weight
//' @param m_weights vector with model weights of m1, m2 and m3
//' @param sd_p standard deviation of parameter change
//' @param self_prob_m probability of drawing model
//' @return vector with weights
//' @export
// [[Rcpp::export]]
NumericVector calc_weights(const NumericMatrix& m1,
                           const NumericMatrix& m2,
                           const NumericMatrix& m3,
                           const NumericMatrix& p,
                           const NumericVector& m_weights,
                           double sd_p,
                           double self_prob_m) {

  statistics::model_dist model_change(self_prob_m);
  statistics::norm_dist_log param_change(sd_p);

  std::vector< std::vector< std::vector< float > > > pop(3);

  pop[0] = convert_matrix2(m1);
  pop[1] = convert_matrix2(m2);
  pop[2] = convert_matrix2(m3);

  std::vector<float> model_weights(m_weights.begin(), m_weights.end());

  std::vector< std::vector< float > > particles = convert_matrix2(p);

  NumericVector output(particles.size());

  std::vector< float > model_probs(3, 0.f);
  for(int i = 0; i < 3; ++i) {
    for(int j = 0; j < 3; ++j) {
      model_probs[i] += model_weights[j] * model_change.pdf(i, j);
    }
  }

  for(int j = 0; j < particles.size(); ++j) {
    std::vector< float > focal_particle = particles[j];
    int model = focal_particle[5]; // check index
    float weight = calc_weight<statistics::model_dist,
                               statistics::norm_dist_log>(particles,
                                                          focal_particle,
                                                          model_probs[model],
                                                          model_change,
                                                          param_change);
    output[j] = weight;
  }

  return output;
}

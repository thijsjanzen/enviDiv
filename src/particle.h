//
//  particle.h
//  abc_smc_new
//
//  Created by Thijs Janzen on 21/10/2020.
//  Copyright © 2020 Thijs Janzen. All rights reserved.
//

#ifndef particle_h
#define particle_h

#include <iostream>

#include "statistics.h"
#include "util.h"
#include <Rcpp.h>



struct particle {
  std::vector<float> parameters;
  std::vector<float> waterlevel_changes;

  std::vector< std::vector< float >> l_table;
  std::vector<float> summary_statistics;


  double weight;
  int model;
  std::string newick_string;

  statistics::model_dist model_change;
  statistics::norm_dist_log param_change;

  particle() {
    weight = 1.0;
    model_change = statistics::model_dist(0.5); // default 0.5 self prob
    param_change = statistics::norm_dist_log(0.05); // default sd = 0.05
    summary_statistics = std::vector<float>(10, 0.f);
  }

  particle(double self_prob, double sd) {
    weight = 1.0;
    model_change = statistics::model_dist(self_prob); // default 0.5 self prob
    param_change = statistics::norm_dist_log(sd);        // default sd = 0.05
    summary_statistics = std::vector<float>(10, 0.f);
  }

  particle& operator=(const particle& other) {
    if (this == &other) return *this;

    weight = other.weight;
    model = other.model;
    newick_string = other.newick_string;
    model_change = other.model_change;
    param_change = other.param_change;
    parameters = other.parameters;
    waterlevel_changes = other.waterlevel_changes;
    l_table = other.l_table;
    summary_statistics = other.summary_statistics;
    return *this;
  }

  particle(const particle& other) :
    weight(other.weight),
    model(other.model),
    newick_string(other.newick_string) {
    model_change = other.model_change;
    param_change = other.param_change;
    parameters = other.parameters;
    waterlevel_changes = other.waterlevel_changes;
    l_table = other.l_table;
    summary_statistics = other.summary_statistics;
  }

  friend std::ostream& operator<<(std::ostream& os, const particle& p) {
    for(auto i : p.parameters) {
      os << i << " ";
    }
    os << p.weight << "\n";
    return os;
  }

  void set_parameters(const std::vector< float > & p,
                      float crown_age,
                      rnd_t& reng) {
    parameters= p;
    update_waterlevel_changes(reng, crown_age);
  }

  void update_waterlevel_changes(rnd_t& reng,
                                 float crown_age) {
    model = parameters[5];
    waterlevel_changes = get_waterlevel_changes(model,
                                                crown_age,
                                                reng);
  }

  void random_parameters(rnd_t& reng,
                         float crown_age) {
    parameters = generateprior(reng);
    waterlevel_changes = get_waterlevel_changes(model,
                                                crown_age,
                                                reng);
  }

  std::vector< float > generateprior(rnd_t& rndgen_) {
    std::vector< float > output(6, 0);
    output[0] = powf(10, (-4 + 6 * rndgen_.uniform()));  // extinction
    output[1] = powf(10, (-4 + 6 * rndgen_.uniform()));  // symp spec high
    output[2] = powf(10, (-4 + 6 * rndgen_.uniform()));  // symp spec low
    output[3] = powf(10, (-4 + 6 * rndgen_.uniform()));  // allo spec
    output[4] = powf(10, (-4 + 6 * rndgen_.uniform()));  // jiggle
    output[5] = static_cast<float>(rndgen_.random_number(3));
    return(output);
  }

  double prior_pdf() {
    double p = 1.0;
    static double range = pow(10, 6);
    for(int i = 0; i < 4; ++i) {
      if (log10(parameters[i]) < -4) {
        return 0.0;
      }
      p *= 1.0 / range;
    }
    p *= 1.0 / 3.0; // each model has 1/3 chance in the prior.
    return p;
  }

  void perturb(rnd_t& reng) {
    for(int i = 0; i < 5; ++i) {
      parameters[i] = param_change.perturb(parameters[i], reng);
    }
    parameters[5] = model_change.perturb(parameters[5], reng);
    model = parameters[5];
  }

  double calc_prob_model(const std::vector<float>& m) {
    double p = 0.0;
    for(int i = 0; i < m.size(); ++i) {
      p += m[i] * model_change.pdf(i, model);
    }
    return p;
  }

  double calc_prob_parameter_vector(const std::vector<float>& other_p) {
    double p = 1.0;
    for(int i = 0; i < other_p.size(); ++i) {
      p *= param_change.pdf(parameters[i], other_p[i]);
    }
    return p;
  }

  double calc_prob_param(const std::vector< particle >& other_particles) {
    double p = 0.0;
    for(const auto& i : other_particles) {
      p += i.weight * calc_prob_parameter_vector(i.parameters);
    }
    return p;
  }

  void calc_weight(const std::vector< particle >& other_particles,
                   const std::vector< float >& model_probs) {

    double s1 = calc_prob_model(model_probs); // s1 = probability of drawing model from models
    double s2 = calc_prob_param(other_particles); // s2 = probability of drawing parameter from parameters

    double S = s1 * s2;
    double prior_prob = prior_pdf(); // prior probability of both model and parameter value.

    weight = prior_prob / S;
    return;
  }

  bool calc_accept(const std::vector<float>& emp_stats,
                   float threshold) {
    return true;
  }

  void update_sum_stats(const std::vector<float>& emp_brts) {
   //    summary_statistics[0] = calc_gamma(l_table);
   // summary_statistics[6] = calc_nltt(l_table, emp_brts);
  }
};

int sample_model(const std::vector<float>& w,
                 rnd_t& reng) {
  // we know there are three weights summing to 1
  double r = reng.uniform();
  for(int i = 0; i < w.size(); ++i) {
    r -= w[i];
    if (r <= 0.0) return i;
  }
  return -1 + w.size();
}

int sample_param(const std::vector< particle >& p,
                 double max_p,
                 rnd_t& reng) {
  while(true) {
    int index =  reng.random_number(static_cast<int>(p.size()));
    double prob = 1.0 * p[index].weight / max_p;
  //  Rcpp::Rcout << index << " " << prob << " " << max_p << " " << p[index].weight << "\n";
    if (reng.uniform() < prob) {
      return index;
    }
  }
}


#endif /* particle_h */

#ifndef statistics_h
#define statistics_h

#include <random>
#include <cmath>
#include "random_thijs.h"

constexpr double pi_ = 3.14159265358979323846;

namespace statistics {

  struct norm_dist {
    double m;
    double s;
    double prefactor;

    std::normal_distribution<double> rdist;

    norm_dist() {
      rdist = std::normal_distribution<double>(0, 1);
      prefactor = 1.f / (s * sqrtf(2 * pi_));
    }

    norm_dist(double sd) : m(0), s(sd) {
      rdist = std::normal_distribution<double>(m, s);
      prefactor = 1.f / (s * sqrtf(2 * pi_));
    }

    double pdf(double a, double b) {
      double a1 = (a - b) / s;
      return prefactor * std::exp(-0.5 * a1 * a1);
    }

    double perturb(double a, rnd_t& rndgen) {
      return a + rdist(rndgen.rndgen);
    }
  };

  struct norm_dist_log {
    double m;
    double s;
    double prefactor;

    std::normal_distribution<double> rdist;

    norm_dist_log() {
      rdist = std::normal_distribution<double>(0, 0.05);
      prefactor = 1.f / (s * sqrtf(2 * pi_));
    }

    norm_dist_log(double sd) : m(0), s(sd) {
      rdist = std::normal_distribution<double>(m, s);
      prefactor = 1.f / (s * sqrtf(2 * pi_));
    }

    double pdf(double a, double b) {
      double a1 = (log10(a) - log10(b)) / s;
      return prefactor * std::exp(-0.5 * a1 * a1);
    }

    double perturb(double a, rnd_t& rndgen) {
      double b = log10(a) + rdist(rndgen.rndgen);
      return( pow(10, b));
    }
  };

  struct model_dist {
    double prob_self;
    std::bernoulli_distribution d;

    model_dist() {
      d = std::bernoulli_distribution(0.5);
    }

    model_dist(double s) : prob_self(s) {
      d = std::bernoulli_distribution(0.5);
    };

    double pdf(int a, int b) {
      return (a == b) ? prob_self : 0.5f * (1 - prob_self);
    }

    double perturb(int a, rnd_t& rndgen) {
      if (rndgen.uniform() < prob_self) {
        return a;
      } else {

        if (d(rndgen.rndgen)) {
          int b = a + 1;
          if (b > 3) b -= 3;
          return b;
        } else {
          int b = a - 1;
          if (b < 1) b+= 3;
          return b;
        }
      }
    }
  };
}

#endif /* statistics_h */


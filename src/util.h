#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <array>
#include <thread>
#include <chrono>
#include <cmath>
#include <random>
#include <numeric>
#include <algorithm>

#include "random_thijs.h"
#include "Rcpp.h"

enum param_type {extinction_rate, sym_high_rate, sym_low_rate, allo_rate, wobble_rate, water_rate, model};

std::vector<double> param_from_prior_cpp(int model);
std::vector<double> param_from_prior_exp_cpp(int model);

std::vector<double> get_waterlevel_cpp(int water_model,
                                      double maximum_time);

std::vector<double> parameters_from_prior(rnd_t& rndgen_);
std::vector<double> parameters_from_prior(rnd_t& rndgen_,
                                         int model);
std::vector<double> get_waterlevel_changes(int water_model,
                                          double maximum_time,
                                          rnd_t& rndgen_,
                                          double rate);



struct ltable_entry {
  double bt;
  double parent;
  double daughter;
  double extant;
  double tend;
  std::string label;

  ltable_entry(double b, double p, double d, double e) : bt(b),
        parent(p), daughter(d), extant(e) {
  }

  ltable_entry() {
  }
};



std::string ltable_to_newick(const std::vector< std::array< double, 4 > >& ltable,
                             double crown_age);

// returns low-entropy 512 bit array for seed sequence
// based on std::chrono::high_resolution_clock.
// ripped from rndutils
inline auto make_low_entropy_seed_array() noexcept->std::array<uint64_t, 8>
{
  // the classic: time, advertised with nano-second resolution.
  const auto e1 = static_cast< uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  // different between invocations from different threads within one app: thread-id
  const auto tid = std::this_thread::get_id();
  const uint64_t e2{ std::hash<typename std::remove_const<decltype(tid)>::type>()(tid) };
  return std::array<uint64_t, 8>{ {
    e1, e2,
    0x000000003c10b019, 0x2bf820b4dd7c1a8a,
    0x9901cf90a40883da, 0x5a3686b2e1de6e51,
    0x000000cc0494d228, 0x000000cc04b66740
    }};
}


// random number generator from low-entropy seed sequence
// ripped from rndutils
template <typename URNG>
inline auto make_random_engine() -> URNG
{
  auto seed_array = make_low_entropy_seed_array();
  std::seed_seq sseq(seed_array.cbegin(), seed_array.cend());
  return URNG(sseq);
}

#endif

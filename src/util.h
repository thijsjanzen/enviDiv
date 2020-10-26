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
#include "Gillespie.h"

#include <Rcpp.h>
using namespace Rcpp;

void make_sleep(size_t ms);

std::vector<float> param_from_prior_cpp();
std::vector<float> param_from_prior_exp_cpp();
std::vector<float> get_waterlevel_cpp(int water_model,
                                      float maximum_time);

// std::vector<float> parameters_from_prior(rnd_t& rndgen_);
std::vector<float> parameters_from_prior(rnd_t& rndgen_,
                                         int model);
std::vector<float> get_waterlevel_changes(int water_model,
                                          float maximum_time,
                                          rnd_t& rndgen_);



struct ltable_entry {
  float bt;
  float parent;
  float daughter;
  float extant;
  float tend;
  std::string label;

  ltable_entry(float b, float p, float d, float e) : bt(b),
        parent(p), daughter(d), extant(e) {
  }

  ltable_entry() {
  }
};

std::string ltable_to_newick(const std::vector< std::vector< float > >& ltable,
                             float crown_age);

std::vector< std::vector< float >> create_l_table_float(
    const std::vector< species > & s1,
    const std::vector< species > & s2);

float calc_nltt(const std::vector< float >& b1,
                const std::vector< std::vector< float > >& ltab);

void force_output(std::string s);

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

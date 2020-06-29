#ifndef RANDOM_THIJS
#define RANDOM_THIJS
#include <random>

struct rnd_t {
  std::mt19937 rndgen;

  rnd_t() {
    std::random_device rd;
    std::mt19937 rndgen_t(rd());
    rndgen = rndgen_t;
  }

  rnd_t(unsigned seed) {
    rndgen = std::mt19937(seed);
  }

  rnd_t(std::mt19937 r) {
    rndgen = r;
  }


  std::uniform_real_distribution<float> unif_dist =
    std::uniform_real_distribution<float>(0.0f, 1.0f);

  std::normal_distribution<float> norm_dist_trunc =
    std::normal_distribution<float>(0.0f, 0.1f);

  int random_number(int n)    {
    if(n <= 1) return 0;
    return std::uniform_int_distribution<> (0, n - 1)(rndgen);
  }

  float uniform()    {
    return unif_dist(rndgen);
  }

  float Expon(float lambda) {
    if (lambda == 0.0) return 1e20f;
    return std::exponential_distribution<float>(lambda)(rndgen);
  }

  float normal(float m, float s) {
    return std::normal_distribution<float>(m, s)(rndgen);
  }

  void set_seed(unsigned seed)    {
    std::mt19937 new_randomizer(seed);
    rndgen = new_randomizer;
  }

  void set_normal_trunc(float m, float s) {
    norm_dist_trunc = std::normal_distribution<float>(m, s);
  }

  float trunc_norm_predif(float trunc) {
    float output = norm_dist_trunc(rndgen);
    while(output > trunc || output < -trunc) {
      output = norm_dist_trunc(rndgen);
    }
    return output;
  }

  float trunc_normal(float m, float s, float t) {
    float output = std::normal_distribution<float>(m, s)(rndgen);
    float diff = output - m;
    if(diff < 0) diff *= -1.f;
    while(diff > t) {
      output = std::normal_distribution<float>(m, s)(rndgen);
      diff = output - m;
      if(diff < 0) diff *= -1.f;
    }
    return output;
  }

};


#endif

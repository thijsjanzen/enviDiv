#include "random_thijs.h"


std::random_device rd;
std::mt19937 rndgen(rd());  //< The one and only random number generator

std::uniform_real_distribution<> unif = std::uniform_real_distribution<>(0.0, 1.0);

int random_number(int n)    {
	return std::uniform_int_distribution<> (0, n-1)(rndgen);
}

float uniform()    {
	return unif(rndgen);
}

float Expon(float lambda) {
    if(lambda == 0.0) return 1e20;
    return std::exponential_distribution<float>(lambda)(rndgen);
}

float normal(float m, float s)   {
  return std::normal_distribution<>(m, s)(rndgen);
}

float trunc_normal(float m, float s, float t) {
  float output = normal(m, s);
  while(output >= t || output <= -t) {
    output = normal(m, s);
  }
  return output;
}


void set_seed(unsigned seed)    {
    std::mt19937 new_randomizer(seed);
    rndgen = new_randomizer;
}

#include "random_thijs.h"


std::random_device rd;
std::mt19937 rndgen(rd());  //< The one and only random number generator

int random_number(int n)    {
	return std::uniform_int_distribution<> (0, n-1)(rndgen);
}

double uniform()    {
	return std::uniform_real_distribution<>(0, 1.0)(rndgen);
}

double normal(double m, double s)   {
    return std::normal_distribution<>(m, s)(rndgen);
}

double Gamma(int k, double lambda)  {
    return std::gamma_distribution<double>(k, lambda)(rndgen);
}

double Expon(double lambda) {
    return std::exponential_distribution<double>(lambda)(rndgen);
}

double trunc_normal(double m, double s, double t) {
  double output = normal(m, s);
  while(output >= t || output <= -t) {
    output = normal(m, s);
  }
  return output;
}


void set_seed(unsigned seed)    {
    std::mt19937 new_randomizer(seed);
    rndgen = new_randomizer;
}

#ifndef RANDOM_THIJS
#define RANDOM_THIJS

#include <random>

int random_number(int n);
double uniform();
double normal(double m, double s);
double Gamma(int k, double lambda);
double Expon(double lambda);

void set_seed(unsigned seed);

#endif 

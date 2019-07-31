#ifndef RANDOM_THIJS
#define RANDOM_THIJS

#include <random>

int random_number(int n);
double uniform();
double Expon(double lambda);
double normal(double m, double s);
double trunc_normal(double m, double s, double t);

void set_seed(unsigned seed);

#endif

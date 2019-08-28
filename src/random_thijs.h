#ifndef RANDOM_THIJS
#define RANDOM_THIJS

#include <random>

int random_number(int n);
float uniform();
float Expon(float lambda);
float normal(float m, float s);
float trunc_normal(float m, float s, float t);

void set_seed(unsigned seed);

#endif

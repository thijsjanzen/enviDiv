#ifndef GILLESPIE_H
#define GILLESPIE_H

#include <utility>    // std::move, std::forward
#include <vector>
#include <string>

#include "random_thijs.h"
#include "species.h"
#include "spec_point.h"
#include "allo_pair.h"

std::string do_run(std::vector<double> parameters,
                   std::vector<double> waterlevel_changes,
                   int maximum_time);

int run(const std::vector<double> parameters,
        const std::vector<double>& W,
        int& id_count,
        std::vector<species>& allSpecies,
        int maximum_time);

int drawEvent(double E, double S, double A);

void updatePairs(std::vector<species>& v,
                  std::vector<allo_pair>& p);

void waterLevelChange(std::vector<species>& v, int& wLevel);

template<typename T>
void removeDuplicates(std::vector<T>& vec);

void waterLevelChange(std::vector<species>& v,
                      int& wLevel);

void extinction(std::vector<species>& v,
                std::vector<species>& extinct_species,
                double time,
                int wLevel);

void Symp_speciation(std::vector<species>& v,
                     int& id_count,
                     std::vector<species>& extinct_species,
                     double time,
                     double waterTime,
                     std::vector<double>& specTimes,
                     int wLevel);

void Allo_speciation(std::vector<species>& v,
                     int& id_count,
                     double time,
                     double water_time,
                     const std::vector<allo_pair>& p,
                     std::vector<double>& specTimes,
                     std::vector<species>& extinct_species);

std::string create_newick_string_r(const std::vector<species>& s1,
                                   const std::vector<species>& s2,
                                   double maximum_time);


void jiggle(std::vector< species > & s1,
            std::vector< species > & s2,
            double maximum_time,
            double jiggle_amount);

void remove_extinct_branches(std::vector<species>& allSpecies);
void merge_single_branches(std::vector<species>& all_species);
#endif

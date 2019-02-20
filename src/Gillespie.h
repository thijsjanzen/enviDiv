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

std::vector<spec_point> run(const std::vector<double> parameters,
                            const std::vector<double>& W,
                            int& id_count,
                            std::vector<species>& allSpecies,
                            int& lins,
                            int maximum_time);

int drawEvent(double E, double S, double A);

void updatePairs2(std::vector<species>& v,
                  std::vector<allo_pair>& p);

void waterLevelChange(std::vector<species>& v, int& wLevel);


std::string create_newick_string( const std::vector<species>& s1,
                                  const std::vector<species>& s2,
                                  const std::vector<spec_point>& b1,
                                  const std::vector<spec_point>& b2);

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
                     std::vector<double>& specTimes);

std::vector<spec_point>  calculateLineages_noextinct(const std::vector<species>& allSp);
std::vector<spec_point>  calculateLineages_withextinct(    std::vector<species>& allSp);





#endif

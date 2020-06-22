#ifndef GILLESPIE_H
#define GILLESPIE_H

#include <vector>
#include <istream>
#include <iostream>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include "Rcpp.h"
#include "random_thijs.h"

struct species
{
  species()
  {
    ID = -1;
    death_time = -1;
  }
  species(int& id);
  species(const species& parent_species,
          int& id_count,
          float b_time);

  species(const species& other):
    birth_time(other.birth_time),
    death_time(other.death_time),
    extant_offspring(other.extant_offspring),
    checked(other.checked),
    ID(other.ID),
    parent(other.parent)
    {	}

  float birth_time; //time of birth;
  float death_time; //time of death;

  bool extant_offspring;
  bool checked;

  void updateHistory(float t);
  bool check_has_viable_offspring(std::vector<species>& v);

  bool operator< (const species& other) const {return       ID < other.ID; }
  bool operator> (const species& other) const {return       ID > other.ID;}
  bool operator<=(const species& other) const {return !operator> (other);}
  bool operator>=(const species& other) const {return !operator< (other);}

  bool operator==(const species& other) const { return ID == other.ID; }
  bool operator!=(const species& other) const {return !((*this) == other);}

  species& operator=(const species& other);

  void set_ID(int id_number) {
    ID = id_number;
  }
  int get_ID() const {
    return(ID);
  }
  void set_parent(int parent_number) {
    parent = parent_number;
  }
  int get_parent() const {
    return parent;
  }

private:
  int ID;
  int parent;
};

struct sort_by_birthtime
{
  inline bool operator() (const species& struct1, const species& struct2)
  {
    return (struct1.birth_time < struct2.birth_time);
  }
};

struct allo_pair
{
  int ID;
  int index_a;
  int index_b;

  allo_pair(int id, int index_A) : ID(id), index_a(index_A) {
    index_b = -1;
  }

  allo_pair(int id, int index_A, int index_B) : ID(id),
                                                index_a(index_A),
                                                index_b(index_B)
  {
  }

};

int drawEvent(float E, float W, float S, float A);

void extinction(std::vector<species>& v,
                std::vector<species>& extinct_species,
                float time, int wLevel);
void extinction2(std::vector<species>& v,
                 std::vector<species>& extinct_species,
                 float time, int waterLevel);

void waterLevelChange(std::vector<species>& v, int& wLevel);
void symp(std::vector<species>& v, int i, std::vector<species>& e,
          int& id_count, float time);

void Symp_speciation(std::vector<species>& v,
                     int& id_count,
                     std::vector<species>& extinct_species,
                     float time,
                     float waterTime,
                     int wLevel);

void Allo_speciation(std::vector<species>& v,
                     int& id_count,
                     float time,
                     float water_time,
                     const std::vector<allo_pair>& p,
                     std::vector<species>& extinct_species);

void updateEpsilonVectors();
std::string generateFileName(int i);
void progressBar(float percent);

template<typename T>
void removeDuplicates(std::vector<T>& vec);

float get_min_time();

std::vector<float> writeMean(const std::vector< std::vector< int > >& outcomes);
void writeAll(const std::vector< std::vector< int > >& outcomes);

void writeWater(const std::vector<std::vector< int> >& W);

template<typename Iterator>
void bubbleSort(Iterator first, Iterator last);

template<typename T>
float calcMean(const std::vector<T> v);

bool sortOnTime(const species& left, const species& right);


std::vector<float> generateWaterLevelChanges(int model);
int find_parent( int ID, const std::vector<species>& v);

int find_Allo(int i, const std::vector<species>& v);
std::vector<int> findOffspring(int ID, const std::vector<species>& v);

template <typename T>
void sortit(std::vector<T>& v);

bool file_exists(const std::string& name);

std::string create_newick_string_r(const std::vector<species>& s1,
                                   const std::vector<species>& s2,
                                   float maximum_time);

std::string writeTREE_3(const std::vector<species> v,
                        float maximum_time);
std::string acquire_offspring_strings(const std::vector<species>& v,
                                      const species& focal,
                                      float maximum_time);

void remove_extinct_branches(std::vector<species>& all_species);
void merge_single_branches(std::vector<species>& all_species);


// new R functions update
std::string do_run_r(const Rcpp::NumericVector& parameters,
                     const Rcpp::NumericVector& waterlevel_changes,
                     float maximum_time,
                     int max_lin,
                     Rcpp::NumericMatrix& l_table,
                     rnd_t& rndgen);

int run(const Rcpp::NumericVector& parameters,
        const Rcpp::NumericVector& W,
        std::vector<species>& allSpecies,
        float maximum_time,
        int max_lin,
        rnd_t& rndgen);

void jiggle(std::vector< species > & s1,
            std::vector< species > & s2,
            float maximum_time,
            float jiggle_amount,
            rnd_t& rndgen);

Rcpp::NumericMatrix create_l_table( const std::vector< species > & s1,
                                    const std::vector< species > & s2);

int find_index_in_species_vector(const std::vector<species>& v,
                                 int id);

bool verify_consistency(const std::vector<species>& pop,
                        const std::vector<species>& extinct_species,
                        const std::string stage);

#endif

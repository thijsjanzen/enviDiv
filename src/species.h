//
//  species.h
//
//
//  Created by Thijs Janzen on 20/02/2019.
//

#ifndef species_h
#define species_h

#include <vector>
#include "Rcpp.h"

// [[Rcpp::plugins(cpp11)]]

struct species
{
    species();
    species(int& id);
    species(const species& parent_species, int& id_count, double b_time);
    species(const species& other);

    // move assignment

    int ID;
    double death_time; //time of death;

    int parent;
    bool has_extant_offspring;
    bool checked;

    bool alloSpeciated;

    void updateHistory(double t);
    bool check_has_viable_offspring(std::vector<species>& v);

    bool operator< (const species& other) const;
    bool operator> (const species& other) const;
    bool operator<=(const species& other) const;
    bool operator>=(const species& other) const;

    bool operator==(const species& other) const;
    bool operator!=(const species& other) const;

    species& operator=(const species& other);

    void set_birth_time(double t) {
      if(t < 0) {
        Rcpp::Rcout << "set_birth_time < 0!\n";
      }
      birth_time = t;
    }
    double get_birth_time() const {
      return birth_time;
    }

  private:
    double birth_time; //time of birth;

};

#endif /* species_h */

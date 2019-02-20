//
//  species.h
//
//
//  Created by Thijs Janzen on 20/02/2019.
//

#ifndef species_h
#define species_h

#include <vector>

// [[Rcpp::plugins(cpp11)]]

struct species
{
    species();
    species(int& id);
    species(const species& parent_species, int& id_count, double b_time);
    //species(species&& other);
    species(const species& other);

    // move assignment
   // species& operator=(species&& other);


    int ID;
    double birth_time; //time of birth;
    double death_time; //time of death;

    int parent;
    bool extant_offspring;
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

    void swap(species& other);
};

namespace std {

    // inject specialization for species into namespace std
    template<>
    inline void swap(species& a, species& b)
    {
        a.swap(b);
    }

}

#endif /* species_h */

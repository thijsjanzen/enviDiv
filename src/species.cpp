#include "species.h"

// [[Rcpp::plugins(cpp11)]]

species::species()  {
    ID = -1;
    death_time = -1;
}

species::species(const species& other): ID(other.ID),
                                        birth_time(other.birth_time),
                                        death_time(other.death_time),
                                        parent(other.parent),
                                        extant_offspring(other.extant_offspring),
                                        checked(other.checked),
                                        alloSpeciated(other.alloSpeciated)
{
}

species::species(int& id) { //constructor

    ID = id;  //every species has a unique ID to keep track of it, also over time
    id++;
    birth_time = 0;
    death_time = -1;
    parent = -1;
    alloSpeciated = false;
    checked = false;
}

species& species::operator=(const species& other) {
    if(this == &other) return *this;

    ID = other.ID;
    birth_time = other.birth_time;
    death_time = other.death_time;
    parent = other.parent;
    alloSpeciated = other.alloSpeciated;
    extant_offspring = other.extant_offspring;
    checked = other.checked;

    return *this;
}

bool species::operator< (const species& other) const {return       ID < other.ID; }
bool species::operator> (const species& other) const {return       ID > other.ID;}
bool species::operator<=(const species& other) const {return !operator> (other);}
bool species::operator>=(const species& other) const {return !operator< (other);}

bool species::operator==(const species& other) const { return ID == other.ID; }
bool species::operator!=(const species& other) const {return !((*this) == other);}

species::species(const species& parent_species, int& id_count, double b_time)
{
    ID = id_count;
    id_count++;
    death_time = -1;
    birth_time = b_time;
    parent = parent_species.ID;
    alloSpeciated = false;
    checked = false;
}

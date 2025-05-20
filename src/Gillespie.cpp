#include "Gillespie.h"
#include <cmath>
#include "random_thijs.h"
#include <string>

#include <Rcpp.h>
using namespace Rcpp;

//' simulate a tree using environmental diversification
//' @param parameters a vector of the four paramaters of the model
//' @param waterlevel_changes a vector that indicates the time points of water level changes
//' @param seed pseudo-random number generator seed
//' @param crown_age crown age of the tree to be simulated
//' @param max_lin maximum number of lineages
//' @return newick string
//' @export
// [[Rcpp::export]]
List create_tree_cpp(std::vector<double> parameters,
                     std::vector<double> waterlevel_changes,
                     double crown_age,
                     int max_lin) {
  // read parameter values

  rnd_t rndgen;

  std::vector< std::array<double, 4>> l_table;

  std::string code = do_run_r(parameters,
                              waterlevel_changes,
                              crown_age,
                              max_lin,
                              l_table,
                              rndgen);

  Rcpp::NumericMatrix l_table_out(l_table.size(), 4);
  for (size_t i = 0; i < l_table.size(); ++i) {
    for (size_t j = 0; j < 4; ++j) {
      l_table_out(i, j) = l_table[i][j];
    }
  }


  return List::create( Named("code") = code,
                       Named("Ltable") = l_table_out);
}


std::string do_run_r(const std::vector<double>& parameters,
                     const std::vector<double>& waterlevel_changes,
                     double maximum_time,
                     int max_lin,
                     std::vector< std::array<double, 4>>& l_table,
                     rnd_t& rndgen)
{
  std::vector< species > s1;

  double jiggle_amount = parameters[4];
  rndgen.set_normal_trunc(0.0f, jiggle_amount);

  //int idCount = 0;
  std::vector < std::vector < double > > l_table1;
  int error_code = run(parameters, waterlevel_changes,
                        s1, maximum_time, max_lin,
                        rndgen);

  if (error_code == 0) {
    return "extinction";
  }

  if (error_code == 1e6) {
    return "overflow";
  }

  std::vector<species> s2;

  int error_code2 = run(parameters,
      waterlevel_changes,
      s2, maximum_time, max_lin,
      rndgen);

  if (error_code2 == 0) {
    return "extinction";
  }
  if (error_code2 == 1e6) {
    return "overflow";
  }

  jiggle(s1, s2, maximum_time, jiggle_amount, rndgen);

  l_table = create_l_table(s1, s2);

  std::string output = "success";
  return output;
}



int drawEvent(double E, double S, double A, rnd_t& rndgen) {
  // this is a rather naive implementation
  // but for such a low number of events it suffices.
  double sum = E + S + A;
  double events[3] = {E/sum, S/sum, A/sum};
  double r = rndgen.uniform();
  for (int i = 0; i < 3; ++i) {
    r -= events[i];
    if(r <= 0.f) return i;
  }
  return 2;
}

bool onlyInstance(const std::vector<species>& v,
                  int i,
                  std::vector< int >& indices) {
  indices.clear();
  if (v.size() < 2) return true;

  int local_ID = v[i].get_ID();

  int count = 0;
  for (std::vector<species>::const_iterator s = v.begin(); s != v.end(); ++s) {
    if((*s).get_ID() == local_ID) {
      indices.push_back(count);
    }
    count++;
  }
  if(indices.size() == 1) return true; //only found self.

  return false;
}

void extinction(std::vector<species>& v,
                std::vector<species>& extinct_species,
                double time,
                int wLevel,
                rnd_t& rndgen) {
  int i = rndgen.random_number((int)v.size());

  v[i].death_time = time;

  if (wLevel == 0) { //low water level, there might be two instances of the same species
    std::vector< int > indices; // not used here
    bool no_other_instance_in_other_pocket = onlyInstance(v,i, indices);

    if(no_other_instance_in_other_pocket) extinct_species.push_back(v[i]);
  }
  if(wLevel == 1) {//high water level, there is only one instance of this species
    extinct_species.push_back(v[i]);
  }

  v[i] = v.back();
  v.pop_back();
  return;
}


void Symp_speciation(std::vector<species>& v,
                     int& id_count,
                     std::vector<species>& extinct_species,
                     double time,
                     double waterTime,
                     int wLevel,
                     rnd_t& rndgen)  {
  if (wLevel == 1) {
    // high water level, this should be easy, just a split
    // pick random parent
    verify_consistency(v, extinct_species, "sym_spec_high_water_before");
    int i = rndgen.random_number((int)v.size());
    species offspring1 = species(v[i], id_count, time);
    species offspring2 = species(v[i], id_count, time);

    // kill parent:
    v[i].death_time = time;
    extinct_species.push_back(v[i]);
    v[i] = v.back();
    v.pop_back();

    v.push_back(offspring1);
    v.push_back(offspring2);

    verify_consistency(v, extinct_species, "sym_spec_high_water_after");

    return;
  }

  if (wLevel == 0) {
    int i = rndgen.random_number((int)v.size());
    std::vector< int > pair_indices;

    bool only_instance = onlyInstance(v, i, pair_indices);

    if (only_instance) {
      // no paired species in other pocket, "simple" diversification
      species offspring1 = species(v[i], id_count, time);
      species offspring2 = species(v[i], id_count, time);

      // kill parent:
      v[i].death_time = time;
      extinct_species.push_back(v[i]);
      v[i] = v.back();
      v.pop_back();

      v.push_back(offspring1);
      v.push_back(offspring2);


      verify_consistency(v, extinct_species, "sym_spec_low_water_only_instance");

      return;
    } else {
      // low water level, with paired species, we get something complicated:
      //       W_t          T
      //       --p3--------- p3
      //  --p1-|
      //       |      ------ c1
      //       --p2--|
      //              ------ c2


      // first, we "redefine" species1 and add it to the extinct species list:
      species parent1 = v[pair_indices[0] ];
      parent1.death_time = waterTime;
      extinct_species.push_back(parent1);

      // we then generate p2 and p3
      species parent2 = species(parent1, id_count, waterTime);
      species parent3 = species(parent1, id_count, waterTime);


      // then we generate children 1 and 2:
      species child1 = species(parent2, id_count, time);
      species child2 = species(parent2, id_count, time);
      // we kill parent 2
      parent2.death_time = time;
      extinct_species.push_back(parent2);

      //add the children and parent3 to the vector:
      v[pair_indices[0] ] = child1;
      v[pair_indices[1] ] = child2;
      v.push_back(parent3);

      verify_consistency(v, extinct_species, "sym_spec_low_water_NOTonly_instance");
      return;
    }
  }
  return;
}

void waterLevelChange(std::vector<species>& v, int& wLevel) {
  if (wLevel == 0) { //waterlevel is low, and rises
    removeDuplicates(v);
  }

  if (wLevel == 1) { //waterlevel is high and drops
    std::vector<species> copy = v;
    v.insert(v.end(),copy.begin(),copy.end());  //all species are distributed over the two pockets, e.g. doubled
    //std::move(copy.begin(),copy.end(),std::back_inserter(v));
  }

  wLevel = 1 - wLevel;
  return;
}

void Allo_speciation(std::vector<species>& v,
                     int& id_count,
                     double time,
                     double water_time,
                     const std::vector<allo_pair>& p,
                     std::vector<species>& extinct_species,
                     rnd_t& rndgen) {
  int i = rndgen.random_number( (int)p.size());
  allo_pair focal_pair = p[i];

  //allright, we have a pair and they are branching:
  //          wT     T
  //          ------ c1
  // -- p1 -- |
  //          ------ c2

  species parent1 = v[focal_pair.index_a]; //these two are identical, shouldn't matter
  parent1.death_time = water_time;
  extinct_species.push_back(parent1);

  species child1 = species(parent1, id_count, water_time);
  species child2 = species(parent1, id_count, water_time);

  v[focal_pair.index_a] = child1;
  v[focal_pair.index_b] = child2;

  return;
}

void updatePairs2(std::vector<species>& v, std::vector<allo_pair>& p) {
  p.clear();
  if (v.size() == 2) {
    if (v[0].get_ID() == v[1].get_ID()) {
      p.push_back(allo_pair(v[0].get_ID(), 0));
      p[0].index_b = 1;
    }
    return;
  }

  std::sort(v.begin(),v.end());

  for (std::size_t i = 1; i < v.size(); ++i) {
    if (v[i].get_ID() == v[i-1].get_ID()) {
      p.push_back( allo_pair(v[i].get_ID(),
                             (int)i-1,
                             (int)i));
      ++i;
    }
  }
  return;
}


int find_index_in_species_vector(const std::vector<species>& v,
                                 int id) {
  for (int i = 0; i < v.size(); ++i) {
    if (v[i].get_ID() == id) {
      return i;
    }
  }
  return -1;
}



bool verify_consistency(const std::vector<species>& pop,
                        const std::vector<species>& extinct_species,
                        const std::string stage)
{
  // pre-emptive return, disable if you want to debug!
   return true;

  std::vector< species > all_species;
  all_species.insert(all_species.end(), pop.begin(), pop.end());
  all_species.insert(all_species.end(), extinct_species.begin(), extinct_species.end());


  for (int i = 0; i < all_species.size(); ++i) {
    int parent_id = all_species[i].get_parent();

    if (parent_id != -1) { // -1 is reserved for the root.

      int index_parent = find_index_in_species_vector(all_species,
                                                      parent_id);
      if (index_parent == -1) {
        Rcout << stage.c_str() << "\n";
        Rcpp::stop(stage);
       // return false;
      }

    }
  }
  return true;
}


int run(const std::vector<double>& parameters,
        const std::vector<double>& W,
        std::vector<species>& allSpecies,
        double maximum_time,
        int max_lin,
        rnd_t& rndgen) {

  double extinction_rate      = parameters[0];
  double sympatric_high_water = parameters[1];
  double sympatric_low_water  = parameters[2];
  double allopatric_spec_rate = parameters[3];

  allSpecies.clear();
  int id_count = 1;

  std::vector<species> pop;
  pop.push_back(species(id_count));

  // int numberExtinctions = 0;

  double time = 0;
  int waterLevel = 1;

  int numberWlevelChanges = 0;
  if (W[0] == 0.f) {
    numberWlevelChanges++;
  }

  std::vector<species> extinct_species;
  std::vector<allo_pair> pairs;

  while (time < maximum_time)  {
    // verify_consistency(pop, extinct_species, "general_beforepairs");

    double Pe = extinction_rate * pop.size();

    double Ps = 0.0;

    if (waterLevel == 0) Ps = sympatric_low_water  * pop.size();
    if (waterLevel == 1) Ps = sympatric_high_water * pop.size();

    if (waterLevel == 0) updatePairs2(pop, pairs); // low water

    double Pa = (1 - waterLevel) * allopatric_spec_rate * (pairs.size() * 2);

    double rate = Pe + Ps + Pa;

    double timestep = rndgen.Expon(rate);

    double prev_t = time;
    time += timestep;

    // verify_consistency(pop, extinct_species, "general_afterpairs");
    if (time >= W[numberWlevelChanges] && prev_t < W[numberWlevelChanges])
    //if (time > W[numberWlevelChanges] && prev_t <= W[numberWlevelChanges])
    {
      time = W[numberWlevelChanges];
      numberWlevelChanges++;
      waterLevelChange(pop, waterLevel);
    } else {

      if (time > maximum_time) break;

      int event_chosen = drawEvent(Pe, Ps, Pa, rndgen);
      double time_of_previous_waterlevelchange = 0.0;
      if(numberWlevelChanges != 0)
        time_of_previous_waterlevelchange = W[numberWlevelChanges - 1];

      switch(event_chosen)
      {
        case 0:
          extinction(pop, extinct_species, time, waterLevel, rndgen);
       //   numberExtinctions++;
          // verify_consistency(pop, extinct_species, "extinction");
          break;
        case 1:
          Symp_speciation(pop, id_count, extinct_species, time, time_of_previous_waterlevelchange, waterLevel, rndgen);
          // verify_consistency(pop, extinct_species, "symp_spec");
          break;
        case 2:
          Allo_speciation(pop, id_count,time, time_of_previous_waterlevelchange,pairs, extinct_species, rndgen);
          // verify_consistency(pop, extinct_species, "allo_spec");
          break;
        }
    }

    if (pop.size() < 1) //everything is extinct
    {
      return 0;
    }

    if (pop.size() > max_lin * 2) { // too many species
      return 1e6;
    }
  }

  verify_consistency(pop, extinct_species, "time limit");
  removeDuplicates(pop); //we remove the duplicates because there might be duplicates due to a low water level stand, these are not interesting (so in effect, we force the simulation to end with high water)

  verify_consistency(pop, extinct_species, "remove_duplicates");

  allSpecies.swap(pop);
  allSpecies.insert(allSpecies.end(), extinct_species.begin(), extinct_species.end());

  return 1;
}


///////////////////////////////////////////////////////////////////////////////////
//////////////////////// MEMBER FUNCTIONS /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

species::species(int& id) { //constructor
  ID = id;  //every species has a unique ID to keep track of it, also over time
  id++;
  birth_time = 0;
  death_time = -1;
  parent = -1;
  checked = false;
}

species& species::operator=(const species& other) {
  if (this == &other) return *this;

  ID = other.ID;
  birth_time = other.birth_time;
  death_time = other.death_time;
  parent = other.parent;
  extant_offspring = other.extant_offspring;
  checked = other.checked;

  return *this;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void removeDuplicates(std::vector<T>& vec) {
  std::sort(vec.begin(), vec.end()); //first we have to sort, and we sort on the ID of species
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); //species with an identical ID are removed (as we assume that they have had the same history)
}

species::species(const species& parent_species, int& id_count, double b_time) {
  ID = id_count;
  id_count++;
  death_time = -1;
  birth_time = b_time;
  parent = parent_species.ID;
  checked = false;
}


void jiggle_species_vector( std::vector< species > & s,
                            double focal_time,
                            double maximum_time,
                            double jiggle_amount,
                            rnd_t& rndgen) {

  for (auto brother = s.begin(); brother != s.end(); ++brother) {
    if ((*brother).birth_time == focal_time ) {

      // find sister species
      auto sister = brother;
      for ( auto jt = brother; jt != s.end(); ++jt) {
        if ((*jt).birth_time == (*brother).birth_time) {
          if ( (*jt).get_parent() == (*brother).get_parent()) {
            if ( (*jt).get_ID() != (*brother).get_ID()) { //because we start at self, we will always have an immediate hit
              sister = jt;
              break;
            }
          }
        }
      }

      double dist_to_upper = 1e6;
      double dist_to_lower = maximum_time - (*brother).birth_time;

      //find parent for upper limit
      //they both have the same parent, so we only have to check against one child
      for (auto p = s.begin(); p != s.end(); ++p) {
        if( (*p).get_ID() == (*brother).get_parent()) {
          dist_to_upper = (*brother).birth_time - (*p).birth_time;
          break;
        }
      }

      //find youngest offspring for lower limit
      for (auto o = s.begin(); o != s.end(); ++o) {
        if ((*o).get_parent() == (*brother).get_ID()) {
          double diff = (*o).birth_time - (*brother).birth_time;
          if(diff < dist_to_lower) dist_to_lower = diff;
        }
        if ((*o).get_parent() == (*sister).get_ID()) {
          double diff = (*o).birth_time - (*sister).birth_time;
          if(diff < dist_to_lower) dist_to_lower = diff;
        }
      }

      //truncation should be shortest distance
      double trunc = dist_to_lower;
      if(dist_to_upper < dist_to_lower) trunc = dist_to_upper;

      double new_birth_time = focal_time + rndgen.trunc_norm_predif(trunc); //trunc_normal(0.0, jiggle_amount, trunc);
      (*brother).birth_time = new_birth_time;
      (*sister).birth_time = new_birth_time;

      // now, if the parent died giving birth to the brother and sister
      // we also have to update the parent
      for (auto p = s.begin(); p != s.end(); ++p) {
        if ( (*p).get_ID() == (*brother).get_parent()) {
          if ( (*p).death_time == focal_time) {
            (*p).death_time = new_birth_time;
          }
        }
      }
    }
  }
  return;
}

//new version
void jiggle(std::vector< species > & s1,
            std::vector< species > & s2,
            double maximum_time,
            double jiggle_amount,
            rnd_t& rndgen) {
  //we have to identify all multiples

  std::vector<double> b_times;
  for (auto it = s1.begin(); it != s1.end(); ++it) { //collect all branching times across both sides of the tree
    b_times.push_back((*it).birth_time);
  }
  for (auto it = s2.begin(); it != s2.end(); ++it) {
    b_times.push_back((*it).birth_time);
  }
  std::sort(b_times.begin(), b_times.end()); //sort them (needed for remove_unique)

  //now we need to find those branching times that occur > 2 times
  std::vector<double> focal_times;
  int counter = 0;
  for (int i = 1; i < (int)b_times.size(); ++i) {
    if (b_times[i] == b_times[i-1]) {
      counter++;
    } else {
      if (counter > 2) {
        focal_times.push_back(b_times[i-1]);
      }
      counter = 0;
    }
  }

  if (focal_times.size() > 0) {
    for (int i = 0; i < (int)focal_times.size(); ++i) {
      double focal_time = focal_times[i];

      if (focal_time > 0) {
        jiggle_species_vector(s1, focal_time, maximum_time, jiggle_amount, rndgen); //jiggle all species with that time
        jiggle_species_vector(s2, focal_time, maximum_time, jiggle_amount, rndgen);
      }
    }
  }
  return;
}



ltab create_l_table(
    const std::vector< species > & s1,
    const std::vector< species > & s2)
{
  ltab output(s1.size() + s2.size());

  int i = 0;
  for (auto it = s1.begin(); it != s1.end(); ++it, ++i) {
    output[i][0] = (*it).birth_time;
    output[i][1] = (*it).get_parent();
    output[i][2] = (*it).get_ID();
    output[i][3] = (*it).death_time;
  }

  // verify l_table
  for (auto it = s2.begin(); it != s2.end(); ++it, ++i) {
    output[i][0] = (*it).birth_time;
    output[i][1] = -1 * ((*it).get_parent());
    output[i][2] = -1 * ((*it).get_ID());
    output[i][3] = (*it).death_time;
  }

  return output;
}


/*  old debug code */

/*

 double max_l_table(const  NumericMatrix& v) {
 double max = -1.f;
 for(int i = 0; i < v.size(); ++i) {
 if(v(i, 2) > max) max = v(i, 2);
 }
 return max;
 }

 double min_l_table(const  NumericMatrix& v) {
 double min = 1e6f;
 for(int i = 0; i < v.size(); ++i) {
 if(v(i, 2) < min) min = v(i, 2);
 }
 return min;
 }

 int which_min_l_table(const  NumericMatrix& v) {
 double min = min_l_table(v);
 for(int i = 0; i < v.size(); ++i) {
 if(v(i, 2) == min) {
 return i;
 }
 }
 return -1;
 }

 int which_max_l_table(const  NumericMatrix& v) {
 double max = max_l_table(v);
 for(int i = 0; i < v.size(); ++i) {
 if(v(i, 2) == max) {
 return i;
 }
 }
 return -1;
 }

 int which_l_table(const  NumericMatrix& v,
 double x) {
 for(int i = 0; i < v.size(); ++i) {
 if(v(i, 2) == x) {
 return i;
 }
 }
 return -1;
 }


 bool verify_l_table(const NumericMatrix& l_table) {
 // follow the trail
 int focal_index = which_max_l_table(l_table);
 int parent = l_table(focal_index, 1);
 while(parent > 1) {
 focal_index = which_l_table(l_table, parent);
 if(focal_index < 0) {
 Rcout << "l_table positive broken!\n";
 return false;
 }
 //   Rcout << l_table(focal_index, 0) << " " << l_table(focal_index, 1) << " " << l_table(focal_index, 2) << " " <<
 //      l_table(focal_index, 3) << "\n";
 parent = l_table(focal_index, 1);
 }
 return true;
 }

 bool verify_l_table2(const NumericMatrix& l_table) {
 // follow the trail
 int focal_index = which_max_l_table(l_table);
 int parent = l_table(focal_index, 1);
 while(parent < 0) {
 focal_index = which_l_table(l_table, parent);
 if(focal_index < -1) {
 Rcout << "l_table negative broken!\n";
 return false;
 }
 //   Rcout << l_table(focal_index, 0) << " " << l_table(focal_index, 1) << " " << l_table(focal_index, 2) << " " <<
 //      l_table(focal_index, 3) << "\n";
 parent = l_table(focal_index, 1);
 }
 return true;
 }
*/

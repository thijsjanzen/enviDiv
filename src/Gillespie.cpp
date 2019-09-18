#include "Gillespie.h"
#include <math.h>
#include "random_thijs.h"

//#include <chrono>
//#include <thread>
#include <string>

#include <Rcpp.h>
using namespace Rcpp;

//' simulate a tree using environmental diversification
//' @param parameters a vector of the four paramaters of the model
//' @param waterlevel_changes a vector that indicates the time points of water level changes
//' @param seed pseudo-random number generator seed
//' @param crown_age crown age of the tree to be simulated
//' @return newick string
//' @export
// [[Rcpp::export]]
List create_tree_cpp(std::vector<float> parameters,
                     std::vector<float> waterlevel_changes,
                     int seed,
                     float crown_age) {
  // read parameter values

  rnd_t rndgen;
  rndgen.set_seed(seed);


  NumericMatrix l_table;

  std::string code = do_run_r(parameters,
                              waterlevel_changes,
                              crown_age,
                              l_table,
                              rndgen);

  return List::create( Named("code") = code,
                       Named("Ltable") = l_table);
}

std::string do_run_r(const std::vector< float >& parameters,
                     const std::vector< float >& waterlevel_changes,
                     float maximum_time,
                     NumericMatrix& l_table,
                     rnd_t& rndgen)
{
  std::vector< species > s1;

  float jiggle_amount = parameters[4];
  rndgen.set_normal_trunc(0.0f, jiggle_amount);

  //int idCount = 0;
  std::vector < std::vector < float > > l_table1;
  int error_code = run(parameters, waterlevel_changes,
                        s1, maximum_time,
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
      s2, maximum_time,
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

int drawEvent(float E, float S, float A, rnd_t& rndgen) {
  // this is a rather naive implementation
  // but for such a low number of events it suffices.
  float sum = E + S + A;
  float events[3] = {E/sum, S/sum, A/sum};
  float r = rndgen.uniform();
  for (int i = 0; i < 3; ++i) {
    r -= events[i];
    if(r <= 0) return i;
  }
  return 0;
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
                float time,
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
                     float time,
                     float waterTime,
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
    v.insert(v.end(),copy.begin(),copy.end());  //all species are distributed over the two pockets, e.g. floatd
    //std::move(copy.begin(),copy.end(),std::back_inserter(v));
  }

  wLevel = 1 - wLevel;
  return;
}

void Allo_speciation(std::vector<species>& v,
                     int& id_count,
                     float time,
                     float water_time,
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


int run(const std::vector<float>& parameters,
        const std::vector<float>& W,
        std::vector<species>& allSpecies,
        float maximum_time,
        rnd_t& rndgen) {

  float extinction_rate      = parameters[0];
  float sympatric_high_water = parameters[1];
  float sympatric_low_water  = parameters[2];
  float allopatric_spec_rate = parameters[3];

  allSpecies.clear();
  int id_count = 1;
  std::vector<species> pop;
  pop.push_back(species(id_count));

  int numberExtinctions = 0;

  float time = 0;
  int waterLevel = 1;

  int iter = 0;

  int numberWlevelChanges = 0;

  std::vector<species> extinct_species;
  std::vector<allo_pair> pairs;

  while (time < maximum_time)  {
    verify_consistency(pop, extinct_species, "general_beforepairs");
    iter ++;

    float Pe = extinction_rate * pop.size();

    float Ps = 0.0;

    if (waterLevel == 0) Ps = sympatric_low_water * pop.size();
    if (waterLevel == 1) Ps = sympatric_high_water * pop.size();

    if (waterLevel == 0) updatePairs2(pop,pairs);

    float Pa = (1 - waterLevel) * allopatric_spec_rate * (pairs.size() * 2);

    float rate = Pe + Ps + Pa;

    float timestep = rndgen.Expon(rate);

    time += timestep;

    verify_consistency(pop, extinct_species, "general_afterpairs");

    if (time >= W[numberWlevelChanges] && (time-timestep) < W[numberWlevelChanges])
    {
      time = W[numberWlevelChanges];
      numberWlevelChanges++;
      waterLevelChange(pop, waterLevel);
    } else {

      if (time > maximum_time) break;

      int event_chosen = drawEvent(Pe, Ps, Pa, rndgen);
      float time_of_previous_waterlevelchange = 0.0;
      if(numberWlevelChanges != 0) time_of_previous_waterlevelchange = W[numberWlevelChanges-1];

      switch(event_chosen)
      {
        case 0:
          extinction(pop, extinct_species, time, waterLevel, rndgen);
          numberExtinctions++;

          verify_consistency(pop, extinct_species, "extinction");
          break;
        case 1:
          Symp_speciation(pop, id_count, extinct_species, time, time_of_previous_waterlevelchange, waterLevel, rndgen);
          verify_consistency(pop, extinct_species, "symp_spec");
          break;
        case 2:
          Allo_speciation(pop, id_count,time, time_of_previous_waterlevelchange,pairs, extinct_species, rndgen);
          verify_consistency(pop, extinct_species, "allo_spec");
          break;
        }
    }

    if (pop.size() < 1) //everything is extinct
    {
      return 0;
    }

    if (pop.size() > 300) //more than 300 species, unlikely to provide a good fit, but slows down the program considerably
    {
      return 1e6;
    }

  }

  verify_consistency(pop, extinct_species, "time limit");
  removeDuplicates(pop); //we remove the duplicates because there might be duplicates due to a low water level stand, these are not interesting (so in effect, we force the simulation to end with high water)

  verify_consistency(pop, extinct_species, "remove_duplicates");

  allSpecies.swap(pop);
  allSpecies.insert(allSpecies.end(), extinct_species.begin(), extinct_species.end());


  // remove_extinct_branches(allSpecies);
  // merge_single_branches(allSpecies);

  return 1;
}

void remove_extinct_branches(std::vector<species>& all_species) {

  std::vector< species > selected_species;
  for (auto it = all_species.begin(); it != all_species.end(); ++it) {

    (*it).check_has_viable_offspring(all_species);
    if ((*it).extant_offspring == true) {
      selected_species.push_back((*it));
    }
  }
  all_species = selected_species;
  return;
}

std::vector<int> find_indices(const std::vector<species>& v, int ID) {
  std::vector<int> indices;
  int count = 0;

  for (auto i = v.begin(); i != v.end(); ++i) {
    if((*i).get_parent() == ID) indices.push_back(count);
    count++;
  }

  return indices;
}

void merge_single_branches(std::vector<species>& all_species) {
  std::vector<species> species_vector = all_species;

  for (int i = 0; i < species_vector.size(); ++i) {
    std::vector<int> offspring = find_indices(species_vector,
                                              species_vector[i].get_ID());

    if (offspring.size() == 1) {
      // we have to do a merge
      species parent = species_vector[i];
      species daughter = species_vector[  offspring[0] ];

      species new_daughter = daughter;
      new_daughter.set_parent( parent.get_parent());
      new_daughter.birth_time = parent.birth_time;

      species_vector[i] = new_daughter;
      species_vector[ offspring[0]] = species_vector.back();
      species_vector.pop_back();

      i = 0;
      // and restart
    }
  }
  all_species = species_vector;
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

species::species(const species& parent_species, int& id_count, float b_time) {
  ID = id_count;
  id_count++;
  death_time = -1;
  birth_time = b_time;
  parent = parent_species.ID;
  checked = false;
}


bool species::check_has_viable_offspring(std::vector<species>& v) {
  if (checked == true) return extant_offspring;
  if (death_time == -1) {
    extant_offspring = true;
    checked = true;
    return extant_offspring;
  }

  std::vector<int> offspring = find_indices(v, ID); //find the positions of the offspring;
  extant_offspring = false;


  for (std::size_t i = 0; i < offspring.size(); ++i) {//ofspring is of size 2 (or 1), so using iterators is useless here
    int index = offspring[i];

    if (v[index].death_time == -1) extant_offspring  = true;  //the offspring is extant
    else   //the offspring has died, but might have given birth to other species that have survived
    {
      if (v[index].checked == true) extant_offspring  = v[index].extant_offspring;
      else extant_offspring = v[index].check_has_viable_offspring(v);
    }

    if (extant_offspring == true) break;
  }

  checked = true;

  return extant_offspring;
}

void jiggle_species_vector( std::vector< species > & s,
                            float focal_time,
                            float maximum_time,
                            float jiggle_amount,
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

      float dist_to_upper = 1e6;
      float dist_to_lower = maximum_time - (*brother).birth_time;

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
          float diff = (*o).birth_time - (*brother).birth_time;
          if(diff < dist_to_lower) dist_to_lower = diff;
        }
        if ((*o).get_parent() == (*sister).get_ID()) {
          float diff = (*o).birth_time - (*sister).birth_time;
          if(diff < dist_to_lower) dist_to_lower = diff;
        }
      }

      //truncation should be shortest distance
      float trunc = dist_to_lower;
      if(dist_to_upper < dist_to_lower) trunc = dist_to_upper;

      float new_birth_time = focal_time + rndgen.trunc_norm_predif(trunc); //trunc_normal(0.0, jiggle_amount, trunc);
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
            float maximum_time,
            float jiggle_amount,
            rnd_t& rndgen) {
  //we have to identify all multiples

  std::vector<float> b_times;
  for (auto it = s1.begin(); it != s1.end(); ++it) { //collect all branching times across both sides of the tree
    b_times.push_back((*it).birth_time);
  }
  for (auto it = s2.begin(); it != s2.end(); ++it) {
    b_times.push_back((*it).birth_time);
  }
  std::sort(b_times.begin(), b_times.end()); //sort them (needed for remove_unique)

  //now we need to find those branching times that occur > 2 times
  std::vector<float> focal_times;
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
      float focal_time = focal_times[i];

      if (focal_time > 0) {
        jiggle_species_vector(s1, focal_time, maximum_time, jiggle_amount, rndgen); //jiggle all species with that time
        jiggle_species_vector(s2, focal_time, maximum_time, jiggle_amount, rndgen);
      }
    }
  }
  return;
}



float max_l_table(const  NumericMatrix& v) {
  float max = -1.f;
  for(int i = 0; i < v.size(); ++i) {
    if(v(i, 2) > max) max = v(i, 2);
  }
  return max;
}

float min_l_table(const  NumericMatrix& v) {
  float min = 1e6f;
  for(int i = 0; i < v.size(); ++i) {
    if(v(i, 2) < min) min = v(i, 2);
  }
  return min;
}

int which_min_l_table(const  NumericMatrix& v) {
  float min = min_l_table(v);
  for(int i = 0; i < v.size(); ++i) {
    if(v(i, 2) == min) {
      return i;
    }
  }
  return -1;
}

int which_max_l_table(const  NumericMatrix& v) {
  float max = max_l_table(v);
  for(int i = 0; i < v.size(); ++i) {
    if(v(i, 2) == max) {
      return i;
    }
  }
  return -1;
}

int which_l_table(const  NumericMatrix& v,
                  float x) {
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

NumericMatrix create_l_table( const std::vector< species > & s1,
                              const std::vector< species > & s2)
{
  NumericMatrix output(s1.size() + s2.size(), 4);

  NumericMatrix temp(s1.size(), 4);

  int i = 0;
  for (auto it = s1.begin(); it != s1.end(); ++it, ++i) {
    output(i, 0) = (*it).birth_time;
    output(i, 1) = (*it).get_parent();
    output(i, 2) = (*it).get_ID();
    output(i, 3) = (*it).death_time;

    temp(i, _) = output(i, _);
  }

  // verify l_table
  verify_l_table(temp);


  NumericMatrix temp2(s2.size(), 4);
  for (auto it = s2.begin(); it != s2.end(); ++it, ++i) {
    output(i, 0) = (*it).birth_time;
    output(i, 1) = -1 * ((*it).get_parent());
    output(i, 2) = -1 * ((*it).get_ID());
    output(i, 3) = (*it).death_time;

    temp2(i - s1.size(), _) = output(i, _);
  }

  verify_l_table2(temp2);

  return output;
}




/*
std::vector< std::vector< float > > create_l_table(
    const std::vector< species > & s1,
    const std::vector< species > & s2)
{
  std::vector< std::array<float, 4 > > output(s1.size() + s2.size());
  std::array<float, 4 > to_add;
  int i = 0;
  for(auto it = s1.begin(); it != s1.end(); ++it, ++i) {
    to_add[0] = (*it).birth_time;
    to_add[1] = 2 + (*it).get_parent();
    to_add[2] = 2 + (*it).get_ID();
    to_add[3] = (*it).death_time;

    output[i] = to_add;
  }
  for(auto it = s2.begin(); it != s2.end(); ++it, ++i) {
    to_add[0] = (*it).birth_time;

    int parent = 2 + (*it).get_parent();
    if(parent > 0) parent *= -1.0;
    to_add[1] = parent;

    int id = 2 + (*it).get_ID();
    if(id > 0) id *= -1.0;
    to_add[2] = id;
    to_add[3] = (*it).death_time;

    output[i] = to_add;
  }

  return output;
}*/

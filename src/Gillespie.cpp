#include "Gillespie.h"
#include <math.h>
#include "random_thijs.h"

#include <Rcpp.h>
using namespace Rcpp;

//' test exponential generator
//' @param lambda input parameter
//' @return number drawn from exponential distribution
//' @export
// [[Rcpp::export]]
float generate_expon(float lambda)
{
  return Expon(lambda);
}

//' simulate a tree using environmental diversification
//' @param parameters a vector of the four paramaters of the model
//' @param waterlevel_changes a vector that indicates the time points of water level changes
//' @param seed pseudo-random number generator seed
//' @param crown_age crown age of the tree to be simulated
//' @return newick string
//' @export
// [[Rcpp::export]]
std::string create_tree_cpp(std::vector<float> parameters,
                            std::vector<float> waterlevel_changes,
                            int seed,
                            float crown_age) {
  // read parameter values
  set_seed(seed);

  std::string newick_tree = do_run_r(parameters,
                                     waterlevel_changes,
                                     crown_age);
  return newick_tree;
}

std::string do_run_r(const std::vector< float >& parameters,
                     const std::vector< float >& waterlevel_changes,
                     float maximum_time)
{
  std::vector<species> s1;

  float jiggle_amount = parameters[4];

  int idCount = 0;
  int error_code = run(parameters, waterlevel_changes,
                       idCount, s1,
                       maximum_time);

  if(error_code == 0) {
    return "extinction";
  }

  if(error_code == 1e6) {
    return "overflow";
  }

  std::vector<species> s2;

  int error_code2 = run(parameters,
      waterlevel_changes,
      idCount, s2,
      maximum_time);

  if(error_code2 == 0) {
    return "extinction";
  }
  if(error_code2 == 1e6) {
    return "overflow";
  }

  jiggle(s1, s2, maximum_time, jiggle_amount);


  std::string output = create_newick_string_r(s1, s2, maximum_time);

  return output;
}

int drawEvent(float E, float S, float A) {
  // this is a rather naive implementation
  // but for such a low number of events it suffices.
  float sum = E + S + A;
  float events[3] = {E/sum, S/sum, A/sum};
  float r = uniform();
  //Rcout << events[0] << " " << events[1] << " " << events[2] << "\n";
  for(int i = 0; i < 3; ++i) {
    r -= events[i];
    if(r <= 0) return i;
  }
  return 0;
}

bool onlyInstance(const std::vector<species>& v,
                  int i,
                  std::vector< int >& indices) {
  indices.clear();
  if(v.size() < 2) return true;

  int local_ID = v[i].get_ID();

  int count = 0;
  for(std::vector<species>::const_iterator s = v.begin(); s != v.end(); ++s) {
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
                int wLevel) {
  int i = random_number((int)v.size());

  v[i].death_time = time;
  std::vector< int > indices; // not used here
  if(wLevel == 0) { //low water level, there might be two instances of the same species

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
                     std::vector< float >& specTimes,
                     int wLevel)  {

  if(wLevel == 1) {
    // high water level, this should be easy, just a split
    // pick random parent
    int i = random_number((int)v.size());
    species offspring1 = species(v[i], id_count, time);
    species offspring2 = species(v[i], id_count, time);

    // kill parent:
    v[i].death_time = time;
    extinct_species.push_back(v[i]);
    v[i] = v.back();
    v.pop_back();

    v.push_back(offspring1);
    v.push_back(offspring2);
    return;
  }

  if(wLevel == 0) {
    int i = random_number((int)v.size());
    std::vector< int > pair_indices;

    bool only_instance = onlyInstance(v, i, pair_indices);

    if(only_instance) {
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
      return;
    } else {
      // low water level, with paired species, we get something complicated:
      //       W_t          T
      //       --p3--------- p3
      //  --p1-|
      //       |      ------ c1
      //       --p2--|
      //              ------ c2

      // we first generate p2 and p3
      species parent2 = species(v[pair_indices[0]], id_count, waterTime);
      species parent3 = species(v[pair_indices[1]], id_count, waterTime);

      // we generate children 1 and 2:
      species child1 = species(parent2, id_count, time);
      species child2 = species(parent2, id_count, time);
      // we kill parent 2
      parent2.death_time = time;
      extinct_species.push_back(parent2);

      //add the children and parent3 to the vector:
      v[pair_indices[0] ] = child1;
      v[pair_indices[1] ] = child2;
      v.push_back(parent3);
      return;
    }
  }
  return;
}

void waterLevelChange(std::vector<species>& v, int& wLevel) {
  if(wLevel == 0) { //waterlevel is low, and rises
    removeDuplicates(v);
  }

  if(wLevel == 1) { //waterlevel is high and drops
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
                     std::vector<float>& specTimes,
                     std::vector<species>& extinct_species)
{
  int i = random_number( (int)p.size());
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
  if(v.size() == 2) {
    if(v[0].get_ID() == v[1].get_ID()) {
      p.push_back(allo_pair(v[0].get_ID(), 0));
      p[0].index_b = 1;
    }
    return;
  }

  std::sort(v.begin(),v.end());

  for(std::size_t i = 1; i < v.size(); ++i) {
    if(v[i].get_ID() == v[i-1].get_ID()) {
      p.push_back( allo_pair(v[i].get_ID(),
                             (int)i-1,
                             (int)i));
      ++i;
    }
  }
  return;
}

int run(const std::vector<float>& parameters,
        const std::vector<float>& W,
        int& id_count,
        std::vector<species>& allSpecies,
        float maximum_time) {

  float extinction_rate      = parameters[0];
  float sympatric_high_water = parameters[1];
  float sympatric_low_water  = parameters[2];
  float allopatric_spec_rate = parameters[3];

  std::vector<float> speciationCompletionTimes;

  allSpecies.clear();

  std::vector<species> pop;
  pop.push_back(species(id_count));

  int numberExtinctions = 0;

  float time = 0;
  int waterLevel = 1;

  int iter = 0;

  int numberWlevelChanges = 0;

  std::vector<species> extinct_species;
  std::vector<allo_pair> pairs;

  while(time < maximum_time)  {

    iter ++;

    float Pe = extinction_rate * pop.size();

    float Ps = 0.0;

    if(waterLevel == 0) Ps = sympatric_low_water * pop.size();
    if(waterLevel == 1) Ps = sympatric_high_water * pop.size();

    if(waterLevel == 0) updatePairs2(pop,pairs);

    float Pa = (1 - waterLevel) * allopatric_spec_rate * (pairs.size() * 2);

    float rate = Pe + Ps + Pa;

    float timestep = Expon(rate);

    time += timestep;



    if(time >= W[numberWlevelChanges] && (time-timestep) < W[numberWlevelChanges])
    {
      time = W[numberWlevelChanges];
      numberWlevelChanges++;
      waterLevelChange(pop, waterLevel);
    } else {

      if(time > maximum_time) break;

      int event_chosen = drawEvent(Pe, Ps, Pa);
      float time_of_previous_waterlevelchange = 0.0;
      if(numberWlevelChanges != 0) time_of_previous_waterlevelchange = W[numberWlevelChanges-1];

      switch(event_chosen)
      {
      case 0:
        extinction(pop, extinct_species, time, waterLevel);
        numberExtinctions++;
        break;
      case 1:
        Symp_speciation(pop, id_count, extinct_species, time, time_of_previous_waterlevelchange,speciationCompletionTimes, waterLevel);
        break;
      case 2:
        Allo_speciation(pop, id_count,time, time_of_previous_waterlevelchange,pairs, speciationCompletionTimes, extinct_species);
        break;
      }
    }

    if(pop.size() < 1) //everything is extinct
    {
      return 0;
    }

    if(pop.size() > 300) //more than 300 species, unlikely to provide a good fit, but slows down the program considerably
    {
      return 1e6;
    }

  }

  removeDuplicates(pop); //we remove the duplicates because there might be duplicates due to a low water level stand, these are not interesting (so in effect, we force the simulation to end with high water)

  allSpecies.clear();
  allSpecies = pop;
  allSpecies.insert(allSpecies.end(), extinct_species.begin(), extinct_species.end());

  remove_extinct_branches(allSpecies);

  allSpecies =  merge_single_branches(allSpecies);

  return 1;
}


///////////////////////////////////////////////////////////////////////////////////
//////////////////////// MEMBER FUNCTIONS /////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

species::species(int& id) //constructor
{
  ID = id;  //every species has a unique ID to keep track of it, also over time
  id++;
  birth_time = 0;
  death_time = -1;
  parent = -1;
  checked = false;
}

species& species::operator=(const species& other)
{
  if(this == &other) return *this;

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
void removeDuplicates(std::vector<T>& vec)
{
  std::sort(vec.begin(), vec.end()); //first we have to sort, and we sort on the ID of species
  vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); //species with an identical ID are removed (as we assume that they have had the same history)
}

species::species(const species& parent_species, int& id_count, float b_time)
{
  ID = id_count;
  id_count++;
  death_time = -1;
  birth_time = b_time;
  parent = parent_species.ID;
  checked = false;
}


std::vector<int> find_indices(const std::vector<species>& v, int ID)
{
  std::vector<int> indices;
  int count = 0;

  for(auto i = v.begin(); i != v.end(); ++i) {
    if((*i).get_parent() == ID) indices.push_back(count);
    count++;
  }

  return indices;
}

bool species::check_has_viable_offspring(std::vector<species>& v)
{

  if(checked == true) return extant_offspring;
  if(death_time == -1) {
    extant_offspring = true;
    checked = true;
    return extant_offspring;
  }

  std::vector<int> offspring = find_indices(v, ID); //find the positions of the offspring;
  extant_offspring = false;


  for(std::size_t i = 0; i < offspring.size(); ++i) //ofspring is of size 2 (or 1), so using iterators is useless here
  {
    int index = offspring[i];

    if(v[index].death_time == -1) extant_offspring  = true;  //the offspring is extant
    else   //the offspring has died, but might have given birth to other species that have survived
    {
      if(v[index].checked == true) extant_offspring  = v[index].extant_offspring;
      else extant_offspring = v[index].check_has_viable_offspring(v);
    }

    if(extant_offspring == true) break;
  }

  checked = true;

  return extant_offspring;
}

std::vector<int> findOffspring(int ID, const std::vector<species>& v)
{
  std::vector<int> output;
  int counter = 0;

  for(auto i = v.begin(); i != v.end(); ++i) {
    if((*i).get_parent() == ID) output.push_back(counter);
    counter++;
  }

  return output;
}

std::string acquire_offspring_strings(const std::vector<species>& v,
                                      const species& focal,
                                      float maximum_time) {

  std::vector<int> offspring = findOffspring(focal.get_ID(), v);

  if(offspring.size() > 0) {
    std::string output = "(";
    for(int i = 0; i < offspring.size(); ++i) {
      output += acquire_offspring_strings(v,
                                          v[offspring[i]],
                                           maximum_time);
      if(i != offspring.size() -1) output += ",";
    }
    output += ")";

    float bl = focal.death_time - focal.birth_time;
    if(focal.death_time == -1) bl = maximum_time - focal.birth_time;
    output += ":" + std::to_string(bl);
    return output;
  } else {
    std::string output = std::to_string(focal.get_ID()) + ":";
    float bl = maximum_time - focal.birth_time;
    if(focal.death_time > 0)    bl = focal.death_time - focal.birth_time;
    output += std::to_string(bl);
    return output;
  }
}


std::string writeTREE_3(const std::vector<species> v,
                        float maximum_time) {
  //find the root
  // std::cout << "looking for root\n";
  species root;
  for(auto it = v.begin(); it != v.end(); ++it) {
    if((*it).get_parent() == -1 || (*it).birth_time == 0) {
      root = (*it);
      // std::cout << "found root, let's dive\n";
      break;
    }
  }

  // we start with the root
  std::string output = "(";
  std::string center = acquire_offspring_strings(v, root, maximum_time);
  output += center;
  output += ")";

  return output;
}

std::string create_newick_string_r(const std::vector<species>& s1,
                                   const std::vector<species>& s2,
                                   float maximum_time) {
  std::string left = writeTREE_3(s1, maximum_time);
  std::string right = writeTREE_3(s2, maximum_time);

  std::string output = "(";
  output += left;
  output += ",";
  output += right;
  output += ");";

  return output;
}

void remove_extinct_branches(std::vector<species>& all_species) {

  std::vector< species > selected_species;
  for(auto it = all_species.begin(); it != all_species.end(); ++it) {
    (*it).check_has_viable_offspring(all_species);
    if((*it).extant_offspring == true) {
      selected_species.push_back((*it));
    }
  }
  all_species = selected_species;
  return;
}

std::vector<species> merge_single_branches(const std::vector<species>& all_species) {
  std::vector<species> species_vector = all_species;

  for(int i = 0; i < species_vector.size(); ++i) {

    std::vector<int> offspring = find_indices(species_vector,
                                              species_vector[i].get_ID());

    if(offspring.size() == 1) {
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
  return species_vector;
}


void jiggle_species_vector( std::vector< species > & s,
                            float focal_time,
                            float maximum_time,
                            float jiggle_amount) {

  for(auto brother = s.begin(); brother != s.end(); ++brother) {
    if((*brother).birth_time == focal_time ) {

      // find sister species
      auto sister = brother;
      for( auto jt = brother; jt != s.end(); ++jt) {
        if((*jt).birth_time == (*brother).birth_time) {
          if( (*jt).get_parent() == (*brother).get_parent()) {
            if( (*jt).get_ID() != (*brother).get_ID()) { //because we start at self, we will always have an immediate hit
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
      for(auto p = s.begin(); p != s.end(); ++p) {
        if( (*p).get_ID() == (*brother).get_parent()) {
          dist_to_upper = (*brother).birth_time - (*p).birth_time;
          break;
        }
      }

      //find youngest offspring for lower limit
      for(auto o = s.begin(); o != s.end(); ++o) {
        if((*o).get_parent() == (*brother).get_ID()) {
          float diff = (*o).birth_time - (*brother).birth_time;
          if(diff < dist_to_lower) dist_to_lower = diff;
        }
        if((*o).get_parent() == (*sister).get_ID()) {
          float diff = (*o).birth_time - (*sister).birth_time;
          if(diff < dist_to_lower) dist_to_lower = diff;
        }
      }

      //truncation should be shortest distance
      float trunc = dist_to_lower;
      if(dist_to_upper < dist_to_lower) trunc = dist_to_upper;

      float new_birth_time = focal_time + trunc_normal(0.0, jiggle_amount, trunc);
      (*brother).birth_time = new_birth_time;
      (*sister).birth_time = new_birth_time;

      // now, if the parent died giving birth to the brother and sister
      // we also have to update the parent
      for(auto p = s.begin(); p != s.end(); ++p) {
        if( (*p).get_ID() == (*brother).get_parent()) {
          if( (*p).death_time == focal_time) {
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
            float jiggle_amount) {
  //we have to identify all multiples

  std::vector<float> b_times;
  for(auto it = s1.begin(); it != s1.end(); ++it) { //collect all branching times across both sides of the tree
    b_times.push_back((*it).birth_time);
  }
  for(auto it = s2.begin(); it != s2.end(); ++it) {
    b_times.push_back((*it).birth_time);
  }
  std::sort(b_times.begin(), b_times.end()); //sort them (needed for remove_unique)

  //now we need to find those branching times that occur > 2 times
  std::vector<float> focal_times;
  int counter = 0;
  for(int i = 1; i < (int)b_times.size(); ++i) {
    if(b_times[i] == b_times[i-1]) {
      counter++;
    } else {
      if(counter > 2) {
        focal_times.push_back(b_times[i-1]);
      }
      counter = 0;
    }
  }

  if(focal_times.size() > 0) {
    for(int i = 0; i < (int)focal_times.size(); ++i) {
      float focal_time = focal_times[i];

      if(focal_time > 0) {
        jiggle_species_vector(s1, focal_time, maximum_time, jiggle_amount); //jiggle all species with that time
        jiggle_species_vector(s2, focal_time, maximum_time, jiggle_amount);
      }
    }
  }
  return;
}

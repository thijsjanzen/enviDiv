#include "Gillespie.h"
#include "species.h"
#include "spec_point.h"
#include "allo_pair.h"
#include "newick_node.h"
#include "find.h"
#include <math.h>

#include <Rcpp.h>
using namespace Rcpp;

//' bogus testing function
//' @param input double number
//' @return square of input
//' @export
// [[Rcpp::export]]
double square_number(double input) {
  return input*input;
}

//' simulate a tree using environmental diversification
//' @param parameters a vector of the four paramaters of the model
//' @param waterlevel_changes a vector that indicates the time points of water level changes
//' @param seed pseudo-random number generator seed
//' @param crown_age crown age of the tree to be simulated
//' @return newick string
//' @export
// [[Rcpp::export]]
std::string create_tree_cpp(std::vector<double> parameters,
                            std::vector<double> waterlevel_changes,
                            int seed,
                            int crown_age) {
	// read parameter values
	set_seed(seed);
	std::string newick_tree = do_run(parameters,
                                   waterlevel_changes,
                                   crown_age);
  return newick_tree;
}

std::string do_run(std::vector<double> parameters,
                   std::vector<double> waterlevel_changes,
                   int maximum_time) {

  std::vector<species> s1;

 	int idCount = 1e5;
	int local_lins = run(parameters,
                       waterlevel_changes,
                       idCount, s1,
                       maximum_time);

	if(local_lins == 0) {
		return "extinction";
	}
	if(local_lins > 300) {
	  return "overflow";
	}

	std::vector<species> s2;

	idCount = 2e5;

	int local_lins2 = run(parameters,
                        waterlevel_changes,
                        idCount, s2,
                        maximum_time);

	if(local_lins == 0 || local_lins2 == 0) {
	  return "extinction";
	}
	if(local_lins > 300 || local_lins2 > 300) {
	  return "overflow";
	}

  int jiggle_amount = parameters[4];

	jiggle(s1, s2, maximum_time, jiggle_amount);

	std::string output =	create_newick_string_r(s1, s2, maximum_time);

	return output;
}

int run(const std::vector<double> parameters,
        const std::vector<double>& W,
        int& id_count,
        std::vector<species>& allSpecies,
        int maximum_time)  {

  double extinction_rate =      parameters[0];
  double sympatric_high_water = parameters[1];
  double sympatric_low_water  = parameters[2];
  double allopatric_spec_rate = parameters[3];

  std::vector<double> speciationCompletionTimes;
  allSpecies.clear();

  std::vector<spec_point> lineages;

  std::vector<species> pop;
  pop.push_back( species(id_count) );

  int numberExtinctions = 0;

  double time = 0;
  int water_level = 1;

  int iter = 0;

  int numberWlevelChanges = 0;
  std::vector<species> extinct_species;
  std::vector<species> allo_species;

  //int alloSpeciations = 0;
  std::vector<allo_pair> pairs;

  while(time < maximum_time)  {
    iter ++;
    double Pe = extinction_rate * pop.size();

    double Ps;
    if(water_level == 1) Ps = sympatric_high_water * pop.size();
    if(water_level == 0) Ps = sympatric_low_water * pop.size();

    if(water_level == 0) updatePairs(pop, pairs);

    double Pa = (1 - water_level) * allopatric_spec_rate * (pairs.size() * 2);

    double rate = Pe + Ps + Pa;

    double timestep = Expon(rate);

    time += timestep;
    if(time > maximum_time) break;

    if(time > W[numberWlevelChanges] && (time-timestep) < W[numberWlevelChanges]) {

      time = W[numberWlevelChanges];
      numberWlevelChanges++;
      waterLevelChange(pop, water_level);
    } else	{

      int event_chosen = drawEvent(Pe, Ps, Pa);

      double time_of_previous_waterlevelchange = 0.0;
      if(numberWlevelChanges != 0) time_of_previous_waterlevelchange = W[numberWlevelChanges-1];

      switch(event_chosen)
      {
        case 0:
          extinction(pop,
                     extinct_species,
                     time,
                     water_level);
          numberExtinctions++;
          break;
        case 1:
          Symp_speciation(pop,
                          id_count,
                          extinct_species,
                          time,
                          time_of_previous_waterlevelchange,
                          speciationCompletionTimes,
                          water_level);
          break;
        case 2:
          Allo_speciation(pop,
                          id_count,
                          time,
                          time_of_previous_waterlevelchange,
                          pairs,
                          speciationCompletionTimes,
                          extinct_species);
          break;
        }
    }

    if(pop.size() < 1) {//everything is extinct
      return pop.size();
    }

    if(pop.size() > 300) {//more than 300 species, unlikely to provide a good fit, but slows down the program considerably
      return pop.size();
    }
  }

  removeDuplicates(pop); //we remove the duplicates because there might be duplicates due to a low water level stand, these are not interesting (so in effect, we force the simulation to end with high water)

  allSpecies.clear();
  allSpecies = pop;
  //std::move(extinct_species.begin(),
  //            extinct_species.end(),
  //          std::back_inserter(allSpecies));
  for(auto it = extinct_species.begin(); it != extinct_species.end(); ++it) {
    allSpecies.push_back((*it));
  }

  if(numberExtinctions > 0) {
    remove_extinct_branches(allSpecies);
  }

  merge_single_branches(allSpecies);

  return pop.size();
}

bool onlyInstance(const std::vector<species>& v,
                  int i,
                  std::vector< int >& indices) {
  indices.clear();
  if(v.size() < 2) return true;

  int local_ID = v[i].ID;

  int count = 0;
  for(std::vector<species>::const_iterator s = v.begin(); s != v.end(); ++s) {
    if((*s).ID == local_ID) {
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
                     double time,
                     double waterTime,
                     std::vector<double>& specTimes,
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

void waterLevelChange(std::vector<species>& v,
                      int& wLevel) {
	if(wLevel == 0) { //waterlevel is low, and rises
		removeDuplicates(v);
	}

	if(wLevel == 1) { //waterlevel is high and drops
		std::vector<species> copy = v;
		v.insert(v.end(), copy.begin(), copy.end());  //all species are distributed over the two pockets, e.g. doubled
		//std::move(copy.begin(), copy.end(), std::back_inserter(v));
	}

	wLevel = 1 - wLevel;
	return;
}

void Allo_speciation(std::vector<species>& v,
                     int& id_count,
                     double time,
                     double water_time,
                     const std::vector<allo_pair>& p,
                     std::vector<double>& specTimes,
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


void updatePairs(std::vector<species>& v,
                  std::vector<allo_pair>& p)
{
	p.clear();
	if(v.size() == 2)
	{
		if(v[0].ID == v[1].ID)
		{
			p.push_back(allo_pair(v[0].ID, 0));
			p[0].index_b = 1;
		}
		return;
	}

	std::sort(v.begin(),v.end());

	for(std::size_t i = 1; i < v.size(); ++i)
	{
		if(v[i].ID == v[i-1].ID)
		{
			p.push_back(allo_pair(v[i].ID,i-1,i));
			++i;
		}
	}
	return;
}

template<typename T>
void removeDuplicates(std::vector<T>& vec)  {
	std::sort(vec.begin(), vec.end()); //first we have to sort, and we sort on the ID of species
	vec.erase(std::unique(vec.begin(), vec.end()), vec.end()); //species with an identical ID are removed (as we assume that they have had the same history)
}

int countLin(const std::vector<species>& v,
             double time) {
	int count = 0;
	for(auto it = v.begin(); it != v.end(); ++it) {
		if( (*it).get_birth_time() <= time) {
			if((*it).death_time == -1 || (*it).death_time > time) {
				count++;
			}
		}
	}

	return count;
}

bool sortOnTime(const species& left,
                const species& right) {
	return left.get_birth_time() < right.get_birth_time();
}

void purgeOutput(std::vector<spec_point>& v)  {
	std::vector<spec_point> new_vec;

	int previous_ID = -1;

	for(auto it = v.begin(); it != v.end(); ++it) {
		if((*it).ID != previous_ID) {
			new_vec.push_back((*it));
			previous_ID = (int)(*it).ID;
		}
	}

	v = new_vec;
	return;
}


int getIndex(int ID, const std::vector<species>& v) {
	for(std::size_t i = 0; i < v.size(); ++i) {
		if(v[i].ID == ID) return i;
	}
	return -1;
}


void updateReferences(int oldID,
                      int newID,
                      std::vector< newick_node >& v)  {
	for(std::size_t i = 0; i < v.size(); ++i)
	{
		if( v[i].parent == oldID) v[i].parent = newID;
	}
}

void updateReferences(int oldID,
                      int newID,
                      std::vector<species>& v,
                      double time)  {
	for(auto it = v.begin(); it != v.end(); ++it)  {
		if( (*it).parent == oldID)	{
			if((*it).get_birth_time() >= time) (*it).parent = newID;
		}
	}
}

std::string to_string_local( double d ) {

  std::ostringstream stm ;
  stm << std::setprecision(std::numeric_limits<double>::digits10) << d ;
  return stm.str() ;
}

std::string construct_string(const newick_node focal_node,
                             const std::vector<newick_node>& others) {

  if(focal_node.extant == true) {
    Rcout << "extant node, return to base\n";
    std::string output = to_string_local(focal_node.ID) + ":" + to_string_local(focal_node.get_branch_length());
    return output;
  }

  // not extant, so we have to find the children.
  Rcout << "not extant, find children\n";
  std::vector< int > children = findOffspring(focal_node.ID, others);
  if(children.size() != 2) {
    Rcout << "not 2 children, wtf: " << children.size() << "\n";
  }
  if(children[0] > others.size()) {
    Rcout << "child[0] > others.size()\n";
  }
  if(children[1] > others.size()) {
    Rcout << "child[1] > others.size()\n";
  }

  assert(children.size() == 2);
  assert(children[0] < others.size());
  assert(children[1] < others.size());

  std::string output =  "(";
  Rcout << "construct child string 1\n";
  output += construct_string(others[ children[0] ], others);
  output += ",";
  Rcout << "construct child string 2\n";
  output += construct_string(others[ children[1] ], others);
  output += ")";
  output += ":";
  output += std::to_string(focal_node.get_branch_length());
  Rcout << "return to base\n";
  return output;
}

std::string create_newick_string_local(const std::vector<newick_node>& v) {

  // find the root:
  newick_node root;
  Rcout << "looking for root\n";
  for(auto it = v.begin(); it != v.end(); ++it) {
    if((*it).parent == -1) {
      root = (*it);
      Rcout << "found root, let's dive\n";
      break;
    }
  }
  if(root.ID == -1) {
    Rcout << "FAILED to find root!\n";
  }

  std::string output = construct_string(root, v);

  return output;
}


bool species::check_has_viable_offspring(std::vector<species>& v) {
	if(checked == true) return has_extant_offspring;
	if(death_time == -1) {
	  has_extant_offspring = true;
	  checked = true;
	  return has_extant_offspring;
	}

	has_extant_offspring = false;
	std::vector<int> offspring = find_indices(v, ID); //find the positions of the offspring;

	for(std::size_t i = 0; i < offspring.size(); ++i) { //ofspring is of size 2 (or 1), so using iterators is useless here
		int index = offspring[i];

		if(v[index].death_time == -1) {
		  has_extant_offspring  = true;  //the offspring is extant
		} else  { //the offspring has died, but might have given birth to other species that have survived
			if(v[index].checked == true) has_extant_offspring  = v[index].has_extant_offspring;
			else has_extant_offspring = v[index].check_has_viable_offspring(v);
		}

		if(has_extant_offspring == true) break;
	}

	checked = true;

	return has_extant_offspring;
}

int countLineages(const std::vector<spec_point>& v, double time)  {
	for(std::size_t i = 0; i < v.size(); ++i) {
		if(v[i].get_time() == time) return (int)v[i].ID;

		if(i > 0) {
			if(v[i].get_time() > time && v[i-1].get_time() < time) return (int)v[i-1].ID;
		}
	}
	return (int)v.back().ID;
}

int drawEvent(double E, double S, double A) {
  double sum = E + S + A;
  double events[3] = {E/sum, S/sum, A/sum};
  double r = uniform();
  for(int i = 0; i < 3; ++i)
  {
    r -= events[i];
    if(r <= 0) return i;
  }
  return 0;
}

void jiggle_species_vector( std::vector< species > & s,
                            double focal_time,
                            double maximum_time,
                            double jiggle_amount) {

  for(auto brother = s.begin(); brother != s.end(); ++brother) {
    if((*brother).get_birth_time() == focal_time ) {

      // find sister species
      auto sister = brother;
      for( auto jt = brother; jt != s.end(); ++jt) {
        if((*jt).get_birth_time() == (*brother).get_birth_time()) {
          if( (*jt).parent == (*brother).parent) {
            if( (*jt).ID != (*brother).ID) { //because we start at self, we will always have an immediate hit
              sister = jt;
              break;
            }
          }
        }
      }

      double dist_to_upper = 1e6;
      double dist_to_lower = maximum_time - (*brother).get_birth_time();

      //find parent for upper limit
      //they both have the same parent, so we only have to check against one child
      for(auto p = s.begin(); p != s.end(); ++p) {
        if( (*p).ID == (*brother).parent) {
          dist_to_upper = (*brother).get_birth_time() - (*p).get_birth_time();
          break;
        }
      }

      //find youngest offspring for lower limit
      for(auto o = s.begin(); o != s.end(); ++o) {
        if((*o).parent == (*brother).ID) {
          double diff = (*o).get_birth_time() - (*brother).get_birth_time();
          if(diff < dist_to_lower) dist_to_lower = diff;
        }
        if((*o).parent == (*sister).ID) {
          double diff = (*o).get_birth_time() - (*sister).get_birth_time();
          if(diff < dist_to_lower) dist_to_lower = diff;
        }
      }

      //truncation should be shortest distance
      double trunc = dist_to_lower;
      if(dist_to_upper < dist_to_lower) trunc = dist_to_upper;

      double new_birth_time = focal_time + trunc_normal(0.0, jiggle_amount, trunc);
      (*brother).set_birth_time(new_birth_time);
      (*sister).set_birth_time(new_birth_time);
    }
  }
  return;
}



//new version
void jiggle(std::vector< species > & s1,
            std::vector< species > & s2,
            double maximum_time,
            double jiggle_amount) {
  //we have to identify all multiples

  std::vector<double> b_times;
  for(auto it = s1.begin(); it != s1.end(); ++it) { //collect all branching times across both sides of the tree
    b_times.push_back((*it).get_birth_time());
  }
  for(auto it = s2.begin(); it != s2.end(); ++it) {
    b_times.push_back((*it).get_birth_time());
  }
  std::sort(b_times.begin(), b_times.end()); //sort them (needed for remove_unique)

  //now we need to find those branching times that occur > 2 times
  std::vector<double> focal_times;
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
      double focal_time = focal_times[i];

      if(focal_time > 0) {
        jiggle_species_vector(s1, focal_time, maximum_time, jiggle_amount); //jiggle all species with that time
        jiggle_species_vector(s2, focal_time, maximum_time, jiggle_amount);
      }
    }
  }
  return;
}

std::string acquire_offspring_strings(const std::vector<species>& v,
                                      const species& focal,
                                      double maximum_time) {

    std::vector<int> offspring = findOffspring(focal.ID, v);

    if(offspring.size() > 0) {
        std::string output = "(";
        for(int i = 0; i < offspring.size(); ++i) {
            output += acquire_offspring_strings(v,
                                                v[offspring[i]],
                                                maximum_time);
            if(i != offspring.size() -1) output += ",";
        }
        output += ")";

        double bl = focal.death_time - focal.get_birth_time();
        if(focal.death_time == -1) bl = maximum_time - focal.get_birth_time();
        output += ":" + std::to_string(bl);
        return output;
    } else {
        std::string output = std::to_string(focal.ID) + ":";
        double bl = maximum_time - focal.get_birth_time();
        if(focal.death_time > 0)    bl = focal.death_time - focal.get_birth_time();
        output += std::to_string(bl);
        return output;
    }
}


std::string writeTREE_3(const std::vector<species> v,
                        double maximum_time) {
    //find the root
    // std::cout << "looking for root\n";
    species root;
    for(auto it = v.begin(); it != v.end(); ++it) {
        if((*it).parent == -1) {
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
                                   double maximum_time) {
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
        if((*it).has_extant_offspring == true) {
            selected_species.push_back((*it));
        }
    }
    all_species = selected_species;
    return;
}

void merge_single_branches(std::vector<species>& all_species) {
    std::vector<int> indices_to_copy;
    for(int i = 0; i < all_species.size(); ++i) {
        std::vector<int> offspring = find_indices(all_species, all_species[i].ID); //find the positions of the offspring;

        if(offspring.size() > 1) {
            indices_to_copy.push_back(i);
        }
        if(offspring.size() == 1) {
            // this species only gave rise to one other species
            // let's see if it is extant
            if(all_species[i].death_time == -1) {
                indices_to_copy.push_back(i);
            } else {
                // all right, so we have to merge the offspring with the parent.
                species parent_branch = all_species[i];
                species offspring_branch = all_species[ offspring[0] ];
                // now the birth time of the offspring has to be set to be the parents
                double parent_birth_time = all_species[i].get_birth_time();
                all_species[ offspring[0] ].set_birth_time(parent_birth_time);
            }
        }
        if(offspring.size() == 0) {
            if(all_species[i].death_time == -1) {
                indices_to_copy.push_back(i);
            } else {
                std::cout << "ERROR in merge branches\n";
            }
        }
    }

    std::vector< species > selected_species;


    for(int i = 0; i < indices_to_copy.size(); ++i) {
        selected_species.push_back(all_species [ indices_to_copy[i] ]);
    }
    all_species = selected_species;
}

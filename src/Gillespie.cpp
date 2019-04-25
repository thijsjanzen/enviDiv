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

 	int idCount = 0;
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


 // if(numberExtinctions == 0) lineages = calculateLineages_noextinct(allSpecies);
//  else
//  {
//   lineages = calculateLineages_withextinct(allSpecies, maximum_time);
//  }
  return pop.size();
}

bool onlyInstance(const std::vector< species >& v, int i) {
	if(v.size() < 2) return true;

	int local_ID = v[i].ID;
	int count = 0;
	for(auto it = v.begin(); it != v.end(); ++it) {
		if((*it).ID == local_ID) count++;

		if(count > 1) return false;
	}
	return true;
}

void extinction(std::vector<species>& v,
                std::vector<species>& extinct_species,
                double time,
                int wLevel) {

	int i = random_number(v.size());

	v[i].death_time = time;

	if(wLevel == 0) {//low water level, there might be two instances of the same species
		if(onlyInstance(v,i)) extinct_species.push_back(v[i]);
	}
	if(wLevel == 1) { //high water level, there is only one instance of this species
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
    int i = random_number(v.size());
    species offspring1 = species(v[i], id_count, time);
    species offspring2 = species(v[i], id_count, time);

    // kill parent:
    v[i].death_time = time;
    extinct_species.push_back(v[i]);
    v[i] = offspring1;
    v.push_back(offspring2);
    return;
  }

  if(wLevel == 0) {
    int i = random_number(v.size());
    bool only_instance = onlyInstance(v, i);

    if(only_instance) {
      // no paired species in other pocket, "simple" diversification
      species offspring1 = species(v[i], id_count, time);
      species offspring2 = species(v[i], id_count, time);

      // kill parent:
      v[i].death_time = time;
      extinct_species.push_back(v[i]);
      v[i] = offspring1;
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
      species parent2 = species(v[i], id_count, waterTime);
      species parent3 = species(v[i], id_count, waterTime);
      // we kill parent 1:
      v[i].death_time = waterTime;
      extinct_species.push_back(v[i]);
      // we replace parent1 with parent 3
      v[i] = parent3;

      // we generate children 1 and 2:
      species child1 = species(parent2, id_count, time);
      species child2 = species(parent2, id_count, time);
      // we kill parent 2
      parent2.death_time = time;
      extinct_species.push_back(parent2);

      //add the children to the vector:
      v.push_back(child1);
      v.push_back(child2);
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
	int i = random_number(p.size());
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

std::vector<newick_node> generateNodeList(const std::vector<species>& v,
                                          double maximum_time)  {
	std::vector<newick_node> node_list;

	std::vector<species> extant;
	std::vector<species> extinct;

	int maxID = 0;

	for(auto it = v.begin(); it != v.end(); ++it) {
		if((*it).ID > maxID) maxID = (*it).ID;

		if((*it).death_time == -1) {
		  extant.push_back((*it));
		} else {
		  extinct.push_back((*it));
		}
	}
	maxID++;

	//we have all extant species, that are connected through the extinct species to the ancestor (parent ID = -1);
	while(!extant.empty())  {
		//we find the instance with the most recent branching time, e.g. the most recent branching moment
		int youngest = findYoungest(extant);
		int other = findOther(youngest,extant);

		if(extant[youngest].parent == -1) //we have reached the root!
		{
			newick_node root(extant[0], maximum_time);
			node_list.push_back(root);
			break;
		}

		if(other != -1) //two species branched off from one other species
		{
			newick_node temp1(extant[youngest], maximum_time);
			newick_node temp2(extant[other], maximum_time);

			node_list.push_back(temp1);
			node_list.push_back(temp2);

			//remove both offspring, and add the parent
			int parent = find_parent(extinct,
                            extant[youngest].parent);

			if(parent != -1)
			{
				extant[youngest] = extinct[parent];
				extant[other] = extant.back();
				extant.pop_back();

				extinct[parent] = extinct.back();
				extinct.pop_back();
			}
			else
			{
				//the parent is already in the extant community
				std::vector<species> temp;
				for(std::size_t i = 0; i < extant.size(); ++i)
				{
					if(i != youngest && i != other) temp.push_back(extant[i]);
				}
				extant = temp;
			}
		}
		else
		{
			//either the species branched off from an extant species, or not:
			int parent = find_parent(extant,
                               extant[youngest].parent);

			if(parent == -1)
			{
				///////////////////////////////////
				/////////   1   /          1 /
				////////       /    -->     /
				////////   2  /          1 /
				///////////////////////////////////////
				parent = find_parent(extinct,
                             extant[youngest].parent);
				int oldID = extant[youngest].ID;
				int newID = extinct[parent].ID;

				double time = extant[youngest].get_birth_time();

				extinct[parent].death_time = extant[youngest].death_time;
				extant[youngest] = extinct[parent];
				extinct[parent] = extinct.back();
				extinct.pop_back();

				updateReferences(oldID,newID,extinct,time);
				updateReferences(oldID,newID,extant,time);
				updateReferences(oldID,newID,node_list);
			}
			else //I don't think this can happen, but just in case
			{
				if(extant[parent].death_time == extant[youngest].get_birth_time())
				{
					///////////////////////////////////
					/////////   1   /
					////////       /
					////////   2  /
					///////////////////////////////////////
					newick_node temp(extant[youngest], maximum_time);
					node_list.push_back(temp);

					extant[youngest] = extant.back();
					extant.pop_back();
				}
				else
				{
					///////////////////////////////////
					/////////   1   /           1  /
					////////       /\    -->      /\
					////////   1  /  \ 2       3 /  \ 2
					///////////////////////////////////////
					newick_node temp1(extant[youngest], maximum_time); //number 2
					newick_node temp2(extant[parent], maximum_time); //future number 3

					temp1.parent = extant[parent].ID;
					temp2.parent = extant[parent].ID; //they both have the same parents

					maxID++;
					temp2.ID = maxID;  //make it number 3
					temp2.set_branch_length(maximum_time - extant[youngest].get_birth_time()); //adjust branch length (shorten it)
					if(extant[parent].death_time != -1)
					  temp2.set_branch_length(extant[parent].death_time -
					                          extant[youngest].get_birth_time());


					extant[parent].death_time = extant[youngest].get_birth_time();

					updateReferences(extant[parent].ID, maxID, extant, extant[youngest].get_birth_time()); //adjust all downstream references
					updateReferences(extant[parent].ID, maxID, extinct, extant[youngest].get_birth_time());
					updateReferences(extant[parent].ID, maxID, node_list);

					node_list.push_back(temp1); //add
					node_list.push_back(temp2);

					extant[youngest] = extant.back();
					extant.pop_back();
				}
			}
		}
	}


	return node_list;
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



std::string writeTREE2(const std::vector<species>& v,
                       double maximum_time)
{
	std::vector<species> extant;
	std::vector<species> extinct;

	int maxID = 0;

	for(auto it = v.begin(); it != v.end(); ++it) {
		if((*it).ID > maxID) maxID = (*it).ID;

		if((*it).death_time == -1) {
		  extant.push_back((*it));
		} else {
		  extinct.push_back((*it));
		}
	}
	maxID++;

	if(extant.size() ==1) { //there is only one, or two species
		std::string s_ID = std::to_string(extant[0].ID);
		std::string s_BL = std::to_string(maximum_time - extant[0].get_birth_time());
		if(atof(s_BL.c_str()) < 0) {
		  Rcpp::Rcout << "extant.size = 1 bl < 0\n";
		}
		std::string output = s_ID;
		return output;
	}
	if(extant.size()==2)  {
		std::string s_ID1 = std::to_string(extant[0].ID);
		std::string s_BL1 = std::to_string(maximum_time - extant[0].get_birth_time());
		std::string s_ID2 =  std::to_string(extant[1].ID);
		std::string s_BL2 = std::to_string(maximum_time - extant[1].get_birth_time());

		if(atof(s_BL1.c_str()) < 0) {
		  Rcpp::Rcout << "extant.size == 2 bl1 < 0\n";
		}

		if(atof(s_BL2.c_str()) < 0) {
		  Rcpp::Rcout << "extant.size == 2 bl2 < 0\n";
		}

		std::string output = "(" + s_ID1 + ":" + s_BL1 + "," + s_ID2 + ":" + s_BL2 + ")";
		return output;
	}

	std::vector< newick_node > node_list = generateNodeList(v, maximum_time);

	//std::string core = node_list.back().composeString(node_list);
	Rcout << "creating newick string\n";
	std::string core = create_newick_string_local(node_list);
	return core;
}

std::vector<spec_point>  calculateLineages_noextinct(const std::vector<species>& allSp) {
	std::vector<spec_point> output;
	std::vector<species> allSpecies = allSp;

	std::sort(allSpecies.begin(), allSpecies.end(), sortOnTime);

	for(auto it = allSpecies.begin(); it != allSpecies.end(); ++it) {
		double time = (*it).get_birth_time();
		double previously_checked_time = -10;
		if(!output.empty()) previously_checked_time = output.back().get_time();

		if(time != previously_checked_time) {
			int lins = countLin(allSpecies,time);
			output.push_back(spec_point(lins, time));
		}
	}

	purgeOutput(output);

	return output;
}

bool species::check_has_viable_offspring(std::vector<species>& v) {
	if(checked == true) return extant_offspring;

	std::vector<int> offspring = find_indices(v, ID); //find the positions of the offspring;
	extant_offspring = false;

	for(std::size_t i = 0; i < offspring.size(); ++i) { //ofspring is of size 2 (or 1), so using iterators is useless here
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

std::vector<spec_point>  calculateLineages_withextinct(std::vector<species>& allSp,
                                                       double maximum_time)  {
	std::vector<spec_point> output;
	std::vector<species> lineages;

	for(auto it = allSp.begin(); it != allSp.end(); ++it)   {
	  if((*it).death_time == -1) lineages.push_back((*it)); //the species is an extant species and is added to the list to create the ltt plot
		else
		{
		  (*it).check_has_viable_offspring(allSp);
			if((*it).extant_offspring == true) lineages.push_back((*it));
		}
	}
	allSp = lineages;
	std::vector<species> filler;
	output = calculateLineages_noextinct(allSp);
	output.push_back(spec_point(output.back().ID, maximum_time));

	return output;
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

std::vector<spec_point> sumBranches(const std::vector<spec_point>& b1,
                                    const std::vector<spec_point>& b2,
                                    double maximum_time)
{
	std::vector<spec_point> output;
	std::vector<double> times;
	for(auto it = b1.begin(); it != b1.end(); ++it) times.push_back((*it).get_time());
	for(auto it = b2.begin(); it != b2.end(); ++it) times.push_back((*it).get_time());

	removeDuplicates(times); //to remove duplicates, the list also gets sorted, how convenient!

	int L1, L2;

	for(auto it = times.begin(); it != times.end(); ++it) {
		L1 = countLineages(b1, (*it));
		L2 = countLineages(b2, (*it));

		spec_point add(L1+L2, (*it));
		output.push_back(add);
	}

	output.push_back(spec_point(output.back().ID, maximum_time));
	removeDuplicates(output);
	return output;
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
  //std::cout << "looking for root\n";
  species root;
  for(auto it = v.begin(); it != v.end(); ++it) {
    if((*it).parent == -1) {
      root = (*it);
   //   std::cout << "found root, let's dive\n";
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


std::string create_newick_string( const std::vector<species>& s1,
                                  const std::vector<species>& s2,
                                  const std::vector<spec_point>& b1,
                                  const std::vector<spec_point>& b2,
                                  double maximum_time)
{
	double b_1, b_2;
	if(b1.size() == 1) b_1 = maximum_time - b1[0].get_time();
	else b_1 = b1[1].get_time();

	if(b2.size() == 1) b_2 = maximum_time - b2[0].get_time();
	else b_2 = b2[1].get_time();

	if(b_1 < 0) {
	  Rcpp::Rcout << "b_1 < 0!\n";
	}
	if(b_2 < 0) {
	  Rcpp::Rcout << "b_2 < 0!\n";
	}

	std::string left = writeTREE2(s1, maximum_time);
	std::string right = writeTREE2(s2, maximum_time);
	std::string BL_left = std::to_string(b_1);
	std::string BL_right = std::to_string(b_2);

	std::string output = "(";
    output += left;
	  output += ":";
  	output += BL_left;
  	output += ",";
  	output += right;
  	output += ":";
  	output += BL_right;
  	output += ");";

	return output;
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

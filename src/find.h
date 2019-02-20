//
//  find.hpp
//
//
//  Created by Thijs Janzen on 20/02/2019.
//

#ifndef find_hpp
#define find_hpp

#include "species.h"
#include "allo_pair.h"
#include "newick_node.h"
#include <vector>




int find_Allo(int i, const std::vector<species>& v);
int findinP(const std::vector<allo_pair>& p, int ID);

std::vector<int> find_indices(const std::vector<species>& v, int ID);

int find_parent(const std::vector<species>& v,
                int parent_id);


int findOther( int youngest,
              const std::vector<species>& v);

int findYoungest(const std::vector<species>& v);

std::vector<int> findOffspring(int ID,
                               const std::vector<species>& v);
std::vector<int> findOffspring(int ID,
                               const std::vector<newick_node>& v);



#endif /* find_hpp */

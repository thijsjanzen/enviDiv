//
//  newick_node.hpp
//
//
//  Created by Thijs Janzen on 20/02/2019.
//

#ifndef newick_node_hpp
#define newick_node_hpp

#include "species.h"
#include <string>
#include <vector>

struct newick_node
{
    newick_node(const species& S,
                double maximum_time);
    newick_node();

    bool extant;
    int ID;
    int parent;
    double branch_length;

    bool checked;
    int left;
    int right;

    std::string composeString(const std::vector<newick_node>& v) const;
};

#endif /* newick_node_hpp */

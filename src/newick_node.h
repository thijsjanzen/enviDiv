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

    bool checked;
    int left;
    int right;

    std::string composeString(const std::vector<newick_node>& v) const;
    double get_branch_length() const;
    void   set_branch_length(double bl);
  private:
    double branch_length;
};

#endif /* newick_node_hpp */

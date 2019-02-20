//
//  newick_node.cpp
//
//
//  Created by Thijs Janzen on 20/02/2019.
//

#include "newick_node.h"
#include "find.h"

newick_node::newick_node()  {
    ID = -1;
    branch_length = -1;
    extant = false;
    parent = -1;
}


std::string newick_node::composeString(const std::vector<newick_node>& v) const
{
    std::string output;
    static std::string haakjeOpen("(");
    static std::string haakjeSluit(")");
    static std::string komma(",");
    static std::string dubbelp(":");

    std::string s_ID = std::to_string(ID);
    std::string s_BL = std::to_string(branch_length);

    if(extant)  {
        output += s_ID;
        output += dubbelp;
        output += s_BL;
    } else {
        std::vector<int> children = findOffspring(ID,v);

        output += haakjeOpen;
        for(std::size_t i = 0; i < children.size(); ++i)  {
            output += v[children[i]].composeString(v);
            if(i != children.size()-1) output += komma;
        }
        output += haakjeSluit;

        if(parent != -1)  {
            output += dubbelp;
            output += s_BL;
        }
    }
    return output;
}

newick_node::newick_node(const species& S,
                         double maximum_time)
{
    if(S.death_time == -1)  {
        extant = true;
        branch_length = maximum_time - S.birth_time;
    } else {
        extant = false;
        branch_length = S.death_time - S.birth_time;
    }

    parent = S.parent;
    ID = S.ID;
}

//
//  find.cpp
//
//
//  Created by Thijs Janzen on 20/02/2019.
//

#include "find.h"
#include "newick_node.h"

// [[Rcpp::plugins(cpp11)]]

int find_Allo(int i, const std::vector<species>& v) {
    int ID = v[i].ID;
    for(std::size_t j = 0; j < v.size(); ++j)   {
        if(v[j].ID == ID && i != j) return j;
    }
    return -1;
}

std::vector<int> find_indices(const std::vector<species>& v, int ID)  {
    std::vector<int> indices;
    int count = 0;

    for(auto it = v.begin(); it != v.end(); ++it) {
        if((*it).parent == ID) indices.push_back(count);

        count++;
    }

    return indices;
}

int find_parent(const std::vector<species>& v,
                int parent_id)    {
    int count = 0;
    for(auto it = v.begin(); it != v.end(); ++it) {
        if((*it).ID == parent_id) return count;
        count++;
    }
    return 0;
}



int findOther( int youngest,
              const std::vector<species>& v) {
    int parent = v[youngest].parent;
    for(std::size_t i = 0; i < v.size(); ++i) {
        if(youngest != i) {
            if(v[i].parent == parent) {
                if(v[i].birth_time == v[youngest].birth_time)
                    return i;
            }
        }
    }
    return -1;
}


int findYoungest(const std::vector<species>& v) {
    double min = -1;
    int index = -1;
    int count = 0;

    for(auto it = v.begin(); it != v.end(); ++it) {
        if((*it).birth_time > min)  {
            min = (*it).birth_time;
            index = count;
        }
        count++;
    }
    return index;
}

std::vector<int> findOffspring(int ID,
                               const std::vector<species>& v) {
    std::vector<int> output;
    int counter = 0;
    for(auto it = v.begin(); it != v.end(); ++it) {
        if((*it).parent == ID) output.push_back(counter);
        counter++;
    }

    return output;
}

std::vector<int> findOffspring(int ID,
                               const std::vector<newick_node>& v)
{
    std::vector<int> output;

    int counter = 0;
    for(auto it = v.begin(); it != v.end(); ++it) {
        if((*it).parent == ID) output.push_back(counter);
        counter++;
    }

    return output;
}

int findinP(const std::vector<allo_pair>& p,
            int ID) {
    int count = 0;
    for(auto it = p.begin(); it != p.end(); ++it)
    {
        if((*it).ID == ID) return count;
        count++;
    }

    return -1;
}

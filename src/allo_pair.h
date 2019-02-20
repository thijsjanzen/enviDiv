//
//  allo_pair.hpp
//  
//
//  Created by Thijs Janzen on 20/02/2019.
//

#ifndef allo_pair_hpp
#define allo_pair_hpp


struct allo_pair
{
    int ID;
    int index_a;
    int index_b;

    allo_pair(int id, int index_A);
    allo_pair(int id, int index_A, int index_B);
};


#endif /* allo_pair_hpp */

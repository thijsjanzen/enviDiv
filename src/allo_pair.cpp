//
//  allo_pair.cpp
//  
//
//  Created by Thijs Janzen on 20/02/2019.
//

#include "allo_pair.h"

allo_pair::allo_pair(int id, int index_A) : ID(id),
                                            index_a(index_A)
{
    index_b = -1;
}

allo_pair::allo_pair(int id, int index_A, int index_B) : ID(id),
                                                         index_a(index_A),
                                                         index_b(index_B)
{
}

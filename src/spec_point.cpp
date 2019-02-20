//
//  struct_point.cpp
//  
//
//  Created by Thijs Janzen on 20/02/2019.
//

#include "spec_point.h"

spec_point::spec_point(double id, double t): ID(id), time(t)
{
}

bool spec_point::operator==(const spec_point& other) const
{
    if(time == other.time) return true;
    else return false;
}

bool spec_point::operator!=(const spec_point& other) const {return !((*this) == other);}
bool spec_point::operator<(const spec_point& other) const { return time < other.time; }

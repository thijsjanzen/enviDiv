//
//  struct_point.cpp
//
//
//  Created by Thijs Janzen on 20/02/2019.
//

#include "spec_point.h"
#include "assert.h"
#include "Rcpp.h"

spec_point::spec_point(double id, double t): ID(id)
{
  set_time(t);
}

bool spec_point::operator==(const spec_point& other) const
{
    if(time == other.time) return true;
    else return false;
}

bool spec_point::operator!=(const spec_point& other) const {return !((*this) == other);}
bool spec_point::operator<(const spec_point& other) const { return time < other.time; }

void spec_point::set_time(double t) {
  assert(t >= 0);
  if(t < 0) {
    Rcpp::Rcout << "spec_point sets negative BL\n";
  }
  time = t;
}

double spec_point::get_time() const {
  return time;
}


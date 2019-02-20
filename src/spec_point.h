//
//  struct_point.hpp
//  
//
//  Created by Thijs Janzen on 20/02/2019.
//

#ifndef struct_point_hpp
#define struct_point_hpp


struct spec_point
{
    spec_point(double id, double t);

    double ID;
    double time;

    bool operator<(const spec_point& other) const;

    bool operator==(const spec_point& other) const;
    bool operator!=(const spec_point& other) const;
};


#endif /* struct_point_hpp */

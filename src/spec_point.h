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

    double get_time() const;
    void set_time(double t);

    bool operator<(const spec_point& other) const;

    bool operator==(const spec_point& other) const;
    bool operator!=(const spec_point& other) const;
  private:
    double time;
};


#endif /* struct_point_hpp */

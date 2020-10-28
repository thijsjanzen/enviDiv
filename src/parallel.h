#ifndef parallel_h
#define parallel_h

#include <vector>

std::string do_run_tbb(const std::vector<float>& parameters,
                       const std::vector<float>& waterlevel_changes,
                       float maximum_time,
                       int max_lin,
                       std::vector< std::vector< float > >& l_table,
                       rnd_t& rndgen);

#endif parallel_h  /* parallel */

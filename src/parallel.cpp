#include <mutex>
#include <atomic>
#include <tbb/tbb.h>

#include <cmath>
#include <string>

#include "random_thijs.h"
#include "util.h"
#include "Gillespie.h"

#include <Rcpp.h>
using namespace Rcpp;

void force_output(std::string s) {
  Rcout << s << "\n";
  R_FlushConsole();
  R_ProcessEvents();
  R_CheckUserInterrupt();
}


std::string do_run_tbb(const std::vector<float>& parameters,
                       const std::vector<float>& waterlevel_changes,
                       float maximum_time,
                       int max_lin,
                       std::vector< std::vector< float > >& l_table,
                       rnd_t& rndgen);

//' function to test table conversion
//' @param model a vector of the four parameters of the model
//' @param crown_age crown age of the tree to be simulated
//' @return list
//' @export
// [[Rcpp::export]]
Rcpp::List test_envidiv_tbb(int model,
                            float crown_age) {

  rnd_t reng;

  std::vector<float> parameters = parameters_from_prior(reng);

  parameters[0] = 0.0;
  parameters[1] = 0.3;
  parameters[5] = model;

  std::vector<float> waterlevel_changes = get_waterlevel_changes(parameters[5],
                                                                 crown_age,
                                                                 reng);

  int max_lin = 1e6;
  std::vector< std::vector< float > > l_table;

  std::string code = do_run_tbb(parameters,
                                waterlevel_changes,
                                crown_age,
                                max_lin,
                                l_table,
                                reng);

  force_output("done simulation, starting on table conversion");

  std::string newick_string = ltable_to_newick(l_table, crown_age);

  force_output("done conversion");
  //  std::string newick_string = "placeholder";
  NumericMatrix ltab;
  if(!l_table.empty()) {
    ltab = NumericMatrix(l_table.size(), l_table[0].size());
    for(int i = 0; i < l_table.size(); ++i) {
      for(int j = 0; j < l_table[i].size(); ++j) {
        ltab(i, j) = l_table[i][j];
      }
    }
  }

  return List::create( Named("newick_string") = newick_string,
                       Named("Ltable") = ltab);
}


//' simulate a tree using environmental diversification
//' @param model a vector of the four paramaters of the model
//' @param num_repl a vector that indicates the time points of water level changes
//' @param crown_age crown age of the tree to be simulated
//' @param min_lin minimum number of lineages
//' @param max_lin maximum number of lineages
//' @param num_threads number of threads
//' @return newick string
//' @export
// [[Rcpp::export]]
List create_ref_table_tbb_serial(int model,
                                 int num_repl,
                                 float crown_age,
                                 int min_lin,
                                 int max_lin,
                                 int num_threads) {

  force_output("welcome");
  std::vector< std::string > trees;

  auto T0 = std::chrono::high_resolution_clock::now();
  int loop_size = num_repl - trees.size();

  while(trees.size() < num_repl) {
    // loop size can be optimized further, depending on the average success rate
    // e.g. loop_size = loop_size * 1.0f / success_rate
    // this is especially interesting once only a few are left.
    //Rcout << loop_size << "\n";
    loop_size = num_repl - trees.size();

   // force_output("loop_size: " + std::to_string(loop_size));

    std::vector< std::string > add(loop_size);
    std::vector< bool > add_flag(loop_size, false);

     for(int i = 0; i < loop_size; ++i) {
     rnd_t reng; // = rnd_t( make_random_engine<std::mt19937>() );

      std::vector<float> parameters = parameters_from_prior(reng);
      std::vector<float> waterlevel_changes = get_waterlevel_changes(parameters[5],
                                                                     crown_age,
                                                                     reng);

      std::vector< std::vector< float > > l_table;

      std::string code = do_run_tbb(parameters,
                                    waterlevel_changes,
                                    crown_age,
                                    max_lin,
                                    l_table,
                                    reng);

      int num_lin = 0;
      for(int k = 0; k < l_table.size(); ++k) {
        if(l_table[k][3] == -1) num_lin++;
      }

      if(num_lin >= min_lin && num_lin <= max_lin) {
        add[i]      = ltable_to_newick(l_table, crown_age);
        add_flag[i] = true;
      }
    }
     for(int j = 0; j < add_flag.size(); ++j) {
      if(add_flag[j]) {
        trees.push_back(add[j]);
      }
    }
  }

  Rcpp::List output(num_repl);

  for(int k = 0; k < trees.size(); ++k) {
    output[k] = trees[k];
  }

  auto T1 = std::chrono::high_resolution_clock::now();
  auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
  Rcout << "computed in: " << elapsed << "ms";
  return output;
}



//' simulate a tree using environmental diversification
//' @param model a vector of the four paramaters of the model
//' @param num_repl a vector that indicates the time points of water level changes
//' @param crown_age crown age of the tree to be simulated
//' @param min_lin minimum number of lineages
//' @param max_lin maximum number of lineages
//' @param num_threads number of threads
//' @return newick string
//' @export
// [[Rcpp::export]]
List create_ref_table_tbb_par(int model,
                              int num_repl,
                              float crown_age,
                              int min_lin,
                              int max_lin,
                              int num_threads) {

  std::vector< std::string > trees;
  std::vector< std::vector < float > > parameter_list;

  auto T0 = std::chrono::high_resolution_clock::now();
  int loop_size = num_repl - trees.size();

  float accept_rate = 1.0;
  int num_tried = 0;
  int num_accepted = 0;

  while(trees.size() < num_repl) {
    // loop size can be optimized further, depending on the average success rate
    // e.g. loop_size = loop_size * 1.0f / success_rate
    // this is especially interesting once only a few are left.
    loop_size = num_repl - trees.size();
    Rcout << "trees remaining: " << loop_size   <<
             " accept rate: "    << accept_rate << "\n";

    if (loop_size < 10) loop_size = 10;
    if (accept_rate > 0) loop_size *= 1.0 / accept_rate;
    if (loop_size > 1e6) loop_size = 1e6;


    std::vector< std::string > add(loop_size);
    std::vector< float > temp_filler(6);
    std::vector< std::vector< float > > add_params(loop_size, temp_filler);
    std::vector< bool > add_flag(loop_size, false);

    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);


    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, loop_size - 1),
      [&](const tbb::blocked_range<unsigned>& r) {
        rnd_t reng;

        for (unsigned i = r.begin(); i < r.end(); ++i) {
          std::vector<float> parameters = parameters_from_prior(reng, model);
          std::vector<float> waterlevel_changes = get_waterlevel_changes(parameters[5],
                                                                         crown_age,
                                                                         reng);

          std::vector< std::vector< float > > l_table;

          std::string code = do_run_tbb(parameters,
                                        waterlevel_changes,
                                        crown_age,
                                        max_lin,
                                        l_table,
                                        reng);

          int num_lin = 0;
          for(int k = 0; k < l_table.size(); ++k) {
            if(l_table[k][3] == -1) num_lin++;
          }

          if(num_lin >= min_lin && num_lin <= max_lin) {
            add[i]        = ltable_to_newick(l_table, crown_age);
            add_params[i] = parameters;
            add_flag[i]   = true;
          }
        }
      });

    for(int j = 0; j < add_flag.size(); ++j) {
      if(add_flag[j]) {
        trees.push_back(add[j]);
        parameter_list.push_back(add_params[j]);
      }
    }

    num_accepted  =  trees.size();
    num_tried    += loop_size;
    accept_rate = 1.0f * num_accepted / num_tried;
  }

  force_output("now converting to R objects\n");

  Rcpp::List output(trees.size());
  for(int k = 0; k < trees.size(); ++k) {
    output[k] = trees[k];
  }

  Rcpp::NumericMatrix parameter_matrix(parameter_list.size(), 6);
  for(int k = 0; k < parameter_list.size(); ++k) {
    for(int j = 0; j < parameter_list[0].size(); ++j) {
      parameter_matrix(k, j) = parameter_list[k][j];
    }
  }


  auto T1 = std::chrono::high_resolution_clock::now();
  auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::seconds>(T1 - T0).count());
  Rcout << "trees simulated in: " << elapsed << "seconds\n";
  return List::create(Named("trees") = trees,
                      Named("parameters") = parameter_matrix);
}


std::vector< int > get_parent_from_linlist(
    const std::vector< ltable_entry >& v,
    int parent) {
  std::vector< int > output;
  for(int i = 0; i < v.size(); ++i) {
    if(v[i].daughter == parent) {
      output.push_back(i);
    }
  }
  return(output);
}


std::string ltable_to_newick(const std::vector< std::vector< float > >& ltable,
                             float crown_age) {

  std::vector< ltable_entry > ltab(ltable.size());
  for(int i = 0; i < ltable.size(); ++i) {
    ltable_entry temp(ltable[i][0],
                      ltable[i][1],
                               ltable[i][2],
                                        ltable[i][3]);
    ltab[i] = temp;
  }

  for(auto& i : ltab) {
    i.bt = crown_age - i.bt;
    if(i.bt < 0) i.bt = 0;
  }
  // sort ltable on column 3 by abs value
  std::sort(ltab.begin(), ltab.end(),
            [](ltable_entry const& a,
               ltable_entry const& b) {return abs(a.daughter) < abs(b.daughter);});

  ltab[0].parent = 0;

  std::vector< ltable_entry > temp;
  for(auto i : ltab) {
    if(i.bt == crown_age) {
      temp.push_back(i);
    }
  }
  if(temp.size() > 1) { // this should always be true!
    bool connected = false;
    if (temp[1].daughter == temp[0].parent)
      connected = true;

    if( temp[0].daughter == temp[1].parent)
      connected = true;

    if (connected == false) {
      int parent_id = ltab[0].daughter;
      for(auto& i : ltab) {
        if(i.parent == -1) {
          i.parent = parent_id;
        }
      }
    }
  }


  std::vector< ltable_entry > linlist;
  float age = ltab[0].bt;
  float tend = age;
  // L[, 1] <- age - L[, 1]
  for(auto& i : ltab) {
    i.bt = age - i.bt;
  }
  ltab[0].bt = -1;

  for(auto& i : ltab) {
    if(i.extant == -1) {
      i.tend = tend;
      i.label = "t" + std::to_string((int)abs(i.daughter));
      linlist.push_back(i);
    }
  }

  bool done = false;
  while(!done) {
    auto m = std::max_element(linlist.begin(), linlist.end(),
                              [](const ltable_entry& a, const ltable_entry& b) {
                                return a.bt < b.bt;
                              });

    int j = std::distance(linlist.begin(), m);

    if (j > linlist.size() || j < 0) {
      return("j out of bounds");
    }

    int parent   = linlist[j].parent;
    std::vector<int> parentj = get_parent_from_linlist(linlist, parent);

    if(parentj.size() == 1) {
      int index = parentj[0];
      if(index > linlist.size() || index < 0) {
        return("index out of bounds");
      }
      std::string spec1 = linlist[index].label + ":" +
        std::to_string(linlist[index].tend - linlist[j].bt);

      std::string spec2 = linlist[j].label + ":" +
        std::to_string(linlist[j].tend - linlist[j].bt);

      linlist[index].label = "(" + spec1 + "," + spec2 + ")";
      linlist[index].tend = linlist[j].bt;
      linlist[j] = linlist.back();
      linlist.pop_back();

    } else {
      std::vector< int > indices = get_parent_from_linlist(ltab, parent);
      if(!indices.empty()) {
        int parent_index = indices[0];

        if(parent_index > ltab.size() || parent_index < 0) {
          return("parent_index out of bounds");
        }

        linlist[j].bt        = ltab[parent_index].bt;
        linlist[j].parent    = ltab[parent_index].parent;
        linlist[j].daughter  = ltab[parent_index].daughter;
      } else {

      }
    }

    if(linlist.size() == 1) {
      done = true;
    }
  }

  std::string newick_tree = linlist[0].label + ":" +
    std::to_string(linlist[0].tend) + ";";

  return newick_tree;
}

//' test multithreaded works!
//' @param n size of vector
//' @param num_threads number of threads
//' @return vector of squared numbers
//' @export
// [[Rcpp::export]]
NumericVector sq_numbers_cpp_tbb(int n,
                                 int num_threads) {
  NumericVector results(n);
  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);

  tbb::parallel_for(
    tbb::blocked_range<unsigned>(0, n),
    [&](const tbb::blocked_range<unsigned>& r) {
      for (unsigned i = r.begin(); i < r.end(); ++i) {
        results[i] = i * sin(i) + cos(i);
      }
    });

  return results;
}


std::vector< std::vector< float >> create_l_table_float(
    const std::vector< species > & s1,
    const std::vector< species > & s2) {

  std::vector< std::vector< float > > output;
  int num_rows = s1.size() + s2.size();
  int num_cols = 4;
  for(int i = 0; i < num_rows; ++i) {
    std::vector< float > row(num_cols);
    output.push_back(row);
  }

  int i = 0;
  for (auto it = s1.begin(); it != s1.end(); ++it, ++i) {
    output[i][0] = (*it).birth_time;
    output[i][1] = (*it).get_parent();
    output[i][2] = (*it).get_ID();
    output[i][3] = (*it).death_time;
  }

  for (auto it = s2.begin(); it != s2.end(); ++it, ++i) {
    output[i][0]  = (*it).birth_time;
    output[i][1]  = -1 * ((*it).get_parent());
    output[i][2]  = -1 * ((*it).get_ID());
    output[i][3]  = (*it).death_time;
  }

  return output;
}




std::string do_run_tbb(const std::vector<float>& parameters,
                       const std::vector<float>& waterlevel_changes,
                       float maximum_time,
                       int max_lin,
                       std::vector< std::vector< float > >& l_table,
                       rnd_t& rndgen)
{
  std::vector< species > s1;

  float jiggle_amount = parameters[4];
  rndgen.set_normal_trunc(0.0f, jiggle_amount);

  //int idCount = 0;
  std::vector < std::vector < float > > l_table1;
  int error_code = run(parameters, waterlevel_changes,
                       s1, maximum_time, max_lin,
                       rndgen);

  if (error_code == 0) {
    return "extinction";
  }

  if (error_code == 1e6) {
    return "overflow";
  }

  std::vector<species> s2;

  int error_code2 = run(parameters,
                        waterlevel_changes,
                        s2, maximum_time, max_lin,
                        rndgen);

  if (error_code2 == 0) {
    return "extinction";
  }
  if (error_code2 == 1e6) {
    return "overflow";
  }

  jiggle(s1, s2, maximum_time, jiggle_amount, rndgen);

  l_table = create_l_table_float(s1, s2);

  std::string output = "success";
  return output;
}

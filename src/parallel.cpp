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
List create_ref_table_tbb(int model,
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

    force_output("loop_size: " + std::to_string(loop_size));

    std::vector< std::string > add(loop_size);
    std::vector< bool > add_flag(loop_size, false);

    //  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);


    for(int i = 0; i < loop_size; ++i) {
  //  tbb::parallel_for(
//      tbb::blocked_range<unsigned>(0, loop_size),
//      [&](const tbb::blocked_range<unsigned>& r) {
        rnd_t thread_local reng = rnd_t( make_random_engine<std::mt19937>() );

  //      for (unsigned i = r.begin(); i < r.end(); ++i) {
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

          force_output("done run");
          int num_lin = 0;
          for(int k = 0; k < l_table.size(); ++k) {
            if(l_table[k][3] == -1) num_lin++;
          }

          if(num_lin >= min_lin && num_lin <= max_lin) {
            force_output(std::to_string(i));

            add[i]      = ltable_to_newick(l_table, crown_age);
            add_flag[i] = true;
            force_output("created newick and added to vector");
          }
        }
 //     });

      force_output("done loop\n");
      for(int j = 0; j < add_flag.size(); ++j) {
        if(add_flag[j]) {
          force_output(std::to_string(j));
          trees.push_back(add[j]);
        }
      }
      force_output("added trees");
  }

  force_output("tbb done\n");

  Rcpp::List output(num_repl);

  for(int k = 0; k < trees.size(); ++k) {
    output[k] = trees[k];
  }

  auto T1 = std::chrono::high_resolution_clock::now();
  auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
  std::cout << "computed in: " << elapsed << "ms";
  return output;
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

  force_output("welcome to ltable_to_newick");
  std::vector< ltable_entry > ltab(ltable.size());
  for(int i = 0; i < ltable.size(); ++i) {
    ltable_entry temp(
                      ltable[i][0],
                      ltable[i][1],
                      ltable[i][2],
                      ltable[i][3]);
    ltab[i] = temp;
  }

  /*
   local_l_table <- input_matrix
   local_l_table[, 1] <- crown_age - local_l_table[, 1]
  */
  for(auto& i : ltab) {
    i.bt = crown_age - i.bt;
    if(i.bt < 0) i.bt = 0;
  }
  // sort ltable on column 3 by abs value
  std::sort(ltab.begin(), ltab.end(),
            [](ltable_entry const& a,
               ltable_entry const& b) {return abs(a.daughter) < abs(b.daughter);});

  ltab[0].parent = 0;
  /*
  * local_l_table <-  local_l_table[order(abs(local_l_table[, 3])), 1:4]

  * local_l_table[1, 2] <- 0
  * local_l_table[which(local_l_table[, 1] < 0), 1] <- 0

  */
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

  /*
   a <- subset(local_l_table, local_l_table[, 1] == crown_age)
   connected <- FALSE
   if (a[2, 3] == a[1, 2]) connected <- TRUE
   if (a[1, 3] == a[2, 2]) connected <- TRUE

   if (connected == FALSE) {
   parent_id <- local_l_table[1, 3]
   local_l_table[which(local_l_table[, 2] == -1), 2] <- parent_id
   }
   */


  force_output("ltab sorted");

  std::vector< ltable_entry > linlist;
  float age = ltab[0].bt;
  float tend = age;
  for(auto& i : ltab) {
    i.bt = age - i.bt;
    if(i.extant == -1) {
      i.tend = tend;
      i.label = "t" + std::to_string((int)abs(i.daughter));
      linlist.push_back(i);
    }
  }
  linlist[0].bt = -1;

  force_output("transported to linlist");

  bool done = false;
  while(!done) {
    force_output(std::to_string(linlist.size()));
    auto m = std::max_element(linlist.begin(), linlist.end(),
                              [](const ltable_entry& a, const ltable_entry& b) {
                                return a.bt < b.bt;
                              });

    int j = std::distance(linlist.begin(), m); //get_max(linlist);
    int parent   = linlist[j].parent;
    std::vector<int> parentj = get_parent_from_linlist(linlist, parent);

    if(parentj.size() == 1) {
      int index = parentj[0];
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
       int parent_index = indices[0];
       linlist[j].bt        = ltab[parent_index].bt;
       linlist[j].parent    = ltab[parent_index].bt;
       linlist[j].daughter  = ltab[parent_index].bt;
    }

    if(linlist.size() == 1) {
      done = true;
    }
  }

  std::string newick_tree = linlist[0].label + ":" +
    std::to_string(linlist[0].tend) + ";";

  return newick_tree;
}


//' test ltable
//' @param input_table input table
//' @export
// [[Rcpp::export]]
std::string test_ltable_to_newick(const NumericMatrix& input_matrix) {
  std::vector< std::vector< float > > ltable;
  for(int i = 0; i < input_matrix.nrow(); ++i) {
    NumericVector v = input_matrix(i, _);
    std::vector< float > temp;
    for(int j = 0; j < v.size(); ++j) {
      temp.push_back(v[j]);
    }
    ltable.push_back(temp);
  }

  // std::string newick = ltable_to_newick(ltable);
  std::string newick = "test";
  return newick;
}

int get_num_lin(const NumericMatrix& l_table) {
  int cnt = 0;
  for(int i = 0; i < l_table.nrow(); ++i) {
    if(l_table(i, 3) == -1) cnt++;
  }
  return(cnt);
}

//' simulate a tree using environmental diversification
//' @param model a vector of the four paramaters of the model
//' @param num_repl a vector that indicates the time points of water level changes
//' @param crown_age crown age of the tree to be simulated
//' @param min_lin minimum number of lineages
//' @param max_lin maximum number of lineages
//' @return newick string
//' @export
// [[Rcpp::export]]
List create_ref_table_cpp_serial(int model,
                                 int num_repl,
                                 float crown_age,
                                 int min_lin,
                                 int max_lin) {

  int num_remaining = num_repl;

  Rcpp::List trees(num_remaining);

  Environment pkg = Environment::namespace_env("enviDiv");
  Function to_newick = pkg["sim_table_to_newick"];

  int cnt = 0;
  auto T0 = std::chrono::high_resolution_clock::now();

  while(cnt < num_remaining) {
    int loop_size = num_remaining - cnt;
    Rcout << loop_size << "\n";

    for (unsigned i = 0; i < loop_size; ++i) {
      std::vector<float> parameters = param_from_prior_cpp();
      int water_model = parameters[5];

      std::vector<float> waterlevel_changes = get_waterlevel_cpp(water_model,
                                                                 crown_age);

      rnd_t thread_local rndgen;
      rndgen.set_seed(cnt * i);

      NumericMatrix l_table;
      std::string code = do_run_r(parameters,
                                  waterlevel_changes,
                                  crown_age,
                                  max_lin,
                                  l_table,
                                  rndgen);

      //  Rcout << code << "\n";
      int num_lin = get_num_lin(l_table);
      if(num_lin >= min_lin && num_lin <= max_lin) {
        auto input = to_newick(l_table, crown_age);
        trees[cnt] = input;
        ++cnt;
      }
    }
  }
  auto T1 = std::chrono::high_resolution_clock::now();
  auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
  std::cout << "computed in: " << elapsed << "ms";

  return trees;
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


//' test for tbb implementation
//' @param loop_size size of task
//' @param num_threads number of threads
//' @return list
//' @export
// [[Rcpp::export]]
List test_tbb(  int loop_size,
                int num_threads) {

  auto T0 = std::chrono::high_resolution_clock::now();

  std::vector< float > for_output;

  tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
  std::mutex mutex;

  tbb::parallel_for(
    tbb::blocked_range<unsigned>(0, loop_size, 1),
    [&](const tbb::blocked_range<unsigned>& r) {
      for (unsigned i = r.begin(); i < r.end(); ++i) {
        try {

          rnd_t thread_local reng = rnd_t( make_random_engine<std::mt19937>() );

          float a = reng.Expon(0.1);
          float b = reng.Expon(0.5);
          float add = sin(a) + cos(b);

          std::lock_guard<std::mutex> _(mutex);
          for_output.push_back(add);
        }
        catch(const std::exception& e) {
          Rcout << "runtime error\n";
        }
      }
    });



  Rcpp::List output(loop_size);
  for(int i = 0; i < for_output.size(); ++i) {
    output[i] = for_output[i];
  }
  auto T1 = std::chrono::high_resolution_clock::now();
  auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
  std::cout << "computed in: " << elapsed << "ms";

  return output;
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

/*


 int get_max(const std::vector< ltable_entry >& v) {
 int index = 0;
 int max = v[0].bt;
 int cnt = 0;
 for(const auto& i : v) {
 if( i.bt >= max) {
 max = i.bt;
 index = cnt;
 }
 cnt++;
 }
 return index;
 }


 void output_ltable(const ltable_entry& l) {
 Rcout << l.bt << " " << l.parent << " " << l.daughter << " ";
 Rcout << l.label << " " << l.tend << "\n";
 R_FlushConsole();
 R_ProcessEvents();
 R_CheckUserInterrupt();
 }

 void force_output(std::string s) {
 Rcout << s << "\n";
 R_FlushConsole();
 R_ProcessEvents();
 R_CheckUserInterrupt();
 }

 std::string ltable_to_newick(const std::vector< std::vector< float > > ltable) {

 std::vector< ltable_entry > ltab;
 for(int i = 0; i < ltable.size(); ++i) {
 ltable_entry temp;
 temp.bt       = ltable[i][0];
 temp.parent   = ltable[i][1];
 temp.daughter = ltable[i][2];
 temp.extant   = ltable[i][3];
 ltab.push_back(temp);
 }

 // sort ltable on column 3 by abs value
 std::sort(ltab.begin(), ltab.end(),
 [](ltable_entry const& a,
 ltable_entry const& b) {return abs(a.daughter) < abs(b.daughter);});

 float age = ltab[0].bt;
 for(auto& i : ltab) {
 i.bt = age - i.bt;
 if (i.extant != -1) {
 i.extant = age - i.extant;
 }
 }
 ltab[0].bt = -1;
 float tend = age;
 std::vector< ltable_entry > linlist;

  // L = L[order(abs(L[, 3])), 1:4]
  // age = L[1, 1]
  // L[, 1] = age - L[, 1]
  // L[1, 1] = -1
  // notmin1 = which(L[, 4] != -1)
  // L[notmin1, 4] = age - L[notmin1, 4]
  // if (dropextinct == T) {
  //  sall = which(L[, 4] == -1)
  //  tend = age
  // } else {
  //  sall = which(L[, 4] >= -1)
  //  tend = (L[, 4] == -1) * age + (L[, 4] > -1) * L[, 4]
  // }
  // L = L[, -4]

  // linlist = cbind(data.frame(L[sall, ]), paste("t", abs(L[sall,
  3]), sep = ""), tend)
  //  linlist[, 4] = as.character(linlist[, 4])
  //  names(linlist) = 1:5

 for(int i = 0;i < ltab.size(); ++i) {
   if(ltab[i].extant == -1) {
     ltab[i].tend = tend;
     ltab[i].label = "t" + std::to_string((int)abs(ltab[i].daughter));
     linlist.push_back(ltab[i]);
   }
 }

 ///// verified correct up until here!!!  ////
 /////////////////////////////////////////////

 bool done = false;
while(!done) {

  int j = get_max(linlist);
  int parent   = linlist[j].parent;
  std::vector<int> parentj = get_parent_from_linlist(linlist, parent);


  // done = 0
  // while (done == 0) {
  // j = which.max(linlist[, 1])
  // daughter = linlist[j, 3]
  // parent = linlist[j, 2]
  // parentj = which(parent == linlist[, 3])
  // parentinlist = length(parentj)

  if(parentj.size() == 1) {
    int index = parentj[0];
    std::string spec1 = linlist[index].label + ":" +
      std::to_string(linlist[index].tend - linlist[j].bt);

    std::string spec2 = linlist[j].label + ":" +
      std::to_string(linlist[j].tend - linlist[j].bt);

    linlist[index].label = "(" + spec1 + "," + spec2 + ")";
    linlist[index].tend = linlist[j].bt;
    linlist[j] = linlist.back();
    linlist.pop_back();

    // spec1 = paste(linlist[parentj, 4], ":", linlist[parentj,
    // 5] - linlist[j, 1], sep = "")
    // spec2 = paste(linlist[j, 4], ":", linlist[j, 5] -
     // linlist[j, 1], sep = "")
     // linlist[parentj, 4] = paste("(", spec1, ",", spec2,
     // ")", sep = "")
     // linlist[parentj, 5] = linlist[j, 1]
     // linlist = linlist[-j, ]

  } else {

    // linlist[j, 1:3] = L[which(L[, 3] == parent), 1:3]
    std::vector< int > indices = get_parent_from_linlist(ltab, parent);
    int parent_index = indices[0];
    linlist[j].bt = ltab[parent_index].bt;
    linlist[j].parent = ltab[parent_index].bt;
    linlist[j].daughter = ltab[parent_index].bt;
  }

  if(linlist.size() == 1) {
    done = true;
  }
}

 // linlist[4] = paste(linlist[4], ":", linlist[5], ";", sep = "")
 // phy = ape::read.tree(text = linlist[1, 4])
 // tree = ape::as.phylo(phy)
 // return(tree)
 //


std::string newick_tree = linlist[0].label + ":" +
  std::to_string(linlist[0].tend) + ";";

return newick_tree;
}

 std::vector< float > parameters = param_from_prior_cpp();
 std::vector< float > waterlevel_changes = get_waterlevel_cpp(parameters[5],
 crown_age);

 rnd_t thread_local rndgen(cnt * i);

 Rcout << i << "\n";
 NumericMatrix l_table;
 std::string code = do_run_r(parameters,
 waterlevel_changes,
 crown_age,
 max_lin,
 l_table,
 rndgen);

 Rcout << code << "\n";
 int num_lin = get_num_lin(l_table);
 if(num_lin >= min_lin && num_lin <= max_lin) {
 auto input = to_newick(l_table, crown_age);
 std::lock_guard<std::mutex> _(mutex);
 trees[cnt] = input;
 ++cnt;
 }

 */

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

  int num_remaining = num_repl;

  std::vector< Rcpp::NumericMatrix > l_tables;

  auto T0 = std::chrono::high_resolution_clock::now();
  int loop_size = num_remaining - l_tables.size();

  while(loop_size > 0) {
    Rcout << loop_size << "\n";

    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);
    std::mutex mutex;

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, loop_size, 1),
      [&](const tbb::blocked_range<unsigned>& r) {
        for (unsigned i = r.begin(); i < r.end(); ++i) {
          try {

            rnd_t thread_local reng = rnd_t( make_random_engine<std::mt19937>() );

            std::vector<float> parameters = parameters_from_prior(reng);
            std::vector<float> waterlevel_changes = get_waterlevel_changes(parameters[5],
                                                        crown_age,
                                                        reng);

            NumericMatrix l_table;

            std::string code = do_run_r(parameters,
                                        waterlevel_changes,
                                        crown_age,
                                        max_lin,
                                        l_table,
                                        reng);

            int num_lin = 0;
            for(int i = 0; i < l_table.nrow(); ++i) {
              if(l_table(i, 3) == -1) num_lin++;
            }

            if(num_lin >= min_lin && num_lin <= max_lin) {
              std::lock_guard<std::mutex> _(mutex);
              l_tables.push_back(l_table);
            }
          }
          catch(const std::exception& e) {
            Rcout << "runtime error\n";
          }
        }
      });
    loop_size = num_remaining - l_tables.size();
  }

  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("enviDiv");
  Rcpp::Function to_newick = pkg["sim_table_to_newick"];

  Rcpp::List trees(num_repl);

  for(int i = 0; i < l_tables.size(); ++i) {
    trees[i] = to_newick(l_tables[i], crown_age);
  }

  auto T1 = std::chrono::high_resolution_clock::now();
  auto elapsed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(T1 - T0).count());
  std::cout << "computed in: " << elapsed << "ms";
  return trees;
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

/*

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

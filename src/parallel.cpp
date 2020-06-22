#include <tbb/tbb.h>
#include "Gillespie.h"
#include "util.h"

#include <cmath>
#include "random_thijs.h"
#include <string>

#include <Rcpp.h>
using namespace Rcpp;


int get_num_lin(const NumericMatrix& l_table) {
  int cnt = 0;
  for(int i = 0; i < l_table.nrow(); ++i) {
    if(l_table(i, 3) == -1) cnt++;
  }
  return(cnt);
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
List create_ref_table_cpp(int model,
                          int num_repl,
                          float crown_age,
                          int min_lin,
                          int max_lin,
                          int num_threads) {

  // Obtaining namespace of Matrix package
  Environment pkg = Environment::namespace_env("enviDiv");
  Function get_prior = pkg["generate_from_prior"];
  Function get_waterlevel_changes = pkg["generate_water"];

  std::vector< NumericMatrix > l_tables;

  int num_remaining = num_repl;
  int cnt = 0;
  while(cnt < num_remaining) {
    int loop_size = num_remaining - cnt;
    tbb::task_scheduler_init _tbb((num_threads > 0) ? num_threads : tbb::task_scheduler_init::automatic);

    tbb::parallel_for(
      tbb::blocked_range<unsigned>(0, loop_size, 1),
      [&](const tbb::blocked_range<unsigned>& r) {

        for (unsigned i = r.begin(); i < r.end(); ++i) {
            NumericVector parameters = get_prior();
            int water_model = parameters[5];
            NumericVector waterlevel_changes = get_waterlevel_changes(water_model,
                                                                      crown_age);
            rnd_t rndgen;
            rndgen.set_seed(cnt * i);

            NumericMatrix l_table;
            std::string code = do_run_r(parameters,
                                        waterlevel_changes,
                                        crown_age,
                                        max_lin,
                                        l_table,
                                        rndgen);

            int num_lin = get_num_lin(l_table);
            if(num_lin >= min_lin && num_lin <= max_lin) {
              l_tables[cnt] = l_table;
              ++cnt;
            }
        }
      });
  }

  return List::create( Named("Ltable") = l_tables);
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
List create_ref_table_cpp_serial(int model,
                          int num_repl,
                          float crown_age,
                          int min_lin,
                          int max_lin,
                          int num_threads) {

  // test code
  NumericVector params = param_from_prior_cpp();
  for(int i = 0; i < params.size(); ++i) {
    Rcout << params[i] << " ";
  }
  Rcout << "\n";

  std::vector<float> w = get_waterlevel_cpp(1, 5);

  for(int i = 0; i < w.size(); ++i) {
    Rcout << w[i] << " ";
  }
  Rcout << "\n";

   w = get_waterlevel_cpp(2, 5);
   for(int i = 0; i < w.size(); ++i) {
     Rcout << w[i] << " ";
   }
   Rcout << "\n";

   w = get_waterlevel_cpp(3, 5);
   for(int i = 0; i < w.size(); ++i) {
     Rcout << w[i] << " ";
   }
   Rcout << "\n";

  std::vector< NumericMatrix > l_tables;
  /*





  int num_remaining = num_repl;
  int cnt = 0;
  while(cnt < num_remaining) {
    int loop_size = num_remaining - cnt;
    Rcout << "num_remaining\n";

    for (unsigned i = 0; i < loop_size; ++i) {
          NumericVector parameters = get_prior();
          int water_model = parameters[5];

          for(int j = 0; j < parameters.size(); ++j) {
            Rcout << parameters[j] << " ";
          }
          Rcout << "\n";

          NumericVector waterlevel_changes = get_waterlevel_changes(water_model,
                                                                    crown_age);
          rnd_t rndgen;
          rndgen.set_seed(cnt * i);

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
            l_tables[cnt] = l_table;
            ++cnt;
          }
    }
  }
*/
  return List::create( Named("Ltable") = l_tables);
}

#' function to perform ABC-SMC
#' @param number_of_replicates number of particles used per iteration of
#'                            the SMC algorithm
#' @param min_tips minimum number of tips
#' @param max_tips maximum number of tips
#' @param model used water model
#' @param emp_tree phy object holding phylogeny of the tree to be fitted on
#' @param crown_age crown age
#' @param write_to_file boolean, if TRUE, results are written to file.
#' @param file_name file name
#' @return a tibble containing the results
#' @export
generate_stack <- function(number_of_replicates = 1000,
                           parameters = NULL,
                           min_tips = 50,
                           max_tips = 150,
                           emp_tree = NULL,
                           crown_age = NULL,
                           write_to_file = FALSE,
                           file_name = NULL) {

  if(!is.null(emp_tree)) {
    crown_age <- max(ape::branching.times(emp_tree))
  }
  if(is.null(crown_age)) {
    stop("Please either provide a reference tree, or provide the crown age")
  }

  if(is.null(parameters)) {
    parameters <- enviDiv::param_from_prior()
  }

  number_accepted <- 0
  remaining_particles <- number_of_replicates - number_accepted

  all_results <- c()

  while (remaining_particles > 0) {
    cat(remaining_particles, "\n")
    sample_size <- max(1000, remaining_particles) #increase if not testing

    param_matrix <- matrix(parameters,
                           nrow = sample_size,
                           ncol = 7, byrow = T)  #6 parameters

    candidate_particles <- param_matrix


    calc_tree_stats <- function(x) {
      stats <- rep(Inf, 15)

      if(x[6] == 4) {
        found_tree <- TreeSim::sim.bd.age(age = crown_age, numbsim = 1, lambda = x[2], mu = x[1], complete = FALSE)[[1]]
      } else {
        found_tree <- enviDiv::sim_envidiv_tree(x, crown_age, abc = TRUE)

      }

      if (is.null(found_tree)) {
        return(stats)
      }

      num_tips <- found_tree$Nnode + 1

      if (num_tips >= min_tips && num_tips <= max_tips) {
        stats <- enviDiv::calc_sum_stats(found_tree, emp_tree)
      }
      return(stats)
    }

    input <- lapply(seq_len(nrow(candidate_particles)),
                    function(i) candidate_particles[i, ])

    stats <- future.apply::future_lapply(input, calc_tree_stats)

    stat_matrix <- matrix(unlist(stats, use.names = FALSE),
                          ncol = 15,
                          byrow = TRUE)

    results <- cbind(candidate_particles, stat_matrix)

    stat_matrix <- stat_matrix[!is.infinite(results[, 8]), ]
    results <- results[!is.infinite(results[, 8]), ]

    if (length(stat_matrix) > 0) {

      num_local_accepted <- nrow(results)
      if (!is.null(num_local_accepted)) {
        number_accepted <- number_accepted + num_local_accepted

        remaining_particles <- number_of_replicates - number_accepted
        colnames(results) <-
          c("extinct", "sym_high", "sym_low", "allo", "jiggle", "model",
            "weight",
            "nltt", "gamma", "mbr", "num_lin",
            "beta", "colless", "sackin", "ladder", "cherries", "ILnumber",
            "pitchforks", "stairs",
            "spectr_eigen", "spectr_asymmetry", "spectr_peakedness")
        results <- tibble::as_tibble(results)
        if (write_to_file) {
          readr::write_tsv(results, path = file_name,
                           append = TRUE)
        } else {
          all_results <- rbind(all_results, results)
        }
      }
    }
  }
  return(all_results)
}

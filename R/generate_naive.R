#' function to perform ABC-SMC
#' @param number_of_particles number of particles used per iteration of
#'                            the SMC algorithm
#' @param min_tips minimum number of tips
#' @param max_tips maximum number of tips
#' @param model used water model
#' @param emp_tree phy object holding phylogeny of the tree to be fitted on
#' @param write_to_file boolean, if TRUE, results are written to file.
#' @param file_name file name
#' @return a tibble containing the results
#' @export
generate_naive <- function(number_of_particles = 1000,
                           min_tips = 50,
                           max_tips = 150,
                           model = -1,
                           emp_tree = NULL,
                           crown_age = NULL,
                           write_to_file = FALSE,
                           file_name) {

  if(!is.null(emp_tree)) {
    crown_age <- max(ape::branching.times(emp_tree))
  }
  if(is.null(crown_age)) {
    stop("Please either provide a reference tree, or provide the crown age")
  }

  number_accepted <- 0
  remaining_particles <- number_of_particles - number_accepted

  all_results <- c()

  while (remaining_particles > 0) {
    cat(remaining_particles, "\n")
    sample_size <- min(10000, remaining_particles)

    param_matrix <- matrix(NA, nrow = sample_size,
                           ncol = 7)  #6 parameters

    candidate_particles <- t(apply(param_matrix, 1, enviDiv::param_from_prior))
    if (model != -1) candidate_particles[, 6] <- model

    calc_tree_stats <- function(x) {
      stats <- rep(Inf, 8)
      found_tree <- enviDiv::sim_envidiv_tree(x, crown_age, abc = TRUE)
      if (is.null(found_tree)) {
        return(stats)
      }

      num_tips <- found_tree$Nnode + 1

      if (num_tips >= min_tips && num_tips <= max_tips) {
        stats <- enviDiv::calc_sum_stats(found_tree, emp_tree)[1:8]
      }
      return(stats)
    }

    input <- lapply(seq_len(nrow(candidate_particles)),
                    function(i) candidate_particles[i, ])

    stats <- future.apply::future_lapply(input, calc_tree_stats)

    stat_matrix <- matrix(unlist(stats, use.names = FALSE),
                          ncol = 8,
                          byrow = TRUE)

    results <- cbind(candidate_particles, stat_matrix)

    stat_matrix <- stat_matrix[!is.infinite(results[, 8]), ]
    results <- results[!is.infinite(results[, 8]), ]

    if (length(stat_matrix) > 0) {

      num_local_accepted <- nrow(results)
      if (!is.null(num_local_accepted)) {
        number_accepted <- number_accepted + num_local_accepted

        remaining_particles <- number_of_particles - number_accepted
        colnames(results) <-
          c("extinct", "sym_high", "sym_low", "allo", "jiggle", "model",
            "weight",
            "nltt", "gamma", "mbr", "num_lin",
            "beta", "colless", "sackin", "ladder")
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

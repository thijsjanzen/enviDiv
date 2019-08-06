#' function to perform ABC-SMC
#' @param number_of_particles number of particles used per iteration of
#'                            the SMC algorithm
#' @param min_tips minimum number of tips
#' @param max_tips maximum number of tips
#' @param model used water model
#' @param emp_tree phy object holding phylogeny of the tree to be fitted on
#' @param file_name file name
#' @return a tibble containing the results
#' @export
generate_naive <- function(number_of_particles,
                           min_tips,
                           max_tips,
                           model,
                           emp_tree,
                           file_name) {

  #oplan <- future::plan()
  #on.exit(future::plan(oplan), add = TRUE)

  #future::plan(future::multiprocess)

  crown_age <- max(ape::branching.times(emp_tree))

  number_accepted <- 0
  remaining_particles <- number_of_particles - number_accepted

  while (remaining_particles > 0) {
    cat(remaining_particles, "\n")
    sample_size <- max(10000, remaining_particles)

    param_matrix <- matrix(NA, nrow = sample_size,
                           ncol = 7)  #6 parameters

    candidate_particles <- t(apply(param_matrix, 1, enviDiv::param_from_prior))
    candidate_particles[, 6] <- model

    calc_tree_stats <- function(x) {
      stats <- rep(Inf, 8)
      found_tree <- enviDiv::sim_envidiv_tree(x, crown_age, abc = TRUE)
      if(is.null(found_tree)) {
        return(stats)
      }

      num_tips = found_tree$Nnode + 1

      if(num_tips >= min_tips && num_tips <= max_tips) {
        stats <- enviDiv::calc_sum_stats(found_tree, emp_tree)[1:8]
      }
      return(stats)
    }

    stat_matrix <- t(future.apply::future_apply(candidate_particles,
                                              1,
                                              calc_tree_stats))

    #stat_matrix <- t(apply(candidate_particles,
    #                     1,
    #                     calc_tree_stats))


    results <- cbind(candidate_particles, stat_matrix)

    results <- results[!is.infinite(results[, 8]), ]

    number_accepted <- number_accepted + length(results[, 1])
    remaining_particles <- number_of_particles - number_accepted
    colnames(results) <-
      c("extinct", "sym_high", "sym_low", "allo", "jiggle", "model",
        "weight",
        "nltt", "gamma", "mbr", "num_lin",
        "beta", "colless", "sackin", "ladder")
    results <- tibble::as_tibble(results)
    readr::write_tsv(results, path = file_name,
                     append = TRUE)
  }

}

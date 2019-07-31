#' function to perform ABC-SMC
#' @param number_of_particles number of particles used per iteration of
#'                            the SMC algorithm
#' @param max_iter maximum number of iterations
#' @param sd_params standard deviation of the paramater perturbation kernel
#' @param emp_tree phy object holding phylogeny of the tree to be fitted on
#' @return nothing
#' @export
infer_params <- function(number_of_particles,
                         max_iter,
                         sd_params,
                         emp_tree) {

  param_matrix <- matrix(NA, nrow = number_of_particles,
                         ncol = 7)  #6 parameters

  emp_stats <- calc_sum_stats(emp_tree, emp_tree)
  crown_age <- max(ape::branching.times(emp_tree))

  # generate from prior:
  previous_par <- t(apply(param_matrix, 1, param_from_prior))

  # now we start SMC
  for (iter in 2:max_iter) {
    cat("iteration: ", iter, "\n")
    local_eps <- 300 * exp(-0.5 * iter)

    next_par <- c()

    remaining_particles <- number_of_particles - length(next_par[, 1])
    while (remaining_particles > 0) {
      cat(remaining_particles, "\n")
      sample_size <- max(100, remaining_particles)
      candidate_indices <- sample(1:length(previous_par[, 1]),
                                  sample_size,
                                  prob = previous_par[, 7])

      candidate_particles <- previous_par[candidate_indices, ]

      candidate_particles <- t(apply(candidate_particles, 1, mutate_params,
                                     sd_params))

      is_within_prior <- apply(candidate_particles, 1, is_within_prior)

      candidate_particles <- candidate_particles[is_within_prior, ]

      found_trees <- apply(candidate_particles, 1, sim_envidiv_tree,
                           crown_age, TRUE)

      if (length(found_trees) > 0) {

        stat_matrix <- matrix(NA, ncol = 4, nrow = length(found_trees))
        for (i in 1:length(found_trees)) {
          stat_matrix[i, ] <- calc_sum_stats(found_trees[[i]], emp_tree)
        }

        results <- cbind(candidate_particles, stat_matrix)
        results <- results[!is.infinite(results[, 8]), ]

        local_fit <- apply(results[, 8:11], 1, calc_fit, emp_stats)
        results <- cbind(results, local_fit)

        results <- results[results[, 12] < local_eps, ]

        next_par <- rbind(next_par, results)
        remaining_particles <- number_of_particles - length(next_par[, 1])
        if (!is.null(dim(results))) {
          cat(iter, " ", remaining_particles, "\t",
              length(results[, 1]), "\t",
              round(length(results[, 1]) / remaining_particles, 2), "\n")
        }
      }
    }

    next_par <- next_par[1:number_of_particles, ]
    # calculate weights
    new_weights <- apply(next_par[, 1:6], 1, calculate_weight,
                         previous_par[, 1:6], previous_par[, 7], sd_params)
    new_weights <- new_weights / sum(new_weights)
    next_par[, 7] <- new_weights
    previous_par <- next_par
    # write next par to file
    colnames(next_par) <-
        c("extinct", "sym_high", "sym_low", "allo", "jiggle", "model",
            "weight",
            "nltt", "gamma", "mbr", "num_lin",
            "fit")
    next_par <- tibble::as_tibble(next_par)
    readr::write_tsv(next_par, path = paste0("iter_", iter, ".txt"))
  }
}

#' function to perform ABC-SMC
#' @param number_of_particles number of particles used per iteration of
#'                            the SMC algorithm
#' @param max_iter maximum number of iterations
#' @param sd_params standard deviation of the paramater perturbation kernel
#' @param emp_tree phy object holding phylogeny of the tree to be fitted on
#' @param write_to_file (boolean) if true, intermediate output is written to file
#' @return a tibble containing the results
#' @export
infer_params <- function(number_of_particles,
                         max_iter,
                         sd_params,
                         emp_tree,
                         write_to_file = TRUE,
                         seed = NULL) {

  if (is.null(seed)) seed <- as.numeric(Sys.time())
  set.seed(seed)

  param_matrix <- matrix(NA, nrow = number_of_particles,
                         ncol = 7)  #6 parameters

  emp_stats <- calc_sum_stats(emp_tree, emp_tree)
  crown_age <- max(ape::branching.times(emp_tree))

  # generate from prior:
  previous_par <- t(apply(param_matrix, 1, param_from_prior))
  previous_par[, 7] <- previous_par[, 7] / sum(previous_par[, 7])
  next_par <- c()
  # now we start SMC
  for (iter in 2:max_iter) {
    cat("iteration: ", iter, "\n")
    local_eps <- 500 * exp(-0.5 * (iter - 2))

    next_par <- c()

    cat("iteration\tremaining\tfound\taccept rate\tmeanfit\n")

    remaining_particles <- number_of_particles - length(next_par[, 1])
    while (remaining_particles > 0) {
      sample_size <- max(100, remaining_particles)
      candidate_indices <- sample(seq_along(previous_par[, 1]),
                                  sample_size,
                                  prob = previous_par[, 7])

      candidate_particles <- previous_par[candidate_indices, ]

      candidate_particles <- t(apply(candidate_particles, 1, mutate_params,
                                     sd_params))

      is_within_prior <- apply(candidate_particles, 1, is_within_prior)

      candidate_particles <- candidate_particles[is_within_prior, ]


      generate_tree <- function(v) {
        candidate_particle <- v[seq_along(candidate_particles[1, ])]
        rand_seed <- seed + v[1 + length(candidate_particles[1,])]
        return(sim_envidiv_tree(candidate_particle, crown_age, TRUE, rand_seed))
      }

      seeds <- seq_along(candidate_particles[, 1]) + sample(1e15, 1)
      found_trees <- apply(cbind(candidate_particles, seeds), 1, generate_tree)

      if (length(found_trees) > 0) {
        stat_matrix <- matrix(NA, ncol = 8, nrow = length(found_trees))
        for (i in seq_along(found_trees)) {
          stat_matrix[i, ] <- calc_sum_stats(found_trees[[i]], emp_tree)[1:8]
        }


        results <- cbind(candidate_particles, stat_matrix)

        stat_matrix <- stat_matrix[!is.infinite(results[, 8]), ]
        results <- results[!is.infinite(results[, 8]), ]


        local_fit <- apply(stat_matrix, 1, calc_fit, emp_stats)
        results <- cbind(results, local_fit)

        results <- results[local_fit < local_eps, ]
        selected_fits <- local_fit[local_fit < local_eps]

        next_par <- rbind(next_par, results)
        remaining_particles <- number_of_particles - length(next_par[, 1])
        if (!is.null(dim(results))) {
          cat(iter, "\t", remaining_particles, "\t",
              length(results[, 1]), "\t",
              round(length(results[, 1]) / remaining_particles, 2),
              mean(selected_fits), mean(local_fit), "\n")
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
            "beta", "colless", "sackin", "ladder",
            "fit")
    next_par <- tibble::as_tibble(next_par)
    if (write_to_file == TRUE) {
      readr::write_tsv(next_par, path = paste0("iter_", iter, ".txt"))
    }
  }
  return(next_par)
}

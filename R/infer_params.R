#' function to perform ABC-SMC
#' @param number_of_particles number of particles used per iteration of
#'                            the SMC algorithm
#' @param max_iter maximum number of iterations
#' @param sd_params standard deviation of the paramater perturbation kernel
#' @param emp_tree phy object holding phylogeny of the tree to be fitted on
#' @param write_to_file (boolean) if true,
#'                      intermediate output is written to file
#' @param seed seed of the pseudo random-number generator
#' @param continue_from_file (boolean) if true, continues a simulation existing
#'                           in the same folder
#' @param fix_model if -1, the model is fitted, if in [1, 2, 3] it is fitted to
#'                  the specified model (e.g. in 1, 2 or 3)
#' @return a tibble containing the results
#' @export
infer_params <- function(number_of_particles,
                         max_iter,
                         sd_params,
                         emp_tree,
                         write_to_file = TRUE,
                         seed = NULL,
                         continue_from_file = FALSE,
                         fix_model = -1) {

  if (is.null(seed)) seed <- as.numeric(Sys.time())
  set.seed(seed)

  param_matrix <- matrix(NA, nrow = number_of_particles,
                         ncol = 7)  #6 parameters

  emp_stats <- calc_sum_stats(emp_tree, emp_tree)

  # generate from prior:
  previous_par <- t(apply(param_matrix, 1, param_from_prior))

  if (fix_model != -1) previous_par[, 6] <- fix_model

  previous_par[, 7] <- previous_par[, 7] / sum(previous_par[, 7])
  next_par <- c()

  start_iter <- 2

  if (continue_from_file == TRUE) {
    for (i in max_iter:0) {
      file_name <- paste0("iter_", i, ".txt")
      if (file.exists(file_name)) {
        previous_par <- readr::read_tsv(file = file_name)
        previous_par[, 7] <- previous_par[, 7] / sum(previous_par[, 7])
        start_iter <- i
        break
      }
    }

    cat("read previous particles from iteration: ", i, "\n")
  }


  # now we start SMC
  for (iter in start_iter:max_iter) {
    cat("iteration: ", iter, "\n")
    local_eps <- 500 * exp(-0.5 * (iter - 2))

    cat("iteration\tremaining\tfound\taccept rate\tmeanfit\ttime_taken\n")

    next_par <- get_next_par(number_of_particles,
                             previous_par,
                             sd_params,
                             fix_model,
                             emp_tree,
                             emp_stats,
                             local_eps,
                             iter)

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


#' @keywords internal
get_next_par <- function(number_of_particles = 1000,
                         previous_par  = NULL,
                         sd_params = 0.1,
                         fix_model = -1,
                         emp_tree = NULL,
                         emp_stats = NULL,
                         local_eps = Inf,
                         iter = 0)
{
  next_par <- c()
  remaining_particles <- number_of_particles
  accept_rate <- 0
  crown_age <- max(ape::branching.times(emp_tree))

  while (remaining_particles > 0) {

    start_time <- Sys.time()

    sample_size <- max(1000, remaining_particles)
    if (is.null(accept_rate) ||
        is.na(accept_rate)   ||
        length(accept_rate) == 0) {
      accept_rate <- 0
    }

    if (accept_rate != 0) {
      sample_size <- remaining_particles * 1 / accept_rate
    }

    # always need 2, otherwise apply doesn't work
    sample_size <- max(sample_size, 2)
    sample_size <- min(1000, sample_size)

    candidate_indices <- sample(seq_along(previous_par[, 1]),
                                sample_size,
                                prob = previous_par[, 7], replace = T)

    candidate_particles <- previous_par[candidate_indices, ]

    candidate_particles <- t(apply(candidate_particles, 1, mutate_params,
                                   sd_params))

    is_within_prior <- apply(candidate_particles, 1, is_within_prior)

    candidate_particles <- candidate_particles[is_within_prior, ]

    if (is.null(nrow(candidate_particles))) next

    if (fix_model != -1) candidate_particles[, 6] <- fix_model

    input <- lapply(seq_len(nrow(candidate_particles)),
                    function(i) candidate_particles[i, ])

    found_trees <- future.apply::future_lapply(input,
                                               sim_envidiv_tree,
                                               crown_age, TRUE)

    stats <- future.apply::future_lapply(found_trees,
                                         calc_sum_stats,
                                         emp_tree)

    stat_matrix <- matrix(unlist(stats, use.names = FALSE),
                          ncol = 8,
                          byrow = TRUE)

    results <- cbind(candidate_particles, stat_matrix)

    stat_matrix <- stat_matrix[!is.infinite(results[, 8]), ]
    results <- results[!is.infinite(results[, 8]), ]

    if (length(stat_matrix) > 0) {

      local_fit <- c()
      if (!is.null(nrow(stat_matrix))) {
        local_fit <- apply(stat_matrix, 1, calc_fit, emp_stats)
      } else {
        local_fit <- calc_fit(stat_matrix, emp_stats)
      }

      selected_fits <- c()
      if (!is.null(nrow(stat_matrix))) {
        results <- cbind(results, local_fit)

        results <- results[local_fit < local_eps, ]
        selected_fits <- local_fit[local_fit < local_eps]
      } else {
        if (local_fit < local_eps) {
          results <- c(results, local_fit)
          selected_fits <- local_fit
        }
      }

      if (length(results) > 0) {
        rows_before <- nrow(next_par)

        next_par <- rbind(next_par, results)

        num_added_particles <- nrow(next_par) - rows_before
        if (is.null(nrow(next_par)) || is.null(num_added_particles)) {
          num_added_particles <- 0
        }

        remaining_particles <- number_of_particles - nrow(next_par)

        accept_rate <- round(num_added_particles / sample_size, 2)
        end_time <- Sys.time()
        diff_time <- difftime(start_time, end_time, units = "secs")

        cat(iter, "\t",
            remaining_particles, "\t",
            num_added_particles, "\t",
            accept_rate, "\t",
            mean(selected_fits), "\t",
            mean(local_fit), "\t",
            round(diff_time, 1), "\n")
      }
    }
  }

  return(next_par)
}

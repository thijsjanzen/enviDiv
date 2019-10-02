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

  crown_age <- max(ape::branching.times(emp_tree))

  param_matrix <- matrix(NA, nrow = number_of_particles,
                         ncol = 7)  #6 parameters

  emp_stats <- calc_sum_stats(emp_tree, emp_tree)

  # generate from prior:
  previous_par <- t(apply(param_matrix, 1, param_from_prior))

  if(fix_model != -1) previous_par[, 6] <- fix_model

  previous_par[, 7] <- previous_par[, 7] / sum(previous_par[, 7])
  next_par <- c()

  start_iter <- 2

  if(continue_from_file == TRUE) {
    for(i in max_iter:0) {
      file_name <- paste0("iter_", i, ".txt")
      if(file.exists(file_name) ) {
        previous_par <- readr::read_tsv(file = file_name)
        previous_par[, 7] <- previous_par[, 7] / sum(previous_par[, 7])
        start_iter <- i
        break
      }
    }

    #emp_tree <- ape::read.nexus("tree.newick")
    #emp_stats <- calc_sum_stats(emp_tree, emp_tree)

    cat("read previous particles from iteration: ", i, "\n")
  }


  # now we start SMC
  for (iter in start_iter:max_iter) {
    cat("iteration: ", iter, "\n")
    local_eps <- 500 * exp(-0.5 * (iter - 2))

    next_par <- c()

    cat("iteration\tremaining\tfound\taccept rate\tmeanfit\ttime_taken\n")

    remaining_particles <- number_of_particles - length(next_par[, 1])
    accept_rate <- 0
    while (remaining_particles > 0) {

      start_time <- Sys.time()

      sample_size <- max(100, remaining_particles)
      if (accept_rate != 0) {
        sample_size <- 0.25 *  remaining_particles * 1 / accept_rate
      }

      candidate_indices <- sample(seq_along(previous_par[, 1]),
                                  sample_size,
                                  prob = previous_par[, 7], replace = T)

      candidate_particles <- previous_par[candidate_indices, ]

      candidate_particles <- t(apply(candidate_particles, 1, mutate_params,
                                     sd_params))

      is_within_prior <- apply(candidate_particles, 1, is_within_prior)

      candidate_particles <- candidate_particles[is_within_prior, ]

      if(fix_model != -1) candidate_particles[, 6] <- fix_model

      found_trees <- c()


      input <- lapply(1:nrow(candidate_particles), function(i) candidate_particles[i,])

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

      if(length(stat_matrix) > 0) {

        local_fit <- apply(stat_matrix, 1, calc_fit, emp_stats)

        results <- cbind(results, local_fit)

        results <- results[local_fit < local_eps, ]
        selected_fits <- local_fit[local_fit < local_eps]

        next_par <- rbind(next_par, results)
        remaining_particles <- number_of_particles - length(next_par[, 1])
        accept_rate <- round(length(results[, 1]) / remaining_particles, 2)
        if (!is.null(dim(results))) {
          end_time <- Sys.time()
          diff_time <- (end_time - start_time)[[1]]


          cat(iter, "\t", remaining_particles, "\t",
              length(results[, 1]), "\t",
              accept_rate, "\t",
              mean(selected_fits), "\t", mean(local_fit), "\t", diff_time, "\n")
        }
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
  return(next_par)
}

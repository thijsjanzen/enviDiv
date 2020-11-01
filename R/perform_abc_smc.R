#' @keywords internal
get_starting_params <- function(num_particles,
                                min_tips,
                                max_tips,
                                crown_age,
                                num_threads = -1,
                                emp_tree) {

  for_sd <- c()
  output_par <- list()
  for (i in 1:3) {
    sim_result <- enviDiv::create_ref_table_tbb_par(model = i,
                                                  num_repl = num_particles / 3,
                                                    crown_age = crown_age,
                                                    min_lin = min_tips,
                                                    max_lin = max_tips,
                                                    num_threads = num_threads)

    `%dopar%` <- foreach::`%dopar%`


    num_cl <- num_threads
    if (num_threads == -1) num_cl <- parallel::detectCores()

    cl <- parallel::makeForkCluster(num_cl)
    doParallel::registerDoParallel(cl)

    # now we split everything up across threads:
    index_matrix <- split_into_blocks(m = length(sim_result$trees),
                                      block_size = 100)

    index_matrix <- tibble::as_tibble(index_matrix)

    do_analysis <- function(newick_strings, indices_matrix, i) {
      output <- list()
      cnt <- 1
      start <- indices_matrix$lower[[i]]
      end   <- indices_matrix$upper[[i]]
      for (j in start:end) {
        phy <- ape::read.tree(text = newick_strings[j])
        stats <- get_stats_in_order(phy, emp_tree)
        output[[cnt]] <- stats
        cnt <- cnt + 1
      }
      return(output)
    }

    indices <- 1:length(index_matrix$upper)

    results <- foreach::foreach(i = indices)  %dopar% {
      do_analysis(sim_result$trees, index_matrix, i)
    }
    parallel::stopCluster(cl)

    stat_matrix <- matrix(unlist(results, use.names = FALSE),
                          ncol = 10,
                          byrow = TRUE)
    for_sd <- rbind(for_sd, stat_matrix)

    sim_result$parameters[, 6] <-  1 / num_particles
    output_par[[i]] <- sim_result$parameters
  }

  sd_vals <- apply(for_sd, 2, sd)
  output <- list("sd" = sd_vals,
                 "params" = output_par)
  return(output)
}


#' perform abc smc
#' @param emp_tree  phylogenetic tree to fit on
#' @param num_particles number of particles
#' @param num_threads number of threads
#' @param sd_p standard deviation of parameter perturbation
#' @param self_prob probability of drawing self model
#' @param min_tips minimum number of tips
#' @param max_tips maximum number of tips
#' @param write_to_file write to file
#' @return vector of weights
#' @export
perform_abc_smc <- function(emp_tree,
                            num_particles = 1000,
                            num_threads = -1,
                            sd_p = 0.05,
                            sd_self = 0.7,
                            min_tips = 40,
                            max_tips = 120,
                            write_to_file = TRUE) {

  emp_stats <- get_stats_in_order(emp_tree, emp_tree)

  crown_age <- beautier::get_crown_age(emp_tree)

  starting_params <- get_starting_params(num_particles,
                                         min_tips,
                                         max_tips,
                                         crown_age,
                                         num_threads,
                                         emp_tree)

  # we need standard deviations of each parameter
  # and we need to isolate the parameter values

  sd_vals <- starting_params$sd
  params <- starting_params$params

  model_weights <- c()
  max_weights <- c()
  for (i in 1:3) {
    model_weights[i] <- sum(params[[i]][, 6])
    max_weights[i] <- max(params[[i]][, 6])
  }
  model_weights <- model_weights / sum(model_weights)

  thresholds <- 10 * exp(-0.25 * (1:20 - 1))

  for (iter in 1:20) {
    new_params <- list(c(), c(), c())
    number_accepted <- 0
    total_num <- c(0, 0, 0)
    accept_rate <- 1.0
    while (number_accepted < num_particles) {
      num_left <- num_particles - number_accepted
      bsize <- min(1000, (num_left) / accept_rate)

      new_batch <- enviDiv::abc_smc_2(m1 = as.matrix(params[[1]]),
                                      m2 = as.matrix(params[[2]]),
                                      m3 = as.matrix(params[[3]]),
                                      m_weights = model_weights,
                                      max_w = max_weights,
                                      batch_size = bsize,
                                      crown_age = crown_age,
                                      min_lin = min_tips,
                                      max_lin = max_tips,
                                      num_threads = num_threads,
                                      sd_p = sd_p,
                                      self_prob_m = sd_self)

      accept_trees <- enviDiv::accept_from_r(emp_stats,
                                             new_batch$trees,
                                             threshold = thresholds[iter],
                                             sd = sd_vals,
                                             emp_tree,
                                             num_threads = num_threads)

      num_fitting <- sum(accept_trees)
      accept_rate <- num_fitting / bsize
      if (num_fitting > 1) {

        particles <- new_batch$parameters[which(accept_trees == 1), ]
        number_accepted <- number_accepted + length(particles[, 1])
        weights <- enviDiv::calc_weights_R(m1 = (as.matrix(params[[1]])),
                                           m2 = (as.matrix(params[[2]])),
                                           m3 = (as.matrix(params[[3]])),
                                           particles,
                                           model_weights,
                                           sd_p = sd_p,
                                           self_prob = sd_self)


        for (i in 1:length(particles[, 1])) {
          new_particle <- particles[i, ]
          model <- new_particle[6]
          new_particle[6] <- weights[i]
          new_params[[model]] <- rbind(new_params[[model]], new_particle)
          total_num[model] <- total_num[model] + 1
        }
      }
    }

    # update weights
    # first we normalize the weights
    sum_w <- 0
    for(model in 1:3) {
      if (total_num[model] > 0)
        sum_w <- sum_w + sum(new_params[[model]][, 6])
    }
    new_model_weights <- c(0, 0, 0)
    new_max_weights <- c(0, 0, 0)
    for (model in 1:3) {
      if (total_num[model] > 0) {
        new_params[[model]][, 6] <- new_params[[model]][, 6] / sum_w
        new_model_weights[model] <- sum(new_params[[model]][, 6])
        new_max_weights[model]   <- max(new_params[[model]][, 6])
      }
    }

    if (write_to_file) {
      cat(iter, model_weights, "\n", file = "model_weights.txt", append = T)
    }

    params        <- new_params
    model_weights <- new_model_weights
    max_weights   <- new_max_weights
    cat(iter, model_weights, "\n")
    best_model <- which.max(model_weights)

    if (write_to_file) {
      for(x in 1:3) {
        output <- params[[x]]
        output <- cbind(output, x)
        output <- as.matrix(output)
        output <- tibble::as_tibble(output)
        colnames(output) <- c("extinct", "sym_high", "sym_low",
                              "allo", "jiggle", "weight", "model")

        readr::write_tsv(output,
                         path = paste0("particles_", iter, ".txt"),
                         append = TRUE)
      }
    }

  }

  best_model <- which.max(model_weights)
  output <- list("model_weights" = model_weights,
                 "found_params"  = params[[best_model]],
                 "best_model"    = best_model)
  return(output)
}

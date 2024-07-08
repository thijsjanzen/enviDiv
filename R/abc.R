calculate_weight <- function(weights,
                             particles,
                             current,
                             sigma,
                             prior_density_function) {
  vals <- c()
  for (i in seq_len(nrow(particles))) {
    vals[i] <- log(weights[i])

    diff <- log(current[1:6]) - log(particles[i, 1:6])
    prob_diff <- stats::dnorm(diff, mean = 0, sd = sigma, log = TRUE)

    prob_diff[7] <- 0
    if (current[7] != particles[i, 6]) {
      prob_diff[7] <- log(1/3)
    }

    vals[i] <- vals[i] + sum(prob_diff)
  }
  vals <- exp(vals)

  numerator <- prior_density_function(current)
  if (is.na(numerator / sum(vals))) {
    cat(numerator, vals, current,"\n")
  }

  if (is.nan(numerator / sum(vals))) {
    cat(numerator, vals, current,"\n")
  }

  return(numerator / sum(vals))
}

#' abc function
#' @param tree an object of class \code{"phylo"}; the tree upon which we want
#'   to fit our diversification model
#' @param statistics A vector containing functions that take a tree
#'   as an argument and return a single scalar value (the statistic).
#' @param simulation_function A function that implements the
#'   diversification model and returns an object of class \code{"phylo"}.
#' @param init_epsilon_values A vector containing the initial threshold values
#'   for the summary statistics from the vector \code{statistics}.
#' @param prior_generating_function Function to generate parameters
#'   from the prior distribution of these parameters (e.g. a function returning
#'   lambda and mu in case of the birth-death model)
#' @param prior_density_function Function to calculate the prior probability
#'   of a set of parameters.
#' @param number_of_particles Number of particles to be used
#'   per iteration of the ABC-SMC algorithm.
#' @param sigma Standard deviation of the perturbance distribution
#'   (perturbance distribution is a gaussian with mean 0).
#' @param stop_rate If the acceptance rate drops below \code{stopRate},
#'   stop the ABC-SMC algorithm  and assume convergence.
#' @param num_iterations num iterations
#' @return A matrix with \code{n} columns,
#'   where \code{n} is the number of parameters you are trying to estimate.
#' @references  Toni, T., Welch, D., Strelkowa, N., Ipsen, A.,
#'   & Stumpf, M.P.H. (2009). Approximate Bayesian computation scheme for
#'   parameter inference and model selection in dynamical systems.
#'   Journal of the Royal Society Interface, 6(31), 187-202.
#' @export
abc_smc <- function(
  ref_tree,
  statistics,
  simulation_function,
  init_epsilon_value,
  prior_generating_function,
  prior_density_function,
  number_of_particles = 1000,
  sigma = 0.05,
  stop_rate = 1e-5,
  num_iterations = 50
) {
  if (!inherits(ref_tree, "phylo")) {
    # Just checking
    stop("abc_smc_nltt: ",
         "tree must be of class 'phylo', ",
         "but was of type '", class(ref_tree), "' instead")
  }

  #just to get the number of parameters to be estimated.
  parameters <- prior_generating_function()
  num_parameters <- length(parameters)

  # compute the observed statistics
  obs_statistics <- statistics(ref_tree)

  #generate a matrix with epsilon values
  #we assume that the SMC algorithm converges within 50 iterations
  epsilon <- init_epsilon_value * exp(-0.5 * 0:num_iterations)

  #store weights
  new_weights <- c()
  new_params <- list(c(seq_along(parameters)))
  previous_weights <- c()
  previous_params  <- list(c(seq_along(parameters)))
  indices <- 1:number_of_particles

  stats <- c()

  #convergence is expected within 50 iterations
  #usually convergence occurs within 20 iterations

  all_res <- list()
  all_wlevel <- list()

  # first we do the initial generation from the prior
  cat("\nGenerating from the prior\n")
  new_params <- enviDiv::initial_draw_from_prior(num_particles = number_of_particles,
                                        crown_age = treestats::crown_age(ref_tree),
                                        min_lin = 5,
                                        max_lin = 500,
                                        verbose = TRUE)
  new_weights <- rep(1, number_of_particles)


  for (gen in 2:num_iterations) {
    cat("\nGenerating Particles for iteration\t", gen, "\n")
    cat("0--------25--------50--------75--------100\n")
    cat("*")
    utils::flush.console()

    print_frequency <- 20
    tried <- 0
    number_accepted <- 0

    water_levels <- list()

    #replace all vectors
    if (gen > 1) {
      #normalize the weights and store them as previous weights.
      previous_weights <- new_weights / sum(new_weights)
      new_weights <- c() #remove all currently stored weights
      previous_params <- new_params # store found params
      new_params <- matrix(nrow = number_of_particles,
                           ncol = num_parameters) #clear new params
    }

    stoprate_reached <- FALSE

    while (number_accepted < number_of_particles) {

      #in this initial step, generate parameters from the prior
      if (gen == 1) {
        parameters <- prior_generating_function()
      } else {
        #if not in the initial step, generate parameters
        #from the weighted previous distribution:
        index <- sample(x = indices, size = 1,
                        replace = TRUE, prob = previous_weights)

        parameters <- previous_params[index, ]

        #only perturb one parameter, to avoid extremely
        #low acceptance rates due to simultaneous perturbation
        to_change <- sample(seq_along(parameters), 1)

        # perturb the parameter a little bit,
        #on log scale, so parameter doesn't go < 0
        if (to_change != 7) {
          eta <- log(parameters[to_change]) + stats::rnorm(1, 0, sigma)
          parameters[to_change] <- exp(eta)
        } else {
          all_models <- 1:4
          prev_model <- parameters[to_change]
          all_models <- all_models[-prev_model]
          parameters[to_change] <- sample(all_models, size = 1)
        }
      }

      #reject if outside the prior
      if (prior_density_function(parameters) > 0) {
        #simulate a new tree, given the proposed parameters
        new_data <- simulation_function(parameters)
        new_tree <- new_data$phy

        if (inherits(new_tree, "phylo")) {

          stats <- statistics(new_tree)
          #check if the summary statistics are sufficiently
          #close to the observed summary statistics
          #
          accept <- TRUE
          if (gen > 1) {
            diff <- stats - obs_statistics
            rel_diff <- (diff * diff) / abs(obs_statistics)
            misses <- rel_diff > epsilon[gen]

            if (sum(misses, na.rm = TRUE) > 0) accept <- FALSE
          }

          if (accept) {
            number_accepted <- number_accepted + 1
            new_params[number_accepted, ] <- parameters
            accepted_weight <- 1
            water_levels[[number_accepted]] <- new_data$water
            #calculate the weight
            if (gen > 1) {
              accepted_weight <- calculate_weight(previous_weights,
                                                  previous_params, parameters,
                                                  sigma, prior_density_function)
            }
            new_weights[number_accepted] <- accepted_weight

            if ((number_accepted) %%
                (number_of_particles / print_frequency) == 0) {
              cat("**")
              utils::flush.console()
            }
          }
        }
      }

      #convergence if the acceptance rate gets too low
      tried <- tried + 1
      if (tried > (1 / stop_rate)) {
        if ((number_accepted / tried) < stop_rate) {
          stoprate_reached <- TRUE
          break
        }
      }
    }

    all_res[[gen - 1]] <- previous_params
    all_wlevel[[gen]] <- water_levels

    if (stoprate_reached) {
      break
    }
  }

  all_res[[length(all_res) + 1]] <- new_params
  all_wlevel[[length(all_wlevel) + 1]] <- water_levels

  return(list("all_parameters" = all_res,
              "all_water" = all_wlevel))
}



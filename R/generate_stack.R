#' function to perform ABC-SMC
#' @param number_of_replicates number of particles used per iteration of
#'                            the SMC algorithm
#' @param parameters parameters, a vector with:
#' 1) extinction rate,
#' 2) (sympatric) speciation rate at high water,
#' 3) sympatric speciation rate at low water,
#' 4) allopatric speciation rate at low water,
#' 5) posterior perturbation and
#' 6) the specific model, where 1) model without water level changes,
#' 2) literature water level changes, 3) extrapolated water level changes,
#' 4) standard birth-death model.
#' @param min_tips minimum number of tips
#' @param max_tips maximum number of tips
#' @param emp_tree phy object holding phylogeny of the tree to be fitted on
#' @param crown_age crown age
#' @param write_to_file boolean, if TRUE, results are written to file.
#' @param file_name file name
#' @param num_threads number of threads
#' @return a tibble containing the results
#' @export
generate_stack <- function(number_of_replicates = 1000,
                           use_exp_prior = FALSE,
                           focal_model = 1,
                           min_tips = 50,
                           max_tips = 150,
                           emp_tree = NULL,
                           crown_age = NULL,
                           write_to_file = FALSE,
                           file_name = NULL,
                           num_threads = 1) {

  if (!is.null(emp_tree)) {
    crown_age <- max(ape::branching.times(emp_tree))
  }
  if (is.null(crown_age)) {
    stop("Please either provide a reference tree, or provide the crown age")
  }

  number_accepted <- 0
  remaining_particles <- number_of_replicates - number_accepted

  all_results <- c()

  while (remaining_particles > 0) {
    cat(remaining_particles, "\n")
    sample_size <- max(1000, remaining_particles) #increase if not testing

    candidate_particles <- list()

    for (i in 1:sample_size) {
      if (use_exp_prior) {
        candidate_particles[[i]] <- param_from_prior_exp_cpp(focal_model)
      } else {
        candidate_particles[[i]] <- param_from_prior_cpp(focal_model)
      }
    }

    dummy_tree <- ape::rphylo(n = 10, birth = 1, death = 0)
    dummy_stats <- treestats::calc_all_stats(dummy_tree)

    calc_tree_and_stats <- function(x) {
      stats <- rep(NA, length(dummy_stats))

      found_tree <- c()
      if (x[6] == 4) {
        found_tree <- TreeSim::sim.bd.age(age = crown_age,
                                          numbsim = 1,
                                          lambda = x[2],
                                          mu = x[1],
                                          mrca = TRUE,
                                          complete = FALSE)[[1]]
        while (is.numeric(found_tree)) {
          found_tree <- TreeSim::sim.bd.age(age = crown_age,
                                            numbsim = 1,
                                            lambda = x[2],
                                            mu = x[1],
                                            mrca = TRUE,
                                            complete = FALSE)[[1]]
        }
      } else {
        found_tree <- enviDiv::sim_envidiv_tree(x, crown_age, abc = TRUE)
      }

      if (is.null(found_tree)) {
        return(c(x, stats))
      }

      num_tips <- treestats::number_of_lineages(found_tree)

      if (num_tips >= min_tips && num_tips <= max_tips) {
        stats <- treestats::calc_all_stats(found_tree)
      }
      return(c(x, as.vector(stats)))
    }


    res <- list()
    if (num_threads == 1) {

      for (i in 1:length(candidate_particles)) {
        res[[i]] <- calc_tree_and_stats(candidate_particles[[i]])
      }
    } else {
      res <- parallel::mclapply(candidate_particles, calc_tree_and_stats,
                                mc.cores = num_threads,
                                mc.preschedule = FALSE)
    }

    results <- matrix(unlist(res, use.names = FALSE),
                          ncol = length(dummy_stats) + length(candidate_particles[[1]]), # 7 parameters
                          byrow = TRUE)

    results <- results[!is.na(results[, 8]), ]

    if (length(results) > 0) {

      num_local_accepted <- nrow(results)
      if (!is.null(num_local_accepted)) {
        number_accepted <- number_accepted + num_local_accepted

        remaining_particles <- number_of_replicates - number_accepted
        colnames(results) <-
          c("extinct", "sym_high", "sym_low", "allo", "jiggle", "model",
            "weight",
            names(dummy_stats))

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

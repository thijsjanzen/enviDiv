#' @keywords internal
prior_prob <- function(focal_particle) {
  model_prob <- 1 / 3
  particle_prob <- 1
  for (i in 1:5) {
    particle_prob <- particle_prob *
            stats::dunif(x = log10(focal_particle[i]), min = -4, max = 2)
  }
  return(particle_prob * model_prob)
}


#' @keywords internal
calc_w <- function(prev_particles,
                   focal_particle,
                   model_prob,
                   sd_p) {

  total_probs <- c(log10(prev_particles[, 6]))
  for (i in 1:5) {
    diff <- log10(prev_particles[, i]) - log10(focal_particle[i])
    probs <- stats::dnorm(diff, sd = sd_p, log = TRUE)
    total_probs <- cbind(total_probs, probs)
  }

  vv <- exp(rowSums(total_probs))
  total_prob <- sum(vv)
  weight <- 0
  if (total_prob > 1e-6) {
    weight <- prior_prob(focal_particle) / total_prob
  }

  return(weight)
}


#' calculate weights
#' @param m1  previous particles of model 1
#' @param m2  previous particles of model 2
#' @param m3  previous particles of model 3
#' @param p   current particles
#' @param m_w previous model weights
#' @param sd_p standard deviation of parameter perturbation
#' @param self_prob probability of drawing self model
#' @param num_threads number of threads
#' @return vector of weights
#' @export
calc_weights_R <- function(m1,
                           m2,
                           m3,
                           p,
                           m_w,
                           sd_p,
                           self_prob,
                           num_threads = -1) {

  prev_m <- list(m1, m2, m3)
  w <- c()

  model_probs <- rep(0, 3)
  for (i in 1:3) {
    for (j in 1:3) {
      prob_m <- self_prob
      if (i != j) prob_m <- 0.5 * (1 - self_prob)
      model_probs[i] <- model_probs[i] + m_w[j] * prob_m
    }
  }

  w <- rep(0, length(p[, 1]))

  num_cl <- num_threads
  if (num_threads == -1) num_cl <- parallel::detectCores()

  cl <- parallel::makeForkCluster(num_cl)
  doParallel::registerDoParallel(cl)

  # now we split everything up across threads:
  index_matrix <- split_into_blocks(m = length(p[, 1]),
                                    block_size = 100)

  index_matrix <- tibble::as_tibble(index_matrix)

  calc_local_weight <- function(p, indices_matrix, i) {
    output <- list()
    cnt <- 1
    start <- indices_matrix$lower[[i]]
    end   <- indices_matrix$upper[[i]]
    for (j in start:end) {
      focal_p <- as.numeric(p[j, ])
      model <- focal_p[[6]]
      output[[cnt]] <- calc_w(prev_m[[model]],
                              focal_p,
                              model_probs[model],
                              sd_p)
      cnt <- cnt + 1
    }
    return(output)
  }

  `%dopar%` <- foreach::`%dopar%`

  indices <- 1:length(index_matrix$upper)
  w <- foreach::foreach(i = indices)  %dopar% {
    calc_local_weight(p, index_matrix, i)
  }
  parallel::stopCluster(cl)
  w <- unlist(w)
  return(w)
}

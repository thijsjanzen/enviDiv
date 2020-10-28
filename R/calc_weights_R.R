#' @keywords internal
prior_prob <- function(focal_particle) {
  model_prob <- 1/3
  particle_prob <- 1
  for(i in 1:5) {
    particle_prob <- particle_prob *
            dunif(x = log10(focal_particle[1]), min = -4, max = 2)
  }
  return(particle_prob * model_prob)
}


#' @keywords internal
calc_w <- function(prev_particles,
                   focal_particle,
                   model_prob,
                   sd_p) {

  total_probs <- c(log(prev_particles[, 6]))
  for(i in 1:5) {
    diff <- log10(prev_particles[, i]) - log10(focal_particle[i])
    probs <- dnorm(diff, sd = sd_p, log = TRUE)
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
#' @return vector of weights
#' @export
calc_weights_R <- function(m1, m2, m3, p, m_w, sd_p, self_prob) {

  prev_m <- list(m1, m2, m3)
  w <- c()

  model_probs <- rep(0, 3)
  for(i in 1:3) {
    for(j in 1:3) {
      prob_m <- self_prob
      if (i != j) prob_m <- 0.5 * (1 - self_prob)
      model_probs[i] = model_probs[i] + m_w[j] * prob_m
    }
  }

  for(i in 1:length(p[, 1])) {
    focal_p <- p[i, ]
    model <- focal_p[6]

    # prev_particles <- prev_m[[model]]
    # focal_particle <- focal_p
    # model_prob <- model_probs[model]
    # sd_p <- sd_p

    w[i] = calc_w(prev_m[[model]],
                  focal_p,
                  model_probs[model],
                  sd_p)
  }
  return(w)
}

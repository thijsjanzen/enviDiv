#' Function to calculate the weight of a particle.
#'   The weight depends on the probability of obtaining the specific parameter
#'   combination, given all the previous parameters. Hence, it is the sum of the
#'   probability for each previous parameter combination to obtain the observed
#'   parameter combination.
#' @param params focal parameters for which to estimate the weight
#' @param other_params previous params from which this parameter
#'                     combination is derived
#' @param weights previous weights
#' @param sigma standard deviation of the perturbation kernel
#' @return weight
#' @export
calculate_weight <- function(params,
                             other_params,
                             weights,
                             sigma) {
  #only the numerical ones
  focal_rates <- log10(params[1:5])
  other_rates <- log10(other_params[, 1:5])
  diff <- t(apply(other_rates, 1, function(x) x - focal_rates))

  if (sum(weights) != 1) {
    weights <- weights / sum(weights)
  }

  vals <- weights * stats::dnorm(diff, mean = 0, sd = sigma, log = FALSE)

  # include weight of model
  model_prob <- params[6] - other_params[, 6]
  model_prob[which(model_prob != 0)] <- 0.25
  model_prob[which(model_prob == 0)] <- 0.5

  final_weight <- 1.0 / (sum(model_prob) * sum(vals) )

  # do something with model weights
  return(final_weight)
}

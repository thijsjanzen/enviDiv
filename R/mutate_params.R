#' perturbate all parameters
#' @param params input vector containing unperturbed parameters
#' @param local_sd standard deviation of the perturbation kernel
#' @return perturbed parameters
#' @export
mutate_params <- function(params, local_sd) {
  output <- c()
  for (i in 1:5) {
    output[i] <- 10 ^ (stats::rnorm(1, log10(params[i]), local_sd))
  }

  probs <- rep(0.25, 3)
  probs[params[6]] <- 0.5
  output[6] <- sample(1:3, 1, prob = probs)

  output[7] <- params[7]
  return(output)
}

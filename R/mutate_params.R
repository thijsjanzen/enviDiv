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

  output[6] <- params[6]
  r <- sample(1:4, 1)
  if (r == 0) output[6] <- params[6] + 1
  if (r == 1) output[6] <- params[6] - 1
  if (output[6] > 2) output[6] <- output[6] - 3
  if (output[6] < 1) output[6] <- output[6] + 3

  output[7] <- params[7]
  return(output)
}

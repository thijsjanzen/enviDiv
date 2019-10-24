#' draw parameter combinations from the prior
#' @param empty_input empty vector, used to facilitate usage of apply
#' @return vector with 7 entries: extinction, sympatric speciation at high
#' water, sympatric speciation at low water, allopatric speciation, amount of
#' perturbation, the chosen water model and the abc-weight
#' @export
param_from_prior <- function(empty_input) {
  output <- c()
  output[1] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #extinction
  output[2] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #symp spec high
  output[3] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #symp spec low
  output[4] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #allo spec
  output[5] <- 10 ^ (-3 + 3 * stats::runif(1, 0, 1))  #jiggle
  output[6] <- sample(1:3, 1)
  output[7] <- 1

  return(output)
}

#' draw parameter combinations from the prior
#' @param empty_input empty vector, used to facilitate usage of apply
#' @return vector with 7 entries: extinction, sympatric speciation at high
#' water, sympatric speciation at low water, allopatric speciation, amount of
#' perturbation, the chosen water model and the abc-weight
#' @export
param_from_prior_exp <- function(empty_input) {
  output <- c()
  output[1] <- stats::rexp(n = 1, rate = 1 / 0.05)  # extinction
  output[2] <- stats::rexp(n = 1, rate = 1 / 0.10)  # symp spec high
  output[3] <- stats::rexp(n = 1, rate = 1 / 0.10)  # symp spec low
  output[4] <- stats::rexp(n = 1, rate = 1 / 0.10)  # allo spec
  output[5] <- stats::rexp(n = 1, rate = 1 / 0.05)  # jiggle
  output[6] <- sample(1:3, 1)
  output[7] <- 1

  return(output)
}

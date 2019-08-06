#' draw parameter combinations from the prior
#' @param model requested water model, default is to sample from 1:3
#' @return vector with 7 entries: extinction, sympatric speciation at high water, sympatric speciation at low water, allopatric speciation, amount of perturbation, the chosen water model and the abc-weight
#' @export
param_from_prior <- function(model = NULL) {
  output <- rep(NA, 7)
  output[1] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #extinction
  output[2] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #symp spec high
  output[3] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #symp spec low
  output[4] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #allo spec
  output[5] <- 10 ^ (-3 + 3 * stats::runif(1, 0, 1))  #jiggle
  if(is.null(model)) output[6] <- sample(1:3, 1) # model
  output[7] <- 1

  return(output)
}

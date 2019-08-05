#' draw parameter combinations from the prior
#' @param model requested water model
#' @return vector with 7 entries: extinction, sympatric speciation at high water, sympatric speciation at low water, allopatric speciation, amount of jiggle and the chosen water model
#' @export
param_from_prior <- function(model = NULL) {
  output <- c()
  output[1] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #extinction
  output[2] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #symp spec high
  output[3] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #symp spec low
  output[4] <- 10 ^ (-3 + 5 * stats::runif(1, 0, 1))  #allo spec
  output[5] <- 10 ^ (-3 + 3 * stats::runif(1, 0, 1))  #jiggle
  output[6] <- model
  if(is.null(model)) output[6] <- sample(1:3, 1) # model
  output[7] <- 1

  return(output)
}

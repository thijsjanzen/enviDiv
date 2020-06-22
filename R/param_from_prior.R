#' draw parameter combinations from the prior
#' @param empty_input empty vector, used to facilitate usage of apply
#' @return vector with 7 entries: extinction, sympatric speciation at high
#' water, sympatric speciation at low water, allopatric speciation, amount of
#' perturbation, the chosen water model and the abc-weight
#' @export
param_from_prior <- function(empty_input) {
  output <- param_from_prior_cpp();

  return(output)
}

#' draw parameter combinations from the prior
#' @param empty_input empty vector, used to facilitate usage of apply
#' @return vector with 7 entries: extinction, sympatric speciation at high
#' water, sympatric speciation at low water, allopatric speciation, amount of
#' perturbation, the chosen water model and the abc-weight
#' @export
param_from_prior_exp <- function(empty_input) {
  output <- param_from_prior_exp_cpp();

  return(output)
}

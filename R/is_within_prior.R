#' function to determine if the parameter combination is within the prior
#' @param params vector of parameters
#' @return true or false
#' @export
is_within_prior <- function(params) {
  for(i in 1:4) {
    if(log10(params[i]) < -3) return(FALSE)
    if(log10(params[i]) >  2) return(FALSE)
  }
  if(log10(params[5]) < -3) return(FALSE)
  if(log10(params[5]) >  0) return(FALSE)

  if(params[6] < 1) return(FALSE)
  if(params[6] > 3) return(FALSE)

  return(TRUE)
}

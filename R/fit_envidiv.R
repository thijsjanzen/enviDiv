#' Maximum-likelihood fit of the envidiv model for a known water history.
#' @description
#' Parameters are optimised on log scale.
#' @param phy a phylo object
#' @param init initial values
#' @param fixed values not optimised, e.g. c(mu = 0) for a model without extinction
#' @param tol tolerance
#' @param maxit maximum iterations of the ML algorithm
#' @param trace verbose output of the optimisation
#' @param return maximum likelihood
#' @export
fit_envidiv <- function(phy, stands,
                        init = c(mu = 0.1, lambda_H = 0.5,
                                 lambda_L = 0.5, alpha = 0.5),
                        fixed = c(),          # e.g. c(alpha = 0) for a no-pump model
                        tol = 1e-6, maxit = 2000, trace = FALSE) {
  free <- setdiff(names(init), names(fixed))
  negll <- function(lp) {
    pars <- c(exp(lp), fixed)[c("mu", "lambda_H", "lambda_L", "alpha")]
    ll <- tryCatch(envidiv_loglik(pars, phy, stands, tol = tol),
                   error = function(e) -Inf)
    if (!is.finite(ll)) return(1e10)
    if (trace) cat(signif(pars, 4), ll, "\n")
    -ll
  }
  opt <- if (requireNamespace("subplex", quietly = TRUE)) {
    subplex::subplex(log(init[free]), negll, control = list(maxit = maxit))
  } else {
    stats::optim(log(init[free]), negll, control = list(maxit = maxit))
  }
  est <- c(exp(opt$par), fixed)[c("mu", "lambda_H", "lambda_L", "alpha")]
  names(est) <- c("mu", "lambda_H", "lambda_L", "alpha")
  list(par = est, loglik = -opt$value, df = length(free),
       AIC = 2 * opt$value + 2 * length(free), convergence = opt$convergence)
}

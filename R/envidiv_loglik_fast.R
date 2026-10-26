# R wrapper around the C++ likelihood in envidiv_loglik.cpp.
# sigma = 0: exact likelihood (nodes within tol of a drop age are pump nodes);
# sigma > 0: mixture likelihood with dating error on pump nodes.

#' Convert a simulator water vector (forward times of level changes, as
#' returned by enviDiv::generate_water) into lowstands in ages.
#' Returns a data.frame with columns rise (younger age) and drop (older age).
water_to_stands <- function(water_changes, crown_age) {
  toggles <- sort(water_changes[water_changes > 0 & water_changes < crown_age])
  if (length(toggles) == 0) {
    return(data.frame(rise = numeric(0), drop = numeric(0)))
  }
  drops <- toggles[seq(1, length(toggles), by = 2)]
  rises <- toggles[seq(2, length(toggles), by = 2)]
  if (length(rises) < length(drops)) rises <- c(rises, crown_age)
  st <- data.frame(rise = crown_age - rises, drop = crown_age - drops)
  st[order(st$rise), , drop = FALSE]
}

#' Log-likelihood (C++; closed forms plus Boost odeint), version 4.
#' @description
#' very fast implementation of the LL, using Rcpp.
#' Written by Claude 5.5
#'
#' @param pars        named vector mu, lambda_H, lambda_L, alpha (rates per unit of *true* time)
#' @param stands      data.frame(rise, drop) in true (geological) ages
#' @param sigma       dating error of pump nodes (true time units); 0 = exact pump nodes
#' @param f           sampling fraction of extant species
#' @param time_scale  c in observed age = c * true age (global calibration factor);
#'                    node ages are divided by c, and the Jacobian c^-(n-2) of the
#'                    non-crown node ages is included
#' @param method      odeint stepper for the parts without closed form:
#'                    "odeint::runge_kutta_cash_karp54" (default),
#'                    "odeint::runge_kutta_fehlberg78", "odeint::runge_kutta_dopri5",
#'                    "odeint::bulirsch_stoer"  r"
#' @param atol        absolute tolerance of the integration
#' @param rtol        relative tolerance of the integration
#' @export
#' @return loglikelihood
envidiv_loglik_fast <- function(pars,
                                phy,
                                stands,
                                sigma = 0,
                                f = 1,
                                time_scale = 1,
                                max_z = 5,
                                tol = 1e-6,
                                method = "odeint::runge_kutta_cash_karp54",
                                atol = 1e-12,
                                rtol = 1e-10) {
  n_tip <- ape::Ntip(phy)
  ages <- ape::node.depth.edgelength(phy)
  ages <- max(ages) - ages
  ages[seq_len(n_tip)] <- 0
  ages <- ages / time_scale
  envidiv_loglik_cpp(pars[c("mu", "lambda_H", "lambda_L", "alpha")], phy$edge, ages,
                     n_tip, stands$rise, stands$drop, sigma, max_z, tol, f,
                     method, atol, rtol) -
    (n_tip - 2) * log(time_scale)
}

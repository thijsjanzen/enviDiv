# Reversible-jump MCMC for the envidiv model with an unknown water-level
# history. Requires envidiv_loglik_fast.R (C++ likelihood with dating error).
#
# State: rates (mu, lambda_H, lambda_L, alpha), sigma, water-change rate r,
# and a set of J non-overlapping lowstands (rise_j < drop_j, ages in (0, T)).
#
# Prior on the water history: level changes form a Poisson process with rate
# r (as in enviDiv water models 3/4), the lake is high at the crown and at
# the present (even number of changes). For J lowstands the density of the
# sorted change times is r^(2J) exp(-rT) / P(even | r).
# Priors on rates, sigma and r are log-uniform within bounds.
#
# Moves: log random walks on rates, sigma and r; shifts of a single drop or
# rise time; shifts of a whole lowstand; birth and death of a lowstand with a
# data-driven proposal for the drop time (near observed node ages).

rj_defaults <- function() list(
  rate_bounds = c(1e-3, 1e2), sigma_bounds = c(1e-3, 0.5), r_bounds = c(0.1, 50),
  move_probs = c(rates = 0.30, sigma = 0.05, r = 0.05, shift = 0.20,
                 shift_stand = 0.10, birth = 0.15, death = 0.15),
  rate_step = 0.3, sigma_step = 0.3, r_step = 0.5,
  shift_scales = c(0.002, 0.02, 0.2),
  birth_sd = 0.01, birth_unif = 0.1, dur_min = 1e-3)

log_prior_water <- function(J, r, T) {
  2 * J * log(r) - r * T - log((1 + exp(-2 * r * T)) / 2)
}

valid_stands <- function(st, T) {
  if (nrow(st) == 0) return(TRUE)
  if (any(st$rise <= 0) || any(st$drop >= T) || any(st$rise >= st$drop)) return(FALSE)
  if (nrow(st) > 1 && any(st$rise[-1] <= st$drop[-nrow(st)])) return(FALSE)
  TRUE
}

#' Reversible-jump MCMC for the envidiv model with an unknown water-level
#' @description
#' Reversible-jump MCMC for the envidiv model with an unknown water-level
#' @param phy phy
#' @param n_iter number of iterations
#' @param thin thinning interval
#' @param init initial values for pars, sigma, r and stands
#' @param control list of control parameters (see rj_defaults)
#' @param prior_only if TRUE, only the prior is evaluated (no likelihood)
#' @param fix character vector of parameters to fix (e.g. c("mu = 0"))
#' @param verbose if TRUE, print progress every 200 samples
#' @return list with trace, stands, acceptance rates, T and control
#' @export
rjmcmc_envidiv <- function(phy, n_iter = 1e5, thin = 50,
                           init = list(pars = c(mu = 0.1, lambda_H = 0.5,
                                                lambda_L = 0.5, alpha = 0.5),
                                       sigma = 0.02, r = 1,
                                       stands = data.frame(rise = numeric(0),
                                                           drop = numeric(0))),
                           control = list(), prior_only = FALSE,
                           fix = character(0), verbose = TRUE) {
  ctl <- modifyList(rj_defaults(), control)
  ages <- ape::node.depth.edgelength(phy); ages <- max(ages) - ages
  T <- max(ages)
  n_tip <- ape::Ntip(phy)
  node_ages <- ages[(n_tip + 2):(n_tip + phy$Nnode)]   # internal, excluding crown
  dur_max <- T
  log_q_drop <- function(a) {
    log((1 - ctl$birth_unif) * mean(dnorm(a, node_ages, ctl$birth_sd)) + ctl$birth_unif / T)
  }
  log_q_dur <- function(d) {
    if (d < ctl$dur_min || d > dur_max) return(-Inf)
    -log(d) - log(log(dur_max / ctl$dur_min))
  }
  loglik <- function(pars, sigma, st) {
    if (prior_only) return(0)
    envidiv_loglik_fast(pars, phy, st, sigma)
  }
  in_bounds <- function(x, b) all(x >= b[1] & x <= b[2])

  pars <- init$pars; sigma <- init$sigma; r <- init$r
  st <- init$stands[order(init$stands$rise), , drop = FALSE]; rownames(st) <- NULL
  ll <- loglik(pars, sigma, st)
  mp <- ctl$move_probs
  if ("r" %in% fix) mp["r"] <- 0
  if ("rates" %in% fix) mp["rates"] <- 0
  if ("sigma" %in% fix) mp["sigma"] <- 0
  mp <- mp / sum(mp)
  p_b <- mp[["birth"]]; p_d <- mp[["death"]]
  moves <- names(mp)

  n_out <- floor(n_iter / thin)
  trace <- data.frame(iter = integer(n_out), loglik = NA_real_, J = NA_integer_,
                      mu = NA_real_, lambda_H = NA_real_, lambda_L = NA_real_,
                      alpha = NA_real_, sigma = NA_real_, r = NA_real_)
  stand_samples <- vector("list", n_out)
  acc <- setNames(numeric(length(moves)), moves); tried <- acc
  t0 <- Sys.time()

  for (it in seq_len(n_iter)) {
    mv <- sample(moves, 1, prob = mp)
    tried[mv] <- tried[mv] + 1
    J <- nrow(st)
    accept <- FALSE
    if (mv == "rates") {
      k <- sample(4, 1); new <- pars
      new[k] <- pars[k] * exp(rnorm(1, 0, ctl$rate_step))
      if (in_bounds(new[k], ctl$rate_bounds)) {
        ll_new <- loglik(new, sigma, st)
        if (log(runif(1)) < ll_new - ll) { pars <- new; ll <- ll_new; accept <- TRUE }
      }
    } else if (mv == "sigma") {
      new <- sigma * exp(rnorm(1, 0, ctl$sigma_step))
      if (in_bounds(new, ctl$sigma_bounds)) {
        ll_new <- loglik(pars, new, st)
        if (log(runif(1)) < ll_new - ll) { sigma <- new; ll <- ll_new; accept <- TRUE }
      }
    } else if (mv == "r") {
      new <- r * exp(rnorm(1, 0, ctl$r_step))
      if (in_bounds(new, ctl$r_bounds) &&
          log(runif(1)) < log_prior_water(J, new, T) - log_prior_water(J, r, T)) {
        r <- new; accept <- TRUE
      }
    } else if (mv == "shift" && J > 0) {
      j <- sample.int(J, 1); col <- sample(c("rise", "drop"), 1)
      new_st <- st
      new_st[j, col] <- st[j, col] + rnorm(1, 0, sample(ctl$shift_scales, 1))
      if (valid_stands(new_st, T)) {
        ll_new <- loglik(pars, sigma, new_st)
        if (log(runif(1)) < ll_new - ll) { st <- new_st; ll <- ll_new; accept <- TRUE }
      }
    } else if (mv == "shift_stand" && J > 0) {
      j <- sample.int(J, 1); new_st <- st
      new_st[j, ] <- st[j, ] + rnorm(1, 0, sample(ctl$shift_scales, 1))
      if (valid_stands(new_st, T)) {
        ll_new <- loglik(pars, sigma, new_st)
        if (log(runif(1)) < ll_new - ll) { st <- new_st; ll <- ll_new; accept <- TRUE }
      }
    } else if (mv == "birth") {
      a <- if (runif(1) < ctl$birth_unif) runif(1, 0, T) else
        sample(node_ages, 1) + rnorm(1, 0, ctl$birth_sd)
      d <- exp(runif(1, log(ctl$dur_min), log(dur_max)))
      new_st <- rbind(st, data.frame(rise = a - d, drop = a))
      new_st <- new_st[order(new_st$rise), , drop = FALSE]; rownames(new_st) <- NULL
      if (valid_stands(new_st, T)) {
        ll_new <- loglik(pars, sigma, new_st)
        log_a <- ll_new - ll + 2 * log(r) +
          log(p_d / (J + 1)) - (log(p_b) + log_q_drop(a) + log_q_dur(d))
        if (log(runif(1)) < log_a) { st <- new_st; ll <- ll_new; accept <- TRUE }
      }
    } else if (mv == "death" && J > 0) {
      j <- sample.int(J, 1)
      a <- st$drop[j]; d <- st$drop[j] - st$rise[j]
      new_st <- st[-j, , drop = FALSE]; rownames(new_st) <- NULL
      ll_new <- loglik(pars, sigma, new_st)
      log_a <- ll_new - ll - 2 * log(r) +
        (log(p_b) + log_q_drop(a) + log_q_dur(d)) - log(p_d / J)
      if (log(runif(1)) < log_a) { st <- new_st; ll <- ll_new; accept <- TRUE }
    }
    if (accept) acc[mv] <- acc[mv] + 1

    if (it %% thin == 0) {
      i <- it / thin
      trace[i, ] <- c(it, ll, nrow(st), pars, sigma, r)
      stand_samples[[i]] <- st
      if (verbose && i %% 200 == 0) {
        cat(sprintf("iter %d  ll %.1f  J %d  (%.0f s)\n", it, ll, nrow(st),
                    as.numeric(Sys.time() - t0, units = "secs")))
      }
    }
  }
  list(trace = trace, stands = stand_samples, acceptance = acc / pmax(tried, 1),
       T = T, control = ctl)
}

#' Posterior probability that a lowstand starts (drop) or is ongoing at each
#' age on a grid; `window` is the half-width for drop-time intensity.
summarise_water <- function(fit, burnin = 0.2, grid = NULL, window = 0.005) {
  keep <- seq(ceiling(burnin * length(fit$stands)) + 1, length(fit$stands))
  if (is.null(grid)) grid <- seq(0, fit$T, length.out = 1001)
  drop_prob <- sapply(grid, function(g) mean(sapply(fit$stands[keep], function(s)
    any(abs(s$drop - g) <= window))))
  low_prob <- sapply(grid, function(g) mean(sapply(fit$stands[keep], function(s)
    any(s$rise < g & g < s$drop))))
  data.frame(age = grid, drop_prob = drop_prob, low_prob = low_prob)
}

# Exact likelihood of a reconstructed tree under the envidiv model, for a
# KNOWN water-level history (water models 1 and 2, or 3/4 conditional on W).
#
# Parameters (same meaning as in the simulator):
#   mu        per-population extinction rate
#   lambda_H  per-population sympatric speciation rate, high water
#   lambda_L  per-population sympatric speciation rate, low water
#   alpha     per-population allopatric rate (a pair resolves at 2 * alpha)
#
# Ages s run backwards from the present (s = 0) to the crown (s = T).
# Lineage types: H (high water), S (one basin), U (unresolved pair).
# All E and D_S / D_H propagation is analytic; only E_U, D_U and the pump
# integrals use one-dimensional quadrature.
#
# Branch states:
#   high water:  logD (= log D_H)
#   low water:   logD (= log D_S) and rho = D_U / D_S

# ---------------------------------------------------------------- helpers --

g_fun <- function(r, dt) {
  # (exp(r dt) - 1) / r, stable for r -> 0
  ifelse(abs(r * dt) < 1e-12, dt, expm1(r * dt) / r)
}

# survival probability F = 1 - E after time dt, constant rates
prop_F <- function(F0, lambda, mu, dt) {
  r <- lambda - mu
  F0 * exp(r * dt) / (1 + lambda * F0 * g_fun(r, dt))
}

# change in log D on a branch with no events, constant rates
prop_logD <- function(F0, lambda, mu, dt) {
  r <- lambda - mu
  r * dt - 2 * log1p(lambda * F0 * g_fun(r, dt))
}

# ------------------------------------------------------- water regimes ----

#' Convert a simulator water vector (forward times of level changes, as
#' returned by enviDiv::generate_water) into lowstands in ages.
#' Returns a data.frame with columns rise (younger age) and drop (older age).
#' @param water_changes water changes in forward time
#' @param crown_age crown age
#' @export
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

# ---------------------------------------------------- global E functions --

#' Precompute the extinction probabilities needed everywhere.
make_env <- function(pars, stands) {
  mu <- pars[["mu"]]; lH <- pars[["lambda_H"]]
  lL <- pars[["lambda_L"]]; al <- pars[["alpha"]]
  theta <- al + lL + mu
  ns <- nrow(stands)
  FH_rise <- numeric(ns)   # 1 - E_H at the rise age of each stand
  FU_drop <- numeric(ns)   # 1 - E_U at the drop age of each stand
  prevF <- 1; prevS <- 0   # start at the present
  for (j in seq_len(ns)) {
    a1 <- stands$rise[j]; a2 <- stands$drop[j]
    FH_rise[j] <- prop_F(prevF, lH, mu, a1 - prevS)
    F0 <- FH_rise[j]
    # 1 - E_U obeys F_U' = -2 theta F_U + 2 theta - g  with g as in the notes
    integrand <- function(u) {
      ES <- 1 - prop_F(F0, lL, mu, u - a1)
      gg <- 2 * mu * ES + 2 * al * ES^2 + 2 * lL * ES^3
      exp(-2 * theta * (a2 - u)) * gg
    }
    EU <- exp(-2 * theta * (a2 - a1)) * (1 - F0) +
      stats::integrate(integrand, a1, a2, rel.tol = 1e-10)$value
    FU_drop[j] <- 1 - EU
    prevF <- FU_drop[j]; prevS <- a2
  }
  list(mu = mu, lH = lH, lL = lL, al = al, theta = theta,
       stands = stands, FH_rise = FH_rise, FU_drop = FU_drop)
}

stand_of <- function(env, s, tol = 1e-9) {
  # index of the lowstand strictly containing age s (a branch starting at a
  # drop age, i.e. above a pump node, is in high water)
  j <- which(s > env$stands$rise + tol & s < env$stands$drop - tol)
  if (length(j) == 0) return(0L)
  j[1]
}

FH_at <- function(env, s) {
  st <- env$stands
  below <- which(st$drop <= s)
  if (length(below) == 0) {
    F0 <- 1; s0 <- 0
  } else {
    j <- max(below); F0 <- env$FU_drop[j]; s0 <- st$drop[j]
  }
  prop_F(F0, env$lH, env$mu, s - s0)
}

FS_at <- function(env, j, s) {
  prop_F(env$FH_rise[j], env$lL, env$mu, s - env$stands$rise[j])
}

# --------------------------------------------------- branch propagation ---

# propagate a low-water state (logD_S, rho) inside stand j from s0 to s1
prop_low <- function(env, j, s0, s1, logD, rho) {
  if (s1 <= s0) return(list(logD = logD, rho = rho))
  F0 <- FS_at(env, j, s0)
  lL <- env$lL; mu <- env$mu; al <- env$al; th <- env$theta
  dl1 <- prop_logD(F0, lL, mu, s1 - s0)        # log D_S(s1)/D_S(s0)
  integrand <- function(u) {
    ES <- 1 - prop_F(F0, lL, mu, u - s0)
    h <- 2 * mu + 4 * al * ES + 6 * lL * ES^2
    dl <- prop_logD(F0, lL, mu, u - s0)        # log D_S(u)/D_S(s0)
    exp(-2 * th * (s1 - u) + dl - dl1) * h
  }
  rho1 <- exp(-2 * th * (s1 - s0) - dl1) * rho +
    stats::integrate(integrand, s0, s1, rel.tol = 1e-10)$value
  list(logD = logD + dl1, rho = rho1)
}

#' Propagate a branch from age s_lo (with given starting state) up to s_hi.
#' state: list(logD, rho) where rho = NA in high water.
#' If stop_at_rise_of is a stand index j, propagation stops at
#' max(s_lo, rise_j) and returns the low-water state there (used below a
#' pump node).
prop_branch <- function(env, s_lo, s_hi, state, stop_at_rise_of = 0L) {
  st <- env$stands
  s <- s_lo
  j <- stand_of(env, s)
  # a branch starting exactly at a rise age is in high water just below it
  while (TRUE) {
    if (j == 0) {
      # high water: next event is entering a stand (rise age) or s_hi
      nxt <- st$rise[st$rise >= s - 1e-12 & st$rise < s_hi]
      nxt <- if (length(nxt)) min(nxt) else s_hi
      state$logD <- state$logD +
        prop_logD(FH_at(env, s), env$lH, env$mu, nxt - s)
      s <- nxt
      if (s >= s_hi - 1e-12) return(state)
      j <- which(abs(st$rise - s) < 1e-12)[1]
      state$rho <- 1            # D_S = D_U = D_H at the rise
      if (stop_at_rise_of == j) return(c(state, s = s))
    } else {
      top <- min(st$drop[j], s_hi)
      if (stop_at_rise_of == j) return(c(state, s = s))
      state <- prop_low(env, j, s, top, state$logD, state$rho)
      s <- top
      if (s >= s_hi - 1e-12) return(state)
      # leaving the stand at its drop age (going older): U -> H
      state <- list(logD = state$logD + log(state$rho), rho = NA)
      j <- 0L
    }
  }
}

# ------------------------------------------------------------ the tree ----

#' Log-likelihood of an ultrametric phylo under envidiv with known stands.
#' Conditioned on survival of both crown lineages (as in the simulator).
#' Nodes within `tol` of a drop age are treated as pump nodes.
#' @param pars pars
#' @param phy phy
#' @param stands stands
#' @param tol tolerance
#' @export
envidiv_loglik <- function(pars, phy, stands, tol = 1e-6) {
  env <- make_env(pars, stands)
  n_tip <- ape::Ntip(phy)
  ages <- ape::node.depth.edgelength(phy)
  ages <- max(ages) - ages
  ages[seq_len(n_tip)] <- 0
  crown <- n_tip + 1
  if (stand_of(env, ages[crown]) > 0) stop("crown inside a lowstand")
  is_pump <- function(v) {
    if (v <= n_tip) return(0L)
    j <- which(abs(env$stands$drop - ages[v]) < tol)
    if (length(j)) j[1] else 0L
  }
  # snap pump nodes exactly onto drop ages
  for (v in (n_tip + 1):(n_tip + phy$Nnode)) {
    j <- is_pump(v); if (j > 0) ages[v] <- env$stands$drop[j]
  }
  children <- split(phy$edge[, 2], phy$edge[, 1])

  # returns the state of the branch *at its lower node* (age of node v)
  node_state <- function(v) {
    if (v <= n_tip) return(list(logD = 0, rho = NA))
    kids <- children[[as.character(v)]]
    x <- ages[v]
    jp <- is_pump(v)
    if (jp > 0) {
      # pump node: both daughters pending until resolution at tau
      a1 <- env$stands$rise[jp]; a2 <- env$stands$drop[jp]
      st <- lapply(kids, function(k) {
        prop_branch(env, ages[k], a2, node_state(k), stop_at_rise_of = jp)
      })
      tmin <- max(st[[1]]$s, st[[2]]$s)
      st <- lapply(st, function(b) {
        if (b$s < tmin) prop_low(env, jp, b$s, tmin, b$logD, b$rho) else b
      })
      F0 <- FS_at(env, jp, tmin)
      lL <- env$lL; mu <- env$mu; al <- env$al; th <- env$theta
      integrand <- function(tau) {
        ES <- 1 - prop_F(F0, lL, mu, tau - tmin)
        dl <- prop_logD(F0, lL, mu, tau - tmin)
        exp(-2 * th * (a2 - tau) + 2 * dl) * (2 * al + 4 * lL * ES)
      }
      val <- stats::integrate(integrand, tmin, a2, rel.tol = 1e-10)$value
      if (tmin > a1 + 1e-12) {
        # older first daughter node lies inside the stand: point mass
        val <- val + exp(-2 * th * (a2 - tmin))
      }
      return(list(logD = st[[1]]$logD + st[[2]]$logD + log(val), rho = NA))
    }
    # regular node
    st <- lapply(kids, function(k) prop_branch(env, ages[k], x, node_state(k)))
    j <- stand_of(env, x)
    if (j == 0) {
      return(list(logD = log(env$lH) + st[[1]]$logD + st[[2]]$logD, rho = NA))
    }
    ES <- 1 - FS_at(env, j, x)
    list(logD = log(env$lL) + st[[1]]$logD + st[[2]]$logD, rho = 2 * ES)
  }

  kids <- children[[as.character(crown)]]
  T <- ages[crown]
  ll <- 0
  for (k in kids) ll <- ll + prop_branch(env, ages[k], T, node_state(k))$logD
  ll - 2 * log(FH_at(env, T))
}

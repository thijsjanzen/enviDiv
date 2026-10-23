# Can the tree locate a lowstand? Remove one true lowstand from the model and
# scan the log-likelihood over candidate drop (and rise) times.
suppressMessages({library(enviDiv); library(ape)})
source("envidiv_loglik_fast.R")
T <- 5; pv <- c(0.2, 0.6, 1.0, 1.5); sig <- 0.02
pars <- c(mu = pv[1], lambda_H = pv[2], lambda_L = pv[3], alpha = pv[4])
W <- generate_water(2, T); st <- water_to_stands(W, T); rownames(st) <- NULL
ntree <- as.integer(commandArgs(TRUE)[1]); set.seed(1)
grid <- seq(0.05, 1.6, by = 0.0025)
res <- list()
for (i in seq_len(ntree)) {
  repeat { r <- create_tree_cpp(c(pv, sig, 2), W, T, 500)
    if (r$code == "success") { phy <- sim_table_to_phy(r$Ltable, T); if (Ntip(phy) >= 60 && Ntip(phy) <= 150) break } }
  for (k in c(2, 4)) {                       # hide the 0.193 or the 0.393 Ma lowstand
    others <- st[-k, ]; dur <- st$drop[k] - st$rise[k]
    ll <- sapply(grid, function(a) {
      cand <- rbind(others, data.frame(rise = a - dur, drop = a)); cand <- cand[order(cand$rise), ]
      if (a - dur < 0 || any(cand$drop[-nrow(cand)] >= cand$rise[-1])) return(NA)
      envidiv_loglik_fast(pars, phy, cand, sig) })
    ll0 <- envidiv_loglik_fast(pars, phy, others, sig)
    res[[length(res) + 1]] <- data.frame(tree = i, hidden = st$drop[k], a = grid, dll = ll - ll0, ntip = Ntip(phy))
  }
  # rise-time scan for the long 1.1 Ma lowstand (drop fixed at 1.1)
  rgrid <- seq(0.40, 1.08, by = 0.01)
  llr <- sapply(rgrid, function(b) { s2 <- st; s2$rise[5] <- b; envidiv_loglik_fast(pars, phy, s2, sig) })
  res[[length(res) + 1]] <- data.frame(tree = i, hidden = -1, a = rgrid, dll = llr - max(llr), ntip = Ntip(phy))
}
out <- do.call(rbind, res); saveRDS(out, "scan_drop.rds")
for (h in c(0.193, 0.393)) {
  d <- out[out$hidden == h, ]
  pk <- sapply(split(d, d$tree), function(x) x$a[which.max(x$dll)])
  gain <- sapply(split(d, d$tree), function(x) max(x$dll, na.rm = TRUE))
  cat(sprintf("hidden lowstand at %.3f Ma: argmax per tree = %s | max dlogL = %s\n", h,
      paste(round(pk, 3), collapse = " "), paste(round(gain, 1), collapse = " ")))
}
d <- out[out$hidden == -1, ]
for (t in unique(d$tree)) { x <- d[d$tree == t, ]; ci <- range(x$a[x$dll > -1.92])
  cat(sprintf("tree %d: rise of 1.1 Ma stand (true 0.55): MLE %.2f, 95%% profile CI %.2f-%.2f\n", t, x$a[which.max(x$dll)], ci[1], ci[2])) }

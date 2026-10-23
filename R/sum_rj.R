source("rjmcmc_envidiv.R")
for (sc in c("pump", "bd")) {
  f <- sprintf("rj_%s.rds", sc); if (!file.exists(f)) next
  fit <- readRDS(f); tr <- fit$trace; keep <- tr$iter > 0.2 * max(tr$iter)
  cat(sprintf("\n=== %s tree (%d tips), %d iterations in %.0f s\n", sc, ape::Ntip(fit$phy), max(tr$iter), fit$secs))
  q <- function(x) sprintf("%.3f [%.3f, %.3f]", median(x), quantile(x, .025), quantile(x, .975))
  for (v in c("mu", "lambda_H", "lambda_L", "alpha", "sigma", "r")) cat(sprintf("  %-9s %s\n", v, q(tr[keep, v])))
  cat("  J posterior:"); print(round(prop.table(table(tr$J[keep])), 3))
  s <- summarise_water(fit, grid = seq(0.005, fit$T - 0.005, by = 0.0025), window = 0.005)
  saveRDS(s, sprintf("rj_%s_summary.rds", sc))
  for (a in c(0.04, 0.193, 0.295, 0.393, 1.1)) cat(sprintf("  P(drop within 5 kyr of %.3f Ma) = %.2f\n", a, max(s$drop_prob[abs(s$age - a) < 0.003])))
  top <- s[order(-s$drop_prob), ]; sel <- c()
  for (i in seq_len(nrow(top))) if (all(abs(top$age[i] - sel) > 0.02)) sel <- c(sel, top$age[i])
  sel <- sel[1:8]
  cat("  highest drop-intensity ages:", paste(sprintf("%.3f (%.2f)", sel, s$drop_prob[match(sel, s$age)]), collapse = ", "), "\n")
  nd <- ape::branching.times(fit$phy)
  cat("  P(low water) at 0.8 Ma:", round(s$low_prob[which.min(abs(s$age - 0.8))], 2), " at 0.45:", round(s$low_prob[which.min(abs(s$age - 0.45))], 2), " at 2.5:", round(s$low_prob[which.min(abs(s$age - 2.5))], 2), "\n")
}

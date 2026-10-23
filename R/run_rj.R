suppressMessages({library(enviDiv); library(ape)})
source("envidiv_loglik_fast.R"); source("rjmcmc_envidiv.R")
args <- commandArgs(TRUE); scen <- args[1]; n_iter <- as.numeric(args[2]); set.seed(as.integer(args[3]))
T <- 5
truth <- if (scen == "pump") c(0.2, 0.6, 1.0, 1.5) else c(0.2, 0.8, 0.8, 0)
model <- if (scen == "pump") 2 else 1
W <- generate_water(model, T)
repeat { r <- create_tree_cpp(c(truth, 0.02, model), W, T, 500)
  if (r$code == "success") { phy <- sim_table_to_phy(r$Ltable, T); if (Ntip(phy) >= 80 && Ntip(phy) <= 140) break } }
cat("tips:", Ntip(phy), "\n")
t0 <- Sys.time()
fit <- rjmcmc_envidiv(phy, n_iter = n_iter, thin = 50, verbose = TRUE)
fit$phy <- phy; fit$truth <- truth; fit$secs <- as.numeric(Sys.time() - t0, units = "secs")
saveRDS(fit, sprintf("rj_%s.rds", scen))
cat("done in", round(fit$secs), "s\n"); print(round(fit$acceptance, 3))

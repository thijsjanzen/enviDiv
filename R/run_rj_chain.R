suppressMessages({library(ape)})
source("envidiv_loglik_fast.R"); source("rjmcmc_envidiv.R")
args <- commandArgs(TRUE); set.seed(as.integer(args[1]))
base <- readRDS("rj_pump.rds")
init <- list(pars = c(mu = 0.5, lambda_H = 0.2, lambda_L = 2, alpha = 0.1), sigma = 0.1, r = 10,
             stands = data.frame(rise = c(1.5, 3.0), drop = c(2.0, 3.6)))   # deliberately poor start
fit <- rjmcmc_envidiv(base$phy, n_iter = 150000, thin = 50, init = init, verbose = FALSE)
fit$phy <- base$phy; saveRDS(fit, sprintf("rj_pump_chain%s.rds", args[1]))

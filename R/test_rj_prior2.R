suppressMessages({library(ape)})
source("envidiv_loglik_fast.R"); source("rjmcmc_envidiv.R")
d <- readRDS("tree_timing.rds"); phy <- d$phy
args <- commandArgs(TRUE); bu <- as.numeric(args[1]); set.seed(as.integer(args[2]))
r <- 1
fit <- rjmcmc_envidiv(phy, n_iter = 4e5, thin = 10, prior_only = TRUE, fix = c("rates", "sigma", "r"), verbose = FALSE,
  control = list(birth_unif = bu),
  init = list(pars = c(mu = .1, lambda_H = .5, lambda_L = .5, alpha = .5), sigma = .02, r = r, stands = data.frame(rise = numeric(0), drop = numeric(0))))
J <- fit$trace$J[-(1:2000)]; jj <- 0:5
p <- dpois(2 * jj, r * fit$T) / ((1 + exp(-2 * r * fit$T)) / 2)
b <- split(J, cut(seq_along(J), 20))   # batch means
m <- sapply(jj, function(j) mean(J == j)); se <- sapply(jj, function(j) sd(sapply(b, function(x) mean(x == j))) / sqrt(20))
cat(sprintf("birth_unif = %.1f\n", bu)); print(round(rbind(analytic = p, mcmc = m, z = (m - p) / se), 3))

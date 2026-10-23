source("rjmcmc_envidiv.R")
fs <- c("rj_pump.rds", list.files(pattern = "^rj_pump_chain.*rds$"))
grid <- seq(0.005, 4.995, by = 0.0025)
res <- lapply(fs, function(f) { fit <- readRDS(f); tr <- fit$trace; keep <- tr$iter > 0.2 * max(tr$iter)
  s <- summarise_water(fit, grid = grid, window = 0.005)
  c(J = median(tr$J[keep]), alpha = median(tr$alpha[keep]), lambda_L = median(tr$lambda_L[keep]), sigma = median(tr$sigma[keep]),
    d0.193 = max(s$drop_prob[abs(s$age - .193) < .003]), d0.295 = max(s$drop_prob[abs(s$age - .295) < .003]),
    d1.083 = max(s$drop_prob[abs(s$age - 1.083) < .003]), d0.427 = max(s$drop_prob[abs(s$age - .427) < .003]),
    low0.8 = s$low_prob[which.min(abs(s$age - .8))], low2.5 = s$low_prob[which.min(abs(s$age - 2.5))]) })
out <- do.call(rbind, res); rownames(out) <- fs; print(round(out, 2))
# pooled summary for plotting
fits <- lapply(fs, readRDS); pooled <- fits[[1]]
pooled$stands <- unlist(lapply(fits, function(f) f$stands[-(1:ceiling(0.2 * length(f$stands)))]), recursive = FALSE)
s <- summarise_water(pooled, burnin = 0, grid = seq(0.0025, 4.9975, by = 0.005), window = 0.0025)
saveRDS(s, "rj_pump_pooled_summary.rds")



while (TRUE) {
  #ref_tree <- sim_func(params = c(0.01, 0.8, 0.1, 0.4, 0.0, 10, 4))
  ref_tree <- sim_func(params = c(0.0, 0.8, 0.1, 0.4, 0.0, 10, 1))
  if (inherits(ref_tree$phy, "phylo")) {
    num_lin <- length(ref_tree$phy$tip.label)
    if (num_lin > 50 && num_lin < 100) {
      break
    }
  }
}

ref_water <- ref_tree$waterlevel
ref_tree <- ref_tree$phy

ca <- treestats::crown_age(ref_tree)
num_lin <- treestats::number_of_lineages(ref_tree)

ca
num_lin


prior_func <- function() {
  vv <- enviDiv::param_from_prior_cpp()
  return(vv)
}

prior_dens_func <- function(params) {
  for (i in 1:6) {
    if (params[i] < 0) return(0)
    x <- log10(params[i])
    if (x < -3 || x > 5) return(0)
  }
  if (params[7] < 1 || params[7] > 4) return(0)

  return(1)
}

test_pars <- prior_func()
test_tree <- sim_func(test_pars)

#stat_func <- function(focal_tree) {
#  res <- treestats::calc_all_stats(focal_tree)
#  index <- which(names(res) == "rquartet")
#  res <- res[-index]
#  index <- which(names(res) == "wiener")
#  res <- res[-index]
#  return(res)
#}

stat_func <- create_statistics_list()
#  t0 <- Sys.time()
  res <- enviDiv::abc_smc_par(ref_tree = ref_tree,
                          statistics = stat_func,
                          simulation_function = sim_func,
                          init_epsilon_value = 1000000,
                          prior_generating_function = prior_func,
                          prior_density_function = prior_dens_func,
                          number_of_particles = 1000,
                          sigma = 0.05,
                          stop_rate = 1e-6,
                          num_iterations = 8,
                          num_threads = 6)
#  t1 <- Sys.time()
#  to_add <- c(nt, difftime(t1, t0, units = "secs")[[1]])
#  cat(to_add , "\n")
#  found <- rbind(found, to_add)
#}



to_plot <- c()
for (r in 1:length(res$all_parameters)) {
  focal_iter <- res$all_parameters[[r]]

  focal_iter <- cbind(focal_iter, r)
  to_plot <- rbind(to_plot, focal_iter)
}

colnames(to_plot) <- c("extinction", "symp_spec_high", "symp_spec_low",
                       "allo_spec", "jiggle", "water", "model", "repl")
require(tidyverse)
to_plot <- as_tibble(to_plot)

to_plot %>%
  gather(key = "parameter", value = "val", -c(repl, model)) %>%
  ggplot(aes(x = repl, y = val, group = interaction(repl,model), fill = as.factor(model))) +
    geom_boxplot() +
    scale_y_log10() +
    facet_wrap(~parameter, scales = "free")

to_plot %>%
  ggplot(aes(x = model)) +
    geom_bar() +
    facet_grid(rows=vars(repl))

to_plot %>%
  filter(repl == max(repl)) %>%
  gather(key = "parameter", value = "val", -c(repl, model)) %>%
  ggplot(aes(x = val)) +
    geom_density(bw = "sj") +
    scale_x_log10() +
    facet_wrap(~parameter, scales = "free")

to_plot %>%
  filter(repl == max(repl)) %>%
  ggplot(aes(x = model)) +
  geom_bar() +
  facet_grid(rows=vars(repl))


wlvls <- c()
for (r in 2:length(res$all_parameters)) {
  focal_iter <- res$all_water[[r]]
  focal_model <- res$all_parameters[[r]][, 7]

  for (i in 1:length(focal_iter)) {
     ws <- focal_iter[[i]]
     to_add <- cbind(focal_model[i], ws, r)
     wlvls <- rbind(wlvls, to_add)
  }
}

colnames(wlvls) <- c("model", "time", "repl")
wlvls <- as_tibble(wlvls)

wlvls %>%
  group_by(model, repl) %>%
  ggplot(aes(x = time, col = as.factor(repl), group = repl)) +
    geom_density(bw = "sj") +
    facet_grid(rows = vars(model), scales = "free")

wlvls %>%
  group_by(model, repl) %>%
  filter(model == 3) %>%
  ggplot(aes(x = time, col = as.factor(repl), group = repl)) +
  geom_density(bw = "sj")


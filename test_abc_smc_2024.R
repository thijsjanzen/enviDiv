#ref_tree <- ape::read.tree("/Users/thijsjanzen/MEGAsync2/Leonel/attempt1/lamprologini_walter.tree")

ca <- 5

sim_func <- function(params) {
  sim_tree <- enviDiv::sim_envidiv_cpp(params[7],
                                       params[1:6],
                                       ca,
                                       500)
  if (sim_tree$code != "failure") {
    phy <- ape::read.tree(text = sim_tree$code)
    if (length(phy$tip.label) < 4) return(list(phy = "failure"))
  } else {
    return(list("phy" = sim_tree$code))
  }
  return(list("phy" = phy,
              "waterlevel" = sim_tree$water))
}

while (TRUE) {
  ref_tree <- sim_func(params = c(0.01, 0.8, 0.1, 0.4, 0.0, 10, 4))
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
    x <- log10(params[i])
    if (x < -3 || x > 5) return(0)
  }
  if (params[7] < 1 || params[7] > 3) return(0)

  return(1)
}

test_pars <- prior_func()
test_tree <- sim_func(test_pars)

stat_func <- function(focal_tree) {
  res <- treestats::calc_all_stats(focal_tree)
  index <- which(names(res) == "rquartet")
  res <- res[-index]
  index <- which(names(res) == "wiener")
  res <- res[-index]
  return(res)
}


res <- enviDiv::abc_smc(ref_tree = ref_tree,
                        statistics = stat_func,
                        simulation_function = sim_func,
                        init_epsilon_value = 10000,
                        prior_generating_function = prior_func,
                        prior_density_function = prior_dens_func,
                        number_of_particles = 1000,
                        sigma = 0.05,
                        stop_rate = 1e-3,
                        num_iterations = 3)

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
    facet_grid(rows = vars(model))
xlim(0, 5)

wlvls %>%
  group_by(model, repl) %>%
  filter(model == 3) %>%
  ggplot(aes(x = time, col = as.factor(repl), group = repl)) +
  geom_density(bw = "sj")


ref_tree = ref_tree
statistics = treestats::calc_all_stats
simulation_function = sim_func
init_epsilon_values = init_epsilon_values
prior_generating_function = prior_func
prior_density_function = prior_dens_func
number_of_particles = 100
sigma = 0.05
stop_rate = 1e-5


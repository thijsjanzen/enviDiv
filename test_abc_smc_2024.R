create_statistics_list <- function() {

  stats <- list()
  stats$rquartet <- treestats::rquartet # moved up to trigger early fail
  stats$colless_quad <- treestats::colless_quad # moved up to trigger early fail
  stats$wiener <- treestats::wiener # moved up to trigger early fail
  stats$gamma <- treestats::gamma_statistic
  stats$sackin <- treestats::sackin
  stats$colless <- treestats::colless
  stats$colless_corr <- treestats::colless_corr
  stats$beta <- treestats::beta_statistic
  stats$blum <- treestats::blum
  stats$pigot_rho <- treestats::pigot_rho

  stats$treeness <- treestats::treeness
  stats$nltt_base <- treestats::nLTT_base
  stats$phylogenetic_div <- treestats::phylogenetic_diversity
  stats$avg_ladder <- treestats::avg_ladder
  stats$max_ladder <- treestats::max_ladder
  stats$cherries <- treestats::cherries
  stats$double_cherries <- treestats::double_cherries
  stats$four_prong <- treestats::four_prong
  stats$il_number <- treestats::ILnumber
  stats$pitchforks <- treestats::pitchforks

  stats$stairs <- treestats::stairs
  stats$imbalance_steps <- treestats::imbalance_steps
  stats$j_one <- treestats::j_one
  stats$b1 <- treestats::b1
  stats$b2 <- treestats::b2
  stats$area_per_pair <- treestats::area_per_pair
  stats$average_leaf_depth <- treestats::average_leaf_depth
  stats$i_stat <- treestats::mean_i
  stats$ew_colless <- treestats::ew_colless
  stats$max_del_width <- treestats::max_del_width

  stats$max_depth <- treestats::max_depth
  stats$avg_vert_depth <- treestats::avg_vert_depth
  stats$max_width <- treestats::max_width
  stats$mw_over_md <- treestats::mw_over_md
  stats$tot_path <- treestats::tot_path_length
  stats$tot_internal_path <- treestats::tot_internal_path
  stats$rogers <- treestats::rogers
  stats$stairs2 <- treestats::stairs2
  stats$tot_coph <- treestats::tot_coph
  stats$var_depth <- treestats::var_leaf_depth

  stats$symmetry_nodes <- treestats::sym_nodes
  stats$mpd <- treestats::mean_pair_dist
  stats$psv <- treestats::psv
  stats$vpd <- treestats::var_pair_dist
  stats$mntd <- treestats::mntd
  stats$j_stat <- treestats::entropy_j
  stats$crown_age <- treestats::crown_age
  stats$tree_height <- treestats::tree_height
  stats$max_betweenness <- treestats::max_betweenness
  stats$diameter <- treestats::diameter

  local_closeness <- function(tree, w, n) {
    return(treestats::max_closeness(tree, weight = w, normalization = ifelse(n ==
                                                                               TRUE, "tips", "none")))
  }
  stats$max_closeness <- function(x) {
    local_closeness(x, FALSE, FALSE)
  }
  stats$max_closenessW <- function(x) {
    local_closeness(x, TRUE, FALSE)
  }

  stats$eigen_centrality <- function(x) {
    return(max(treestats::eigen_centrality(x, weight = FALSE)$eigenvector))
  }
  stats$eigen_centralityW <- function(x) {
    return(max(treestats::eigen_centrality(x, weight = TRUE)$eigenvector))
  }
  stats$mean_branch_length <- treestats::mean_branch_length
  stats$var_branch_length <- treestats::var_branch_length
  stats$mean_branch_length_int <- treestats::mean_branch_length_int
  stats$mean_branch_length_ext <- treestats::mean_branch_length_ext
  stats$var_branch_length_int <- treestats::var_branch_length_int
  stats$var_branch_length_ext <- treestats::var_branch_length_ext


  stats$root_imbalance <- treestats::root_imbalance
  stats$number_of_lineages <- treestats::number_of_lineages

  get_minmax_lapl <- function(phylo) {
    out <- list()
    temp_stats <- treestats::minmax_laplace(phylo, TRUE)

    if (length(temp_stats) >= 2) {
      out$min_laplace <- temp_stats$min
      out$max_laplace <- temp_stats$max
    }
    else {
      out$min_laplace <- NA
      out$max_laplace <- NA
    }
    return(out)
  }
  stats$minmax_lapl <- get_minmax_lapl

  get_minmax_adj <- function(phylo) {
    out <- list()
    temp_stats <- treestats::minmax_adj(phylo, TRUE)

    if (length(temp_stats) >= 2) {
      out$min_adj <- temp_stats$min
      out$max_adj <- temp_stats$max
    }
    else {
      out$min_adj <- NA
      out$max_adj <- NA
    }
    return(out)
  }
  stats$minmax_adj <- get_minmax_adj

  get_laplacian_dist <- function(phylo) {
    temp_stats <- tryCatch(expr = {
      treestats::laplacian_spectrum(phylo)
    }, error = function(e) {
      return(NA)
    })

    out <- list()
    if (length(temp_stats) == 5) {
      out$laplace_spectrum_a <- temp_stats$asymmetry
      out$laplace_spectrum_p <- temp_stats$peakedness
      out$laplace_spectrum_e <- log(temp_stats$principal_eigenvalue)
      out$laplace_spectrum_g <- temp_stats$eigengap[[1]]
    }
    else {
      out$laplace_spectrum_a <- NA
      out$laplace_spectrum_p <- NA
      out$laplace_spectrum_e <- NA
      out$laplace_spectrum_g <- NA
    }
    return(out)
  }
  stats$laplacian_d <- get_laplacian_dist

  return(stats)
}

ref_tree <- ape::read.tree("/Users/thijsjanzen/MEGAsync2/Leonel/attempt1/lamprologini_walter.tree")

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


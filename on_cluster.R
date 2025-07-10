sim_func <- function(params) {
  sim_tree <- enviDiv::sim_envidiv_tree_new(params,
                                              model = params[7],
                                              crown_age =  6.170882,
                                              max_lin = 500)
  if (is.null(sim_tree$phy)) {
    return(list("phy" = sim_tree$code))
  }
  return(list("phy" = sim_tree$phy,
              "waterlevel" = sim_tree$water))
}


param_grid <- expand.grid(model = 1:3,
                          repl = 1:100)

args = commandArgs(trailingOnly = TRUE)

arg_number <- as.numeric(args[[1]])

focal_model <- param_grid$model[[arg_number]]

while (TRUE) {
  sim_params <- enviDiv::param_from_prior_cpp(model = focal_model)
  ref_tree <- sim_func(params = sim_params)
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
    if (params[i] < 0) return(-Inf)
    x <- log10(params[i])
    if (x < -3 || x > 5) return(-Inf)
  }
  if (params[7] < 1 || params[7] > 3) return(-Inf)

  return(1)
}

test_pars <- prior_func()
test_tree <- sim_func(test_pars)

stat_func <- create_statistics_list()

res <- enviDiv::abc_smc_par(ref_tree = ref_tree,
                            statistics = stat_func,
                            simulation_function = sim_func,
                            init_epsilon_value = 100000,
                            prior_generating_function = prior_func,
                            prior_density_function = prior_dens_func,
                            number_of_particles = 1000,
                            sigma = 0.01,
                            stop_rate = 1e-10,
                            num_iterations = 3,
                            num_threads = 8,
                            write_to_file = TRUE,
                            file_name_start = "out_")



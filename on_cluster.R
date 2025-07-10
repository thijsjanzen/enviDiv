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

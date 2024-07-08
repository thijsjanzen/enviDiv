v <- enviDiv::initial_draw_from_prior(num_particles = 1000,
                                      crown_age = 5,
                                      min_lin = 5,
                                      max_lin = 500,
                                      verbose = TRUE)

sim_func <- function(params, ca) {
  sim_tree <- enviDiv::sim_envidiv_cpp(params[6],
                                       params[1:5],
                                       ca,
                                       500)
  if (sim_tree != "failure") {
    phy <- ape::read.tree(text = sim_tree)
    if (length(phy$tip.label) < 4) return("failure")
  } else {
    return(sim_tree)
  }
  return(phy)
}



ca <- 5
focal_params <- c(0.0, 0.0, 0.5, 0, 0.0, 4)

num_repl <- 1000
num_trees <- 0
found <- c()
pb <- txtProgressBar(max = num_repl, style = 3)
for (r in 1:num_repl) {
  focal_tree <- sim_func(focal_params, ca)
  if (inherits(focal_tree, "phylo")) {
    brts <- treestats::branching_times(focal_tree)
    found <- c(found, brts)
    num_trees <- num_trees + 1
  }
}

plot(density(found, bw = "sj"), xlim = c(0, ca))
#hist(found)

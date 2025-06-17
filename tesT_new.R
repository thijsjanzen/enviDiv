extinct_rate <- 0
sym_high <- 0.3
sym_low <- 0
allo_low <- 0
wiggle <- 0
model <- 1

params <- c(extinct_rate, sym_high, sym_low, allo_low)

crown_age <- 4

found_tree <- enviDiv::sim_envidiv_tree_new(params = params,
                                            model = 1,
                                            crown_age = crown_age,
                                            max_lin = 500,
                                            seed = 1)
cat(length(found_tree$ltable[, 1]), "\n")

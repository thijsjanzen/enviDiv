extinct_rate <- 0
sym_high <- 0.5
sym_low <- 0
allo_low <- 0
wiggle <- 0
model <- 1

params <- c(extinct_rate, sym_high, sym_low, allo_low, wiggle, model)

crown_age <- 3

found_tree <- enviDiv::sim_envidiv_tree_new(params = params,
                                            crown_age = crown_age,
                                            max_lin = 500)
plot(found_tree$phy)
found_tree$ltable
tab2 <- correct_ltable_internal(found_tree$ltable, crown_age)

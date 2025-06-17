extinct_rate <- 0
sym_high <- 0.5
sym_low <- 0
allo_low <- 0
wiggle <- 0
model <- 1
water_rate <- 3.1

params <- c(extinct_rate, sym_high, sym_low, allo_low, wiggle, water_rate)

crown_age <- 5

found_tree <- enviDiv::sim_envidiv_tree_new(params = params,
                                            model = 3,
                                            crown_age = crown_age,
                                            max_lin = 500,
                                            seed = 2)
plot(found_tree$phy)
found_tree$water

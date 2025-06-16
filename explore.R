extinct_rate <- 0
sym_high <- 0.5
sym_low <- 0
allo_low <- 0
wiggle <- 0
model <- 3

params <- c(extinct_rate, sym_high, sym_low, allo_low, wiggle, model)

crown_age <- 3

found_tree <- enviDiv::sim_envidiv_tree2(params, crown_age)

plot(found_tree$phy, type = "cladogram", direction = "up")
library(hybridBD)
library(tidyverse)
corrected_ltable <- correct_ltable(found_tree$ltable, crown_age)
plot(treestats::l_to_phylo(corrected_ltable), type = "cladogram", direction = "up")

hybridBD::plot_ltable(corrected_ltable)


found_tree$water

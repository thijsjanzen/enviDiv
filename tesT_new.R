extinct_rate <- 0
sym_high <- 0.5
sym_low <- 0
allo_low <- 0
wiggle <- 0
model <- 1
water_rate <- 3.3

params <- c(extinct_rate, sym_high, sym_low, allo_low, wiggle, water_rate)

crown_age <- 5
for (s in 1:100) {
#s <- 99
  #  cat(s, "\n")
found_tree <- enviDiv::sim_envidiv_tree_new(params = params,
                                            model = 3,
                                            crown_age = crown_age,
                                            max_lin = 500,
                                            seed = s)
    if(!is.null(found_tree$phy)) {
      plot(found_tree$phy, main = s)
        ltab <- found_tree$ltable
        vv <- diff(ltab[, 1])
        vv <- vv[2:length(vv)]
        if (length(which(vv > 0)))
          cat(s, "\n")
    }

}

found_tree$water
ltab <- found_tree$ltable
ltab[, 1] <- crown_age - ltab[, 1]
ltab
plot(found_tree$phy)
cat("done\n")

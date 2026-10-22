ref_tree <- ape::read.tree("/Users/thijsjanzen/MEGAsync2/Leonel/attempt1/lamprologini_walter.tree")
ca <- treestats::crown_age(ref_tree)
require(tidyverse)
found_stack <- c()
require(enviDiv)
for (model in c(1, 2, 4)) {
  local_stack <- generate_stack_simple(number_of_replicates = 100000,
                                       focal_model = model,
                                       min_tips = 80,
                                       max_tips = 80,
                                       crown_age = ca,
                                       num_threads = 8)
  found_stack <- rbind(found_stack, local_stack)
}

write_tsv(found_stack, file = "focal_stack_100k_nowobble.txt")

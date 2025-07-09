

for (model in 1:3) {
  enviDiv::generate_stack_simple(number_of_replicates = 100000,
                 focal_model = model,
                 min_tips = 75,
                 max_tips = 85,
                 crown_age = 6.2,
                 write_to_file = TRUE,
                 file_name = paste0("res_restrict_", model, ".txt"),
                 num_threads = 8)
}

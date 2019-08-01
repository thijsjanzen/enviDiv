context("infer_params")

test_that("usage", {
  #skip_on_cran()

  params <- c(0, 0.5, 0, 0, 0, 1)

  emp_tree <- sim_envidiv_tree(params, crown_age = 5, seed = 5)

  do_smc <- infer_params(number_of_particles = 1000,
                                  max_iter = 15,
                                  sd_params = 0.01,
                                  emp_tree,
                                  write_to_file = FALSE)

  v <- mean(do_smc$sym_high)
  v2 <- table(do_smc$model)
  testthat::expect_true(which.max(v2) == 1)
  q1 <- quantile(do_smc$sym_high, c(0.025, 0.975))
  testthat::expect_true(q1[[1]] < 0.5 && q1[[2]] > 0.5)
})

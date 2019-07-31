context("calc_fit")

test_that("usage", {
  params <- c(0, 0.5, 0, 0, 0, 1)

  emp_tree <- sim_envidiv_tree(params, crown_age = 5)
  while(is.null(emp_tree)) emp_tree <- sim_envidiv_tree(params, crown_age = 5)

  v1 <- calc_sum_stats(emp_tree, emp_tree)

  focal_fit <- calc_fit(v1[1:4], v1[1:4])
  testthat::expect_equal(focal_fit, 0)
})

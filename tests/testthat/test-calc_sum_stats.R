context("calculate sum stats")

test_that("usage", {
  params <- c(0, 0.5, 0, 0, 0, 1)
  emp_tree <- sim_envidiv_tree(params, crown_age = 5, seed = 42)

  while (is.null(emp_tree))
    emp_tree <- sim_envidiv_tree(params, crown_age = 5)

  emp_tree2 <- sim_envidiv_tree(params, crown_age = 5, seed = 43)

  while (is.null(emp_tree2))
    emp_tree2 <- sim_envidiv_tree(params, crown_age = 5)

  v1 <- calc_sum_stats(emp_tree)
  v2 <- calc_sum_stats(emp_tree, emp_tree2)

  a <- sum(v1 != v2)
  testthat::expect_true(a == 1)
})

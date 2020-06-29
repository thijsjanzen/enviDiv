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
  v2 <- calc_sum_stats(emp_tree2)

  a <- sum(v1 != v2)
  testthat::expect_true(a != 0)
  testthat::expect_equal(v1[[11]], 5)
  testthat::expect_equal(v2[[11]], 5)
})

test_that("abuse", {
  params <- c(0, 0.5, 0, 0, 0, 1)
  phy <- NULL
  testthat::expect_warning(calc_sum_stats(phy))
  testthat::expect_warning(calc_sum_stats(5))

  params <- c(0, 0.01, 0, 0, 0, 1)
  emp_tree <- sim_envidiv_tree(params, crown_age = 5, seed = 1)
  vx <- calc_sum_stats(emp_tree)
  testthat::expect_true(is.nan(vx[2]))
  testthat::expect_true(is.infinite(vx[6]))
  testthat::expect_true(is.infinite(vx[7]))
})

context("generate_naive")

test_that("usage", {

  params <- c(0, 0.5, 0, 0, 0, 1)
  emp_tree <- sim_envidiv_tree(params, crown_age = 5, seed = 5)

  # skip
  if(1 == 2) {
    vx <- generate_naive(number_of_particle = 10,
                         min_tips = 50,
                         max_tips = 150,
                         emp_tree = emp_tree)

    testthat::expect_true(length(vx$extinct) >= 10)
    testthat::expect_true(min(vx$num_lin) >= 50)
    testthat::expect_true(max(vx$num_lin) <= 150)
  }
})

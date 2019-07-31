context("mutate params")

test_that("use", {
  set.seed(42)
  params <- param_from_prior()

  new_params <- mutate_params(params, local_sd = 0)

  testthat::expect_equal(new_params, params)

  other_params <- mutate_params(params, local_sd = 0.1)

  # if they are identical, the sum of != is equal to 0
  testthat::expect_gt(sum(other_params != params), 0)
  testthat::expect_equal(sum(new_params != params), 0)

  for(r in 1:100) {
    new_params <- mutate_params(params, local_sd = 0.1)
  }
})

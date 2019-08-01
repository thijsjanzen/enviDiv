context("is within prior")

test_that("use", {
  params <- param_from_prior()
  testthat::expect_true(is_within_prior(params))

  bogus_params <- c(1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 2, 1)
  testthat::expect_false(is_within_prior(bogus_params))

  params <- param_from_prior()
  test_params <- params
  test_params[6] <- 5
  testthat::expect_false(is_within_prior(test_params))

  test_params <- params
  test_params[6] <- 0
  testthat::expect_false(is_within_prior(test_params))

  for (i in 5:1) {
    test_params <- params
    test_params[i] <- 10 ^ -4
    testthat::expect_false(is_within_prior(test_params))

    test_params <- params
    test_params[i] <- 10 ^ 3
    testthat::expect_false(is_within_prior(test_params))
  }
})

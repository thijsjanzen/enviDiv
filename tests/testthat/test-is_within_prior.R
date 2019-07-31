context("is within prior")

test_that("use", {
  params <- param_from_prior()
  testthat::expect_true(is_within_prior(params))

  bogus_params <- c(1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 2, 1)
  testthat::expect_false(is_within_prior(bogus_params))
})

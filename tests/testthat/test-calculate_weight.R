context("calculate weight")

test_that("usage", {
  param1 <- param_from_prior()

  bogus_previous_params <- rbind(param1, param1)

  new_params1 <- mutate_params(param1, local_sd = 0.5)
  new_weight1 <- calculate_weight(new_params1, bogus_previous_params,
                                 c(1, 1)/2, sigma = 0.1)

  new_weight2 <-  calculate_weight(param1, bogus_previous_params,
                                   c(1, 1)/2, sigma = 0.1)


  # less weight is added to particles that are easy to obtain:
  testthat::expect_lt(new_weight2, new_weight1)
})

test_that("abuse", {
  param1 <- param_from_prior()

  bogus_previous_params <- rbind(param1, param1)

  testthat::expect_error(
    new_weight1 <- calculate_weight(new_params1, bogus_previous_params,
                                  c(1, 1), sigma = 0.1)
  )

  new_weight1 <- calculate_weight(new_params1, bogus_previous_params,
                                  c(1, 1) / 2, sigma = 0.0)

  testthat::expect_true(is.infinite(new_weight1))
})

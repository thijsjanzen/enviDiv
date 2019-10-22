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

  for (r in 1:100) {
    new_params <- mutate_params(new_params, local_sd = 0.1)
  }

  param_table <- matrix(nrow = 1e4, ncol = 7)
  for (i in 1:1e4) param_table[i, ] <- param_from_prior()


  for (model in 1:3) {
      param_table[, 6] <- model
      param_table2 <- matrix(nrow = 1e4, ncol = 7)
      for (i in 1:1e4) param_table2[i, ] <- mutate_params(param_table[i, ],
                                                          local_sd = 0.1)
      vx <- table(param_table2[, 6])
      vx <- vx / sum(vx)
      a <- vx[model]
      b <- vx[-model]

      testthat::expect_equal(a[[1]], 0.9, tolerance = 0.1)
      testthat::expect_equal(b[[1]], 0.05, tolerance = 0.1)
      testthat::expect_equal(b[[2]], 0.05, tolerance = 0.1)
  }


})

context("get from prior")

test_that("get_from_prior", {

  for (r in 1:100) {
    v1 <- enviDiv::param_from_prior()
    extinct <- v1[1]
    sym_high <- v1[2]
    sym_low <- v1[3]
    allo <- v1[4]
    jiggle <- v1[5]
    model <- v1[6]

    for (i in 1:6) {
      testthat::expect_gt(v1[i], 0)
    }

    testthat::expect_true(model %in% 1:3)
  }
})

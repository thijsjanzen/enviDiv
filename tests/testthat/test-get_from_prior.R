context("get from prior")

test_that("get_from_prior", {

  for(r in 1:100) {
    v1 <- enviDiv::param_from_prior()
    extinct <- v1[1]
    sym_high <- v1[2]
    sym_low <- v1[3]
    allo <- v1[4]
    jiggle <- v1[5]
    model <- v1[6]

    testthat::expect_gt(max(v1), 0)

    testthat::expect_true((extinct) > 0)
    testthat::expect_true(log10(extinct) < 2)

    testthat::expect_true((sym_high) > 0)
    testthat::expect_true(log10(sym_high) < 2)

    testthat::expect_true((sym_low) > 0)
    testthat::expect_true(log10(sym_low) < 2)

    testthat::expect_true((allo) > 0)
    testthat::expect_true(log10(allo) < 2)

    testthat::expect_true((jiggle) > 0)
    testthat::expect_true(log10(jiggle) < 1)

    testthat::expect_true(model %in% 1:3)
  }
})

context("generate water")

test_that("generate water", {

  crown_age <- 5

  w1 <- generate_water(1, crown_age)
  w2 <- generate_water(2, crown_age)
  w3 <- generate_water(3, crown_age)

  testthat::expect_true(length(w1) < length(w2))
  testthat::expect_true(length(w2) < length(w3))

  testthat::expect_lte(max(w1), 2 * crown_age)
  testthat::expect_lte(max(w2), crown_age)
  testthat::expect_lte(max(w3), crown_age)
})

test_that("abuse", {

  crown_age <- -1

  testthat::expect_warning(
    w1 <- generate_water(1, crown_age),
    "maximum time has to be > 0")
  testthat::expect_warning(
    w2 <- generate_water(2, crown_age),
    "maximum time has to be > 0")

  testthat::expect_warning(
    w3 <- generate_water(3, crown_age),
    "maximum time has to be > 0")
})

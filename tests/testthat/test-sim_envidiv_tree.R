context("sim envidiv tree ")

test_that("usage", {

  params <- c(0, 1, 0, 0, 0, 1)

  crown_age <- 5

  t1 <- sim_envidiv_tree(params, crown_age, abc = FALSE, seed = 1)
  testthat::expect_true(class(t1) == "phylo")

  # without seed:
  t1 <- sim_envidiv_tree(params, crown_age, abc = FALSE)
})

test_that("other models", {
  params <- c(0, 0.5, 1, 1, 0.1, 2)

  crown_age <- 5

  t1 <- sim_envidiv_tree(params, crown_age, abc = FALSE, seed = 1)

  params <- c(0, 0.5, 0.1, 0.1, 0.1, 3)

  crown_age <- 5

  t1 <- sim_envidiv_tree(params, crown_age, abc = FALSE, seed = 1)
})

test_that("abuse", {
  params <- c(0, 0.5, 0, 0, 0, 1)

  crown_age <- -5
  testthat::expect_warning(
    t1 <- sim_envidiv_tree(params, crown_age, abc = FALSE, seed = 1),
    "crown age should be larger than zero"
  )

  params <- c(10, 0.5, 0, 0, 0, 1)
  crown_age <- 5
  testthat::expect_warning(
    t1 <- sim_envidiv_tree(params, crown_age, abc = FALSE, seed = 1),
    "Tree went extinct, returning NULL"
  )

  params <- c(0, 10, 0, 0, 0, 1)
  crown_age <- 5
  testthat::expect_warning(
    t1 <- sim_envidiv_tree(params, crown_age, abc = FALSE, seed = 1),
    "Tree too big, returning NULL"
  )

  params <- c(0, 0.1, 0, 0, 0, 1)
  crown_age <- 5
  testthat::expect_warning(
    t1 <- sim_envidiv_tree(params, crown_age, abc = FALSE, seed = 1),
    "tree has only two tips"
  )
  testthat::expect_silent(
    t1 <- sim_envidiv_tree(params, crown_age, abc = TRUE, seed = 1)
  )
})

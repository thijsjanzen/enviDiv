context("generate_tbb")

test_that("usage", {

  v1 <- generate_trees_tbb(number_of_trees = 100,
                           min_tips = 2,
                           max_tips = 200,
                           model = 1,
                           crown_age = 4,
                           num_threads = 1)

  testthat::expect_true(length(v1$extinct) >= 100)
  testthat::expect_true(min(v1$num_lin) >= 2)
  testthat::expect_true(max(v1$num_lin) <= 200)

  v2 <- generate_trees_tbb(number_of_trees = 100,
                           min_tips = 2,
                           max_tips = 200,
                           model = -1,
                           crown_age = 4,
                           num_threads = 3,
                           write_to_file = FALSE)

  testthat::expect_true(length(v1$extinct) >= 100)
  testthat::expect_true(min(v1$num_lin) >= 2)
  testthat::expect_true(max(v1$num_lin) <= 200)
})

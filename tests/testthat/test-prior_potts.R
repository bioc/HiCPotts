test_that("gamma_prior_value returns a finite scalar in [0, 1]", {
  x <- replicate(100, gamma_prior_value())
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0 & x <= 1))
  expect_gt(length(unique(x)), 1L)       # stochastic
})

test_that("size_prior rejects invalid inputs", {
  expect_error(size_prior(-1, 1),  "positive")
  expect_error(size_prior(NA, 1),  "finite")
  expect_error(size_prior(1, 4),   "component")
  expect_true(is.finite(size_prior(3.0, 2)))
})
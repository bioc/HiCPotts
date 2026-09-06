test_that("size_prior computes the correct log-probability for valid inputs", {
  # Test case 1: Valid input for component 1
  result <- size_prior(size_value = 10, component = 1)
  expected <- dgamma(10, shape = 3, rate = 1, log = TRUE)
  expect_equal(result, expected)

  # Test case 2: Valid input for component 2
  result <- size_prior(size_value = 5, component = 2)
  expected <- dgamma(5, shape = 2, rate = 0.2, log = TRUE)
  expect_equal(result, expected)

  # Test case 3: Valid input for component 3
  result <- size_prior(size_value = 3, component = 3)
  expected <- dgamma(3, shape = 2, rate = 0.2, log = TRUE)
  expect_equal(result, expected)
})


# Test invalid component
test_that("size_prior handles invalid inputs gracefully", {
  # Test invalid component less than 1
  expect_error(size_prior(10, 0), "Invalid component specified")

  # Test invalid component greater than 3
  expect_error(size_prior(10, 4), "Invalid component specified")
})

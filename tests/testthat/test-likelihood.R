make_pred <- function(log_mu) {
  force(log_mu)
  function(params, z, x_vars, component, N) {
    rep(log_mu, sum(z == component))
  }
}

test_that("Poisson likelihood (component 1) equals sum(dpois(..., log=TRUE))", {
  N <- 3
  z <- matrix(c(1,1,2,3,1,2,3,1,2), 3, 3, byrow = TRUE)
  y <- matrix(c(0,3,5,1,0,4,6,2,0), 3, 3, byrow = TRUE)
  params <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x_vars <- replicate(4, list(matrix(1, N, N)), simplify = FALSE)
  
  log_mu <- 1.5
  pred <- make_pred(log_mu)
  yc  <- y[z == 1]
  expected <- sum(dpois(yc, lambda = exp(log_mu), log = TRUE))
  
  got <- likelihood_combined(
    pred_combined = pred, params = params, z = z, y = y,
    x_vars = x_vars, component = 1,
    theta = 0, size = 1, N = N, dist = "Poisson"
  )
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("ZIP zero branch uses exp(-mu), not exp(-eta) (regression test)", {
  N <- 2
  z <- matrix(1, N, N)               # all component 1
  y <- matrix(c(0, 0, 0, 0), N, N)   # all zeros
  params <- rep(0, 5)
  x_vars <- replicate(4, list(matrix(1, N, N)), simplify = FALSE)
  
  log_mu <- 2; mu <- exp(log_mu); theta <- 0.3
  pred <- make_pred(log_mu)
  
  ## Every cell in the complete N x N matrix contributes once.
  pair_w <- rep(1, 4)
  expected <- sum(pair_w * rep(log(theta + (1 - theta) * exp(-mu)), 4))
  wrong    <- sum(pair_w * rep(log(theta + (1 - theta) * exp(-log_mu)), 4))
  
  got <- likelihood_combined(
    pred_combined = pred, params = params, z = z, y = y,
    x_vars = x_vars, component = 1,
    theta = theta, size = 1, N = N, dist = "ZIP"
  )
  expect_equal(got, expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(got, wrong)))
})

test_that("Component 2 NB likelihood uses ALL data.", {
  N <- 3
  z <- matrix(c(1,1,2, 2,2,3, 3,3,2), 3, 3, byrow = TRUE)   # three 2s
  y <- matrix(c(0,1,5, 3,2,4, 6,2,9), 3, 3, byrow = TRUE)
  params <- rep(0, 5)
  x_vars <- replicate(4, list(matrix(1, N, N)), simplify = FALSE)
  
  log_mu <- 0.7; size <- 4
  pred <- make_pred(log_mu)
  
  yc <- y[z == 2]
  expect_length(yc, 4)   # four 2s in z
  expected_all <- sum(dnbinom(yc, size = size, mu = exp(log_mu), log = TRUE))
  
  got <- likelihood_combined(
    pred_combined = pred, params = params, z = z, y = y,
    x_vars = x_vars, component = 2,
    theta = 0, size = size, N = N, dist = "NB"
  )
  expect_equal(got, expected_all, tolerance = 1e-12)
  
  ## The result must not match a single-observation likelihood, confirming
  ## the function sums over all of component 2's observations.
  single_obs <- dnbinom(yc[2], size = size, mu = exp(log_mu), log = TRUE)
  expect_false(isTRUE(all.equal(got, single_obs)))
})

test_that("Empty component returns 0 instead of erroring", {
  N <- 2
  z <- matrix(1, N, N)                   # no cells in component 2
  y <- matrix(1, N, N)
  params <- rep(0, 5)
  x_vars <- replicate(4, list(matrix(1, N, N)), simplify = FALSE)
  pred <- make_pred(0.5)
  
  expect_equal(
    likelihood_combined(
      pred_combined = pred, params = params, z = z, y = y,
      x_vars = x_vars, component = 2,
      theta = 0, size = 1, N = N, dist = "Poisson"
    ),
    0
  )
})

test_that("Invalid distribution raises", {
  N <- 2
  z <- matrix(1, N, N); y <- matrix(1, N, N)
  x_vars <- replicate(4, list(matrix(1, N, N)), simplify = FALSE)
  pred <- make_pred(0)
  expect_error(
    likelihood_combined(pred, rep(0, 5), z, y, x_vars, 1, 0, 1, N, "Nonsense"),
    "Invalid distribution"
  )
})

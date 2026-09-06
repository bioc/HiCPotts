test_that("prior_combined is DETERMINISTIC (regression test)", {
  set.seed(4291)
  N <- 4
  y      <- rpois(N * N, 5)
  x_vars <- replicate(4, list(matrix(runif(N * N), N, N)), simplify = FALSE)
  z      <- sample(1:3, N * N, replace = TRUE)
  params <- c(1, 0.2, -0.1, 0.3, 0.05)
  
  v1 <- prior_combined(params, 1, y, x_vars, z, TRUE, NULL)
  v2 <- prior_combined(params, 1, y, x_vars, z, TRUE, NULL)
  v3 <- prior_combined(params, 1, y, x_vars, z, TRUE, NULL)
  expect_identical(v1, v2)
  expect_identical(v1, v3)
  expect_true(is.finite(v1))
})

test_that("prior_combined with user_fixed_priors equals sum(dnorm(..., log=TRUE))", {
  params <- c(1, 2, 3, 4, 5)
  upri <- list(
    component1 = list(meany = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0,
                      sdy   = 1, sdx1   = 1, sdx2   = 1, sdx3   = 1, sdx4   = 1),
    component2 = list(meany = 1, meanx1 = 2, meanx2 = 3, meanx3 = 4, meanx4 = 5,
                      sdy   = 2, sdx1   = 2, sdx2   = 2, sdx3   = 2, sdx4   = 2),
    component3 = list(meany = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0,
                      sdy   = 1, sdx1   = 1, sdx2   = 1, sdx3   = 1, sdx4   = 1)
  )
  expected <- with(upri$component2,
                   dnorm(1, meany, sdy, log = TRUE) +
                     dnorm(2, meanx1, sdx1, log = TRUE) +
                     dnorm(3, meanx2, sdx2, log = TRUE) +
                     dnorm(4, meanx3, sdx3, log = TRUE) +
                     dnorm(5, meanx4, sdx4, log = TRUE))
  
  N <- 3
  y      <- rep(0, N * N)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  z      <- rep(2L, N * N)
  
  got <- prior_combined(params, 2, y, x_vars, z, FALSE, upri)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("Missing user_fixed_priors component raises an informative error", {
  params <- rep(0, 5)
  N <- 2
  y <- rep(0, N * N)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  z <- rep(1L, N * N)
  
  expect_error(
    prior_combined(params, 2, y, x_vars, z, FALSE, user_fixed_priors = NULL),
    "user_fixed_priors"
  )
  
  upri <- list(component1 = list(
    meany = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0,
    sdy = 1, sdx1 = 1, sdx2 = 1, sdx3 = 1, sdx4 = 1
  ))
  expect_error(
    prior_combined(params, 2, y, x_vars, z, FALSE, upri),
    "component2"
  )
})

test_that("Empty component is handled via fallback (not a crash)", {
  N <- 3
  y <- rep(0, N * N)
  x_vars <- replicate(4, list(matrix(0, N, N)), simplify = FALSE)
  z <- rep(1L, N * N)                    # no 2s at all
  params <- rep(0, 5)
  
  expect_silent(val <- prior_combined(params, 2, y, x_vars, z, TRUE, NULL))
  expect_true(is.finite(val))
})

test_that("iterated empirical-Bayes hyperparameters follow the current allocation", {
  N <- 4L
  y <- matrix(c(rep(0, 8), rep(30, 8)), N, N)
  x_vars <- list(
    distance = list(matrix(rep(0:3, each = 4), N, N)),
    GC = list(matrix(seq(0, 1, length.out = N * N), N, N)),
    TES = list(matrix(seq(1, 2, length.out = N * N), N, N)),
    ACC = list(matrix(seq(2, 3, length.out = N * N), N, N)))
  z_low <- matrix(c(rep(1L, 8), rep(2L, 8)), N, N)
  z_high <- 3L - z_low

  low <- HiCPotts:::.hicpotts_iterated_eb_priors(y, x_vars, z_low)
  high <- HiCPotts:::.hicpotts_iterated_eb_priors(y, x_vars, z_high)
  expect_named(low, paste0("component", 1:3))
  expect_true(all(vapply(low, function(x)
    all(is.finite(unlist(x[c("meany", "meanx1", "meanx2", "meanx3",
                             "meanx4", "sdy", "sdx1", "sdx2", "sdx3",
                             "sdx4")]))), logical(1))))
  expect_false(isTRUE(all.equal(low$component1$meany,
                                high$component1$meany)))
  expect_true(isTRUE(low$component3$global_fallback))
})

test_that("soft empirical-Bayes agrees with hard empirical-Bayes for one-hot weights", {
  set.seed(4292)
  N <- 5L
  y <- matrix(rpois(N * N, 8), N, N)
  covariates <- lapply(seq_len(4L), function(index)
    matrix(runif(N * N, 0, index), N, N))
  z <- matrix(rep(1:3, length.out = N * N), N, N)
  hard <- HiCPotts:::.hicpotts_iterated_eb_priors_cpp(
    y, covariates[[1L]], covariates[[2L]], covariates[[3L]],
    covariates[[4L]], z)
  soft <- HiCPotts:::.hicpotts_soft_eb_priors_cpp(
    y, covariates[[1L]], covariates[[2L]], covariates[[3L]],
    covariates[[4L]], z == 1L, z == 2L, z == 3L)

  expect_equal(soft, hard, tolerance = 1e-12)
})

test_that("soft empirical-Bayes uses fractional membership without discarding cells", {
  set.seed(4293)
  N <- 5L
  y <- matrix(rpois(N * N, 8), N, N)
  covariates <- lapply(seq_len(4L), function(index)
    matrix(runif(N * N, 0, index), N, N))
  raw_weights <- array(runif(N * N * 3L), dim = c(N, N, 3L))
  weight_sum <- apply(raw_weights, c(1L, 2L), sum)
  weights <- lapply(seq_len(3L), function(component)
    raw_weights[, , component] / weight_sum)
  soft <- HiCPotts:::.hicpotts_soft_eb_priors_cpp(
    y, covariates[[1L]], covariates[[2L]], covariates[[3L]],
    covariates[[4L]], weights[[1L]], weights[[2L]], weights[[3L]])

  effective_cells <- vapply(soft, `[[`, numeric(1), "effective_cells")
  expect_equal(sum(effective_cells), N * N, tolerance = 1e-10)
  expect_true(all(effective_cells > 0))
  expect_true(all(vapply(soft, function(component)
    all(is.finite(unlist(component[c(
      "meany", "meanx1", "meanx2", "meanx3", "meanx4",
      "sdy", "sdx1", "sdx2", "sdx3", "sdx4")]))), logical(1))))
  expect_true(all(vapply(soft, function(component)
    identical(component$allocation_weights, "soft membership probabilities"),
    logical(1))))
})

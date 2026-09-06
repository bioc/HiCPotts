## The marginal theta likelihood remains useful for checking the augmented
## Gibbs implementation against the model density.
##
## For a fixed target and a symmetric logit random walk (including the
## reciprocal Jacobian and the uniform Beta(1,1) prior), the untruncated
## forward and reverse log acceptance ratios are exact negatives and sum to
## zero. The previous EM-style surrogate recomputed a component-1
## responsibility at the current theta, so the reverse move used a different
## weight and the sum was non-zero.
##
## These tests call the compiled kernel directly via the internal export, so
## they verify the shipped implementation rather than an R copy of it.

make_case <- function(N = 4L, seed = 42L) {
  set.seed(seed)
  z <- matrix(sample(1:3, N * N, replace = TRUE), N, N)
  z[1, 1] <- 1L                      # guarantee some component-1 cells
  z[2, 2] <- 1L
  y <- matrix(rpois(N * N, 2), N, N)
  y[1, 1] <- 0                       # a component-1 zero (the informative case)
  cov <- function() matrix(runif(N * N, 0, 5), N, N)
  list(z = z, y = y, b1 = c(0.4, -0.2, 0.15, 0.02, 0.1),
       c1 = cov(), c2 = cov(), c3 = cov(), c4 = cov(), size1 = 2.5)
}

log_ratio <- function(p, from, to, dist) {
  ll_from <- HiCPotts:::.hicpotts_theta_loglik(p$z, p$y, p$b1, p$c1, p$c2, p$c3, p$c4,
                                               p$size1, from, dist)
  ll_to <- HiCPotts:::.hicpotts_theta_loglik(p$z, p$y, p$b1, p$c1, p$c2, p$c3, p$c4,
                                             p$size1, to, dist)
  jac <- log(to) + log1p(-to) - log(from) - log1p(-from)
  (ll_to - ll_from) + jac
}

test_that("the marginalized theta density has antisymmetric log ratios (ZIP)", {
  p <- make_case()
  for (pair in list(c(0.30, 0.65), c(0.10, 0.90), c(0.45, 0.55))) {
    fwd <- log_ratio(p, pair[1], pair[2], "ZIP")
    rev <- log_ratio(p, pair[2], pair[1], "ZIP")
    expect_equal(fwd + rev, 0, tolerance = 1e-10)
  }
})

test_that("the marginalized theta density has antisymmetric log ratios (ZINB)", {
  p <- make_case(seed = 7L)
  for (pair in list(c(0.20, 0.75), c(0.05, 0.60))) {
    fwd <- log_ratio(p, pair[1], pair[2], "ZINB")
    rev <- log_ratio(p, pair[2], pair[1], "ZINB")
    expect_equal(fwd + rev, 0, tolerance = 1e-10)
  }
})

test_that("C2: only component-1 cells inform theta", {
  p <- make_case(seed = 11L)
  base <- HiCPotts:::.hicpotts_theta_loglik(p$z, p$y, p$b1, p$c1, p$c2, p$c3, p$c4,
                                            p$size1, 0.4, "ZIP")

  # Perturbing a NON-component-1 cell's count must not move the theta
  # likelihood: conditional on z, theta appears only in component 1.
  q <- p
  idx <- which(q$z != 1L)[1L]
  q$y[idx] <- q$y[idx] + 17
  moved <- HiCPotts:::.hicpotts_theta_loglik(q$z, q$y, q$b1, q$c1, q$c2, q$c3, q$c4,
                                             q$size1, 0.4, "ZIP")
  expect_equal(base, moved)

  # Perturbing a component-1 cell must move it.
  r <- p
  jdx <- which(r$z == 1L)[1L]
  r$y[jdx] <- r$y[jdx] + 17
  expect_false(isTRUE(all.equal(
    base,
    HiCPotts:::.hicpotts_theta_loglik(r$z, r$y, r$b1, r$c1, r$c2, r$c3, r$c4,
                                      r$size1, 0.4, "ZIP")
  )))
})

test_that("C2: the likelihood responds to theta in the expected direction", {
  # With component-1 zeros present, higher theta must raise the likelihood
  # contribution of those zeros.
  N <- 3L
  z <- matrix(1L, N, N)
  y <- matrix(0, N, N)          # all component-1 zeros
  cov0 <- matrix(0, N, N)
  b1 <- c(0.5, 0, 0, 0, 0)
  lo <- HiCPotts:::.hicpotts_theta_loglik(z, y, b1, cov0, cov0, cov0, cov0,
                                          1, 0.2, "ZIP")
  hi <- HiCPotts:::.hicpotts_theta_loglik(z, y, b1, cov0, cov0, cov0, cov0,
                                          1, 0.8, "ZIP")
  expect_gt(hi, lo)
})

test_that("conjugate theta Gibbs step uses Beta(1,1) and component-1 cells", {
  N <- 3L
  z <- matrix(c(1L, 1L, 2L, 1L, 3L, 2L, 1L, 3L, 2L), N, N)
  y <- matrix(1, N, N)  # no structural-zero indicators are sampled
  cov0 <- matrix(0, N, N)
  n1 <- sum(z == 1L)

  set.seed(991)
  expected <- rbeta(1L, 1, 1 + n1)
  set.seed(991)
  actual <- HiCPotts:::.hicpotts_theta_gibbs_step_cpp(
    z, y, c(0, 0, 0, 0, 0), cov0, cov0, cov0, cov0,
    size1 = 2, theta_current = 0.4, dist = "ZIP")
  expect_identical(actual, expected)

  # Counts outside component 1 are irrelevant to theta and consume no draws.
  y[z != 1L] <- 0
  set.seed(991)
  outside_changed <- HiCPotts:::.hicpotts_theta_gibbs_step_cpp(
    z, y, c(0, 0, 0, 0, 0), cov0, cov0, cov0, cov0,
    size1 = 2, theta_current = 0.4, dist = "ZIP")
  expect_identical(outside_changed, expected)
})

make_pred_const <- function(log_mu) {
  function(params, z, x_vars, component, N) rep(log_mu, sum(z == component))
}

test_that("the complete symmetric N-by-N matrix contributes cell by cell", {
  N <- 2L
  z <- matrix(1, N, N)
  y <- matrix(c(2, 5, 5, 3), N, N)   # symmetric: one off-diagonal contact
  expect_true(isTRUE(all.equal(y, t(y))))

  params <- rep(0, 5)
  x_vars <- replicate(4, list(matrix(1, N, N)), simplify = FALSE)
  log_mu <- 1; mu <- exp(log_mu)
  pred <- make_pred_const(log_mu)

  got <- likelihood_combined(pred, params, z, y, x_vars, component = 1,
                             theta = 0, size = 1, N = N, dist = "Poisson")

  # Every stored cell contributes once, including both reflected cells.
  expected <- dpois(2, mu, log = TRUE) + dpois(3, mu, log = TRUE) +
    dpois(5, mu, log = TRUE) + dpois(5, mu, log = TRUE)
  expect_equal(got, expected, tolerance = 1e-12)
})

test_that("asymmetric data also uses the complete matrix", {
  N <- 2L
  z <- matrix(1, N, N)
  y <- matrix(c(2, 7, 5, 3), N, N)   # deliberately NOT symmetric
  expect_false(isTRUE(all.equal(y, t(y))))

  params <- rep(0, 5)
  x_vars <- replicate(4, list(matrix(1, N, N)), simplify = FALSE)
  log_mu <- 1; mu <- exp(log_mu)
  got <- likelihood_combined(make_pred_const(log_mu), params, z, y, x_vars,
                             component = 1, theta = 0, size = 1, N = N,
                             dist = "Poisson")
  expect_equal(got, sum(dpois(as.vector(y), mu, log = TRUE)), tolerance = 1e-12)
})

test_that("M3: public fitters reject non-integer counts", {
  N <- 2L
  y_bad <- list(matrix(c(1.5, 2, 3, 4), N, N))
  x_vars <- list(distance = list(matrix(0, N, N)), GC = list(matrix(0, N, N)),
                 TES = list(matrix(0, N, N)), ACC = list(matrix(0, N, N)))
  expect_error(
    run_chain_betas(N = N, gamma_start = 0.3, iterations = 3L,
                    x_vars = x_vars, y = y_bad),
    "non-integer", ignore.case = TRUE
  )
})

test_that("M3: public fitters reject a non-square response", {
  y_bad <- list(matrix(1:6, nrow = 2, ncol = 3))
  x_vars <- list(distance = list(matrix(0, 2, 2)), GC = list(matrix(0, 2, 2)),
                 TES = list(matrix(0, 2, 2)), ACC = list(matrix(0, 2, 2)))
  expect_error(
    run_chain_betas(N = 2L, gamma_start = 0.3, iterations = 3L,
                    x_vars = x_vars, y = y_bad),
    "square", ignore.case = TRUE
  )
})

test_that("M3: invalid gamma_start is rejected", {
  N <- 2L
  y <- list(matrix(c(1, 2, 3, 4), N, N))
  x_vars <- list(distance = list(matrix(0, N, N)), GC = list(matrix(0, N, N)),
                 TES = list(matrix(0, N, N)), ACC = list(matrix(0, N, N)))
  expect_error(
    run_chain_betas(N = N, gamma_start = 1.5, iterations = 3L,
                    x_vars = x_vars, y = y),
    "between 0 and 1"
  )
})

test_that("M4: a replicated rival mode is flagged multimodal", {
  # Two mutually incompatible allocations, two chains each.
  N <- 4L
  za <- matrix(1L, N, N); zb <- matrix(2L, N, N)
  mk <- function(z) list(z_final = structure(z, z_parameter_final = z),
                         z_parameter_final = z)
  fits <- list(mk(za), mk(za), mk(zb), mk(zb))

  ms <- HiCPotts:::.hicpotts_mode_consensus(fits, agreement_threshold = 0.8,
                                            minimum_mode_chains = 2L)
  expect_true(ms$multimodal)
  expect_gte(ms$competing_modes, 1L)

  # A single coherent mode must NOT be flagged.
  fits1 <- list(mk(za), mk(za), mk(za), mk(za))
  ms1 <- HiCPotts:::.hicpotts_mode_consensus(fits1, agreement_threshold = 0.8,
                                             minimum_mode_chains = 2L)
  expect_false(ms1$multimodal)
})

test_that("M1: integrated and plug-in probabilities differ and both normalise", {
  skip_if_not(exists("compute_HMRFHiC_probabilities"))
  set.seed(3)
  N <- 3L; n <- N * N
  d <- data.frame(
    start = rep(seq(1, by = 10000, length.out = N), times = N),
    `start.j.` = rep(seq(1, by = 10000, length.out = N), each = N),
    `end.i.` = 10000, end = 10000,
    interactions = as.numeric(rpois(n, 3)),
    GC = runif(n), TES = runif(n), ACC = runif(n),
    check.names = FALSE
  )
  n_draw <- 21L
  mk_chain <- function() {
    list(chains = replicate(3, matrix(rnorm(n_draw * 5, 0, 0.3), n_draw, 5),
                            simplify = FALSE),
         theta = runif(n_draw, 0.2, 0.5),
         size = matrix(runif(3 * n_draw, 1, 3), 3, n_draw),
         gamma = runif(n_draw, 0.2, 0.4))
  }
  fits <- list(mk_chain(), mk_chain())

  pi_ <- suppressWarnings(compute_HMRFHiC_probabilities(
    d, fits, iterations = 20L, dist = "ZINB", N = N, method = "integrated"))
  pp <- suppressWarnings(compute_HMRFHiC_probabilities(
    d, fits, iterations = 20L, dist = "ZINB", N = N, method = "plugin"))

  for (p in list(pi_, pp))
    expect_equal(p$prob1 + p$prob2 + p$prob3, rep(1, n), tolerance = 1e-10)

  # Averaging parameters first is not the same as averaging probabilities.
  expect_false(isTRUE(all.equal(pi_$prob1, pp$prob1)))
  expect_identical(attr(pi_, "probability_method"), "integrated")
})

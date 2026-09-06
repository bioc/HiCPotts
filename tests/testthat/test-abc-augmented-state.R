## C3: the gamma update must be a proper augmented-state ABC-MCMC kernel.
##
## State is the pair (gamma, x) with x the auxiliary dataset summarised by T.
## On rejection the WHOLE pair is retained, including the auxiliary draw. The
## previous implementation re-simulated the current gamma's auxiliary every
## iteration and discarded it on rejection, which is noisy MCMC rather than
## ABC-MCMC and does not target the stated ABC posterior.

abc_setup <- function(N = 6L, seed = 11L) {
  set.seed(seed)
  y <- matrix(rpois(N * N, 3), N, N)
  y[lower.tri(y)] <- t(y)[lower.tri(y)]
  list(N = N, y = y, x_vars = list(
    distance = list(abs(row(y) - col(y))),
    GC = list(matrix(0.5, N, N)),
    TES = list(matrix(1, N, N)),
    ACC = list(matrix(1, N, N))))
}

abc_run <- function(s, seed, iterations = 80L) {
  set.seed(seed)
  suppressWarnings(run_metropolis_MCMC_betas(
    s$N, 0.3, iterations, s$x_vars, s$y, TRUE, NULL, "Poisson",
    mcse_stop = FALSE, gamma_method = "abc"))
}

reference_potts_simulator <- function(N, gamma, sweeps) {
  z <- matrix(0, N, N)
  for (i in seq_len(N)) for (j in seq_len(N))
    z[i, j] <- floor(runif(1L, 1, 4))
  for (sweep in seq_len(sweeps)) {
    for (colour in 0:1) {
      for (i in seq_len(N)) {
        for (j in seq_len(N)) {
          if (((i - 1L + j - 1L) %% 2L) != colour) next
          neighbours <- c(
            if (i > 1L) z[i - 1L, j],
            if (i < N) z[i + 1L, j],
            if (j > 1L) z[i, j - 1L],
            if (j < N) z[i, j + 1L])
          counts <- tabulate(as.integer(neighbours), nbins = 3L)
          log_probability <- gamma * counts
          probability <- exp(log_probability - max(log_probability))
          draw <- runif(1L, 0, sum(probability))
          z[i, j] <- if (draw < probability[[1L]]) 1 else
            if (draw < probability[[1L]] + probability[[2L]]) 2 else 3
        }
      }
    }
  }
  z
}

test_that("compact ABC simulator exactly reproduces the reference transition", {
  for (seed in c(1L, 29L, 907L)) {
    set.seed(seed)
    reference <- reference_potts_simulator(5L, 0.37, 7L)
    set.seed(seed)
    compact <- HiCPotts:::.hicpotts_simulate_potts_labels_cpp(5L, 0.37, 7L)
    expect_identical(compact, reference)
  }
})

test_that("C3: the auxiliary draw is retained on rejection", {
  s <- abc_setup()
  fit <- abc_run(s, 1L)
  g <- as.numeric(fit$gamma)
  ## The conditioning value is S(z), the sufficient statistic of the sampler's
  ## own latent field, which moves every sweep. The retained
  ## auxiliary DRAW is therefore pinned by its summary T_sim, not by its
  ## distance: the draw is never re-simulated on rejection, but its distance to
  ## the (moving) S(z) is legitimately re-scored. T_sim is the invariant that
  ## actually encodes "the auxiliary was retained".
  d <- attr(fit$gamma, "T_sim")

  ## Index k is the state AFTER the transition into k. The first slot of the
  ## distance chain is never written (it has no preceding transition), so only
  ## consider transitions where both endpoints are finite.
  step <- seq.int(2L, length(g))
  usable <- step[is.finite(d[step]) & is.finite(d[step - 1L])]
  rejected <- usable[g[usable] == g[usable - 1L]]
  skip_if(length(rejected) < 5L, "too few rejections to test")

  # gamma did not move, so the auxiliary must not have moved either.
  expect_true(all(abs(d[rejected] - d[rejected - 1L]) < 1e-12))

  # Sanity: the test would catch a violation - accepted steps DO change it.
  accepted <- usable[g[usable] != g[usable - 1L]]
  skip_if(length(accepted) < 1L, "no acceptances observed")
  expect_true(any(abs(d[accepted] - d[accepted - 1L]) > 1e-12))
})

test_that("C3: T_sim and abc_distance describe the same auxiliary draw", {
  s <- abc_setup()
  fit <- abc_run(s, 2L)
  t_obs <- attr(fit$gamma, "T_obs")
  t_sim <- attr(fit$gamma, "T_sim")
  d <- attr(fit$gamma, "abc_distance")
  expect_equal(abs(t_sim - t_obs), d, tolerance = 1e-10)
})

test_that("C3: the ABC tolerance does not depend on the chain seed", {
  s <- abc_setup()
  eps <- vapply(c(1L, 2L, 3L), function(sd)
    attr(abc_run(s, sd, iterations = 40L)$gamma, "abc_epsilon"), numeric(1))
  expect_equal(length(unique(signif(eps, 10))), 1L)
  expect_true(is.finite(eps[1]) && eps[1] > 0)
})

test_that("C3: fixing the calibration stream leaves the sampler random", {
  # The calibration substream must not leak into the chain: different seeds
  # must still give different chains, and the RNG stream must be restored.
  s <- abc_setup()
  f1 <- abc_run(s, 1L, 40L)
  f2 <- abc_run(s, 2L, 40L)
  expect_false(isTRUE(all.equal(as.numeric(f1$gamma), as.numeric(f2$gamma))))
  expect_false(isTRUE(all.equal(f1$chains[[1]][, 1], f2$chains[[1]][, 1])))

  # Same seed still reproduces exactly.
  expect_equal(as.numeric(abc_run(s, 7L, 40L)$gamma),
               as.numeric(abc_run(s, 7L, 40L)$gamma))

  # And the caller's stream continues predictably after the run.
  set.seed(5); invisible(abc_run(s, 5L, 40L)); a <- runif(3)
  set.seed(5); invisible(abc_run(s, 5L, 40L)); b <- runif(3)
  expect_equal(a, b)
})

test_that("C3: a user-supplied epsilon is honoured exactly", {
  s <- abc_setup()
  set.seed(3)
  fit <- suppressWarnings(run_metropolis_MCMC_betas(
    s$N, 0.3, 40L, s$x_vars, s$y, TRUE, NULL, "Poisson",
    epsilon = 0.25, mcse_stop = FALSE, gamma_method = "abc"))
  expect_equal(attr(fit$gamma, "abc_epsilon"), 0.25)
  expect_identical(attr(fit$gamma, "abc_calibration_reps_executed"), 0L)
  expect_identical(attr(fit$gamma, "abc_simulator_engine"),
                   "compact_lookup_v1")
})

test_that("repeated chains calibrate epsilon once outside sampler workers", {
  s <- abc_setup(N = 5L)
  fits <- suppressWarnings(run_hicpotts_chains(
    N = s$N, gamma_prior = 0.3, iterations = 8L,
    x_vars = s$x_vars, y = s$y, dist = "Poisson",
    seeds = 31:32, initialization = "random", relabel = FALSE,
    mcse_stop = FALSE, abc_potts_sweeps = 3L, abc_sim_reps = 1L))
  epsilon <- vapply(fits, function(fit)
    attr(fit$gamma, "abc_epsilon"), numeric(1))
  expect_equal(length(unique(signif(epsilon, 12L))), 1L)
  expect_true(all(vapply(fits, function(fit)
    isTRUE(attr(fit$gamma, "abc_epsilon_shared_calibration")), logical(1))))
  expect_true(all(vapply(fits, function(fit)
    identical(attr(fit$gamma, "abc_calibration_reps_executed"), 0L),
    logical(1))))
  expect_true(all(vapply(fits, function(fit)
    identical(attr(fit$gamma, "abc_shared_calibration_reps"), 60L),
    logical(1))))
})

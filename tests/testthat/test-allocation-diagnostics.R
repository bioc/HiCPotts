## Item 6: batched latent-state frequencies give per-cell MCSE, allocation ESS
## and between-chain disagreement, without changing the three-way MAP label.

make_run <- function(seed, N = 8L, iterations = 400L) {
  set.seed(seed)
  y <- matrix(rpois(N * N, 4), N, N)
  x_vars <- list(
    distance = list(abs(row(y) - col(y))),
    GC = list(matrix(runif(N * N), N, N)),
    TES = list(matrix(runif(N * N), N, N)),
    ACC = list(matrix(runif(N * N), N, N)))
  suppressWarnings(run_metropolis_MCMC_betas(
    N = N, gamma_prior = 0.4, iterations = iterations,
    x_vars = x_vars, y = y, use_data_priors = TRUE, dist = "Poisson",
    gamma_update_interval = 10L))
}

test_that("the sampler stores batched latent-state frequencies", {
  fit <- make_run(101)
  b <- fit$z_probability_batches
  expect_true(is.list(b))
  d <- dim(b$batch_means)
  expect_length(d, 4L)
  expect_identical(d[3L], 3L)
  expect_gte(d[4L], 2L)
  expect_identical(length(b$batch_draws), d[4L])
  expect_true(all(b$batch_draws > 0))
  ## Batch draw counts must account for every retained draw.
  expect_identical(sum(b$batch_draws), fit$z_probability_draws)
  ## Each batch mean is a probability.
  expect_true(all(b$batch_means >= 0 & b$batch_means <= 1))
})

test_that("batch means reproduce the pooled probabilities", {
  fit <- make_run(102)
  b <- fit$z_probability_batches
  w <- as.numeric(b$batch_draws) / sum(b$batch_draws)
  d <- dim(b$batch_means)
  for (k in 1:3) {
    recombined <- matrix(0, d[1L], d[2L])
    for (i in seq_len(d[4L]))
      recombined <- recombined + b$batch_means[, , k, i] * w[i]
    expect_equal(recombined, fit$z_probabilities[[k]], tolerance = 1e-10)
  }
})

test_that("allocation_diagnostics reports per-cell MCSE and ESS", {
  fit <- make_run(103)
  dg <- allocation_diagnostics(fit)
  expect_named(dg$mcse, c("noise", "signal", "false signal"))
  expect_true(all(vapply(dg$mcse, function(m) all(m >= 0), logical(1))))
  expect_true(all(vapply(dg$ess, function(m) all(m >= 0), logical(1))))
  ## ESS can never exceed the number of draws actually taken.
  expect_lte(dg$worst_cell_ess, fit$z_probability_draws)
  expect_true(is.finite(dg$max_mcse))
  expect_null(dg$between_chain_disagreement)
})

test_that("between-chain disagreement is reported for multiple chains", {
  fits <- list(make_run(201), make_run(202))
  dg <- allocation_diagnostics(fits)
  d <- dg$between_chain_disagreement
  expect_false(is.null(d))
  expect_identical(d$n_chain_pairs, 1L)
  expect_true(d$mean >= 0 && d$mean <= 1)
  expect_true(d$max >= d$mean)
})

test_that("a component 2/3 exchange is counted as biological disagreement", {
  fit <- make_run(204)
  flipped <- fit
  flipped$z_probabilities <- fit$z_probabilities[c(1L, 3L, 2L)]
  dg <- allocation_diagnostics(list(fit, flipped))
  expect_gt(dg$between_chain_disagreement$mean, 1e-4)
})

test_that("diagnostics do not alter the three-way MAP classification", {
  fit <- make_run(105)
  before <- classify_hicpotts(fit, relabel = FALSE, min_draws = 10L)
  invisible(allocation_diagnostics(fit))
  after <- classify_hicpotts(fit, relabel = FALSE, min_draws = 10L)
  expect_identical(before$map_component, after$map_component)
  expect_identical(before$classification, after$classification)
})

test_that("a near-degenerate run is refused rather than classified", {
  ## z_probability_burnin_arg is validated to be < iterations, so at least one
  ## draw always accumulates and the C++ zero-draw guard is a defensive
  ## backstop rather than a reachable path. The user-visible guarantee is that
  ## a run this short cannot be turned into a classification: there is no
  ## substitution of the single final latent state, and the draw floor rejects
  ## it outright.
  set.seed(7); N <- 4L
  y <- matrix(rpois(N * N, 3), N, N)
  x_vars <- list(distance = list(abs(row(y) - col(y))),
                 GC = list(matrix(runif(N * N), N, N)),
                 TES = list(matrix(runif(N * N), N, N)),
                 ACC = list(matrix(runif(N * N), N, N)))
  fit <- suppressWarnings(run_metropolis_MCMC_betas(
    N = N, gamma_prior = 0.4, iterations = 6L, x_vars = x_vars, y = y,
    use_data_priors = TRUE, dist = "Poisson",
    z_probability_burnin_arg = 5L))
  expect_identical(fit$z_probability_draws, 1L)
  expect_error(classify_hicpotts(fit, relabel = FALSE),
               "below the required minimum")
})

test_that("whole 2/3 branch prior ratio is an exact involution", {
  beta1 <- c(0, -0.1, 0.2, -0.2, 0.1)
  beta2 <- c(3, -0.1, 0.2, -0.2, 0.1)
  beta3 <- c(3, -1.2, 1.4, 0.9, -1.0)
  sds <- c(0.7, 0.2, 0.2, 0.2)

  forward <- HiCPotts:::.hicpotts_branch_prior_logratio_cpp(
    beta1, beta2, beta3, sds)
  reverse <- HiCPotts:::.hicpotts_branch_prior_logratio_cpp(
    beta1, beta3, beta2, sds)

  expect_gt(forward, 0)
  expect_equal(forward + reverse, 0, tolerance = 1e-12)
})

test_that("gamma mixture proposal ratio is antisymmetric", {
  forward <- HiCPotts:::.hicpotts_gamma_proposal_logratio_cpp(
    current = 0.17, proposed = 0.83, local_step = 0.31,
    large_probability = 0.2, large_multiplier = 3,
    independence_probability = 0.1,
    prior_shape1 = 1.5, prior_shape2 = 2.5)
  reverse <- HiCPotts:::.hicpotts_gamma_proposal_logratio_cpp(
    current = 0.83, proposed = 0.17, local_step = 0.31,
    large_probability = 0.2, large_multiplier = 3,
    independence_probability = 0.1,
    prior_shape1 = 1.5, prior_shape2 = 2.5)

  expect_true(is.finite(forward))
  expect_equal(forward + reverse, 0, tolerance = 1e-12)
})

test_that("known-gamma Potts truth is reproducible and records its generator", {
  a <- simulate_hicpotts_potts_truth(
    N = 4, gamma = 0.4, potts_sweeps = 8, seed = 91)
  b <- simulate_hicpotts_potts_truth(
    N = 4, gamma = 0.4, potts_sweeps = 8, seed = 91)

  expect_identical(a$z_true, b$z_true)
  expect_identical(a$y, b$y)
  expect_equal(a$settings$gamma, 0.4)
  expect_match(a$settings$spatial_model, "Potts")
  expect_true(all(a$z_true %in% 1:3))
  expect_equal(sum(a$settings$proportions), 1)
})

test_that("sampler reports branch, block and gamma-mixture movement", {
  N <- 3L
  y <- matrix(c(0, 1, 2, 1, 8, 10, 2, 11, 18), N, N)
  distance <- abs(row(y) - col(y))
  mk <- function(offset) list(matrix(seq(0.05, 0.85, length.out = N * N) +
                                      offset, N, N))
  x_vars <- list(distance = list(distance), GC = mk(0),
                 TES = mk(0.03), ACC = mk(0.06))
  z_start <- matrix(c(1, 1, 2, 1, 2, 2, 3, 3, 3), N, N)

  set.seed(2026)
  fit <- suppressWarnings(run_metropolis_MCMC_betas(
    N = N, gamma_prior = 0.3, iterations = 90,
    x_vars = x_vars, y = y, use_data_priors = TRUE,
    dist = "ZIP", theta_start = 0.3, z_start = z_start,
    mcse_stop = FALSE, z_probability_burnin_arg = 45,
    branch_swap_interval = 5, signal_block_move_interval = 5,
    gamma_update_interval = 15, abc_potts_sweeps_arg = 2,
    abc_sim_reps = 1))

  expect_true(is.list(fit$branch_mixing))
  expect_true(isTRUE(fit$branch_mixing$enabled))
  expect_gt(fit$branch_mixing$attempts, 0)
  expect_length(fit$branch_mixing$branch_gap, 91)
  expect_true(is.list(fit$signal_block_mixing))
  expect_gt(fit$signal_block_mixing$attempts, 0)

  proposal_counts <- c(
    attr(fit$gamma, "gamma_local_proposals"),
    attr(fit$gamma, "gamma_large_proposals"),
    attr(fit$gamma, "gamma_independence_proposals"))
  expect_equal(sum(proposal_counts), attr(fit$gamma, "gamma_proposals"))
  expect_true(is.finite(attr(fit$gamma, "gamma_acceptance_rate")))
  expect_identical(attr(fit$gamma, "gamma_method"), "abc")
  expect_identical(attr(fit$gamma, "gamma_update"),
                   "ABC_potts_label_agreement")
  expect_identical(attr(fit$theta, "theta_update"),
                   "Beta(1,1)-Bernoulli conjugate Gibbs")
})

test_that("the removed exchange gamma method is rejected", {
  N <- 2L
  x_vars <- list(
    distance = list(matrix(0, N, N)), GC = list(matrix(0, N, N)),
    TES = list(matrix(0, N, N)), ACC = list(matrix(0, N, N)))
  expect_error(
    run_metropolis_MCMC_betas(
      N, 0.3, 2, x_vars, matrix(1, N, N), TRUE,
      dist = "Poisson", mcse_stop = FALSE,
      z_probability_burnin_arg = 0, abc_potts_sweeps_arg = 1,
      abc_sim_reps = 1, gamma_method = "exchange"),
    "removed")
})

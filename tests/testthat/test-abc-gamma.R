test_that("gamma uses join-count-excess ABC and preserves result structure", {
  set.seed(134)
  N <- 4L
  y <- matrix(rpois(N * N, 3), N, N)
  x_vars <- list(
    distance = list(abs(row(y) - col(y))),
    GC = list(matrix(runif(N * N), N, N)),
    TES = list(matrix(runif(N * N), N, N)),
    ACC = list(matrix(runif(N * N), N, N))
  )
  fit <- run_chain_betas(
    N = N, gamma_start = 0.45, iterations = 20L,
    x_vars = x_vars, y = list(y), use_data_priors = TRUE,
    dist = "Poisson", epsilon = NULL, mc_cores = 1L
  )[[1]]

  expect_true(all(c("chains", "gamma", "theta", "size",
                      "z_final", "z_checkpoints", "z_probabilities",
                      "z_probability_draws", "z_probability_burnin",
                      "branch_mixing", "signal_block_mixing",
                      "dispersion_mixing", "comp23_barrier",
                      "noise_relationship", "regression_prior",
                      "dispersion_prior", "pair_weighting",
                      "z_probability_batches", "performance",
                      "sampler_settings", "data_fingerprint",
                      "component_definition", "provenance") %in% names(fit)))
  expect_identical(attr(fit$gamma, "gamma_update"), "ABC_potts_label_agreement")
  expect_true(all(is.finite(fit$gamma)))
  expect_true(all(fit$gamma > 0 & fit$gamma < 1))
  expect_true(all(is.finite(attr(fit$gamma, "abc_distance")[-1L])))
  expect_true(all(is.finite(attr(fit$gamma, "T_obs")[-1L])))
  expect_true(all(is.finite(attr(fit$gamma, "T_sim")[-1L])))
  expect_gt(attr(fit$gamma, "abc_epsilon"), 0)
  expect_gte(attr(fit$gamma, "gamma_acceptance_rate"), 0)
  expect_lte(attr(fit$gamma, "gamma_acceptance_rate"), 1)
})

test_that("invalid ABC bandwidth is rejected", {
  N <- 3L
  y <- matrix(1, N, N)
  x_vars <- list(
    distance = list(matrix(0, N, N)),
    GC = list(matrix(0, N, N)),
    TES = list(matrix(0, N, N)),
    ACC = list(matrix(0, N, N))
  )
  expect_error(
    run_chain_betas(N = N, gamma_start = 0.4, iterations = 2L,
                    x_vars = x_vars, y = list(y), use_data_priors = TRUE,
                    dist = "Poisson", epsilon = 0),
    "epsilon must be NULL or a finite positive"
  )
})

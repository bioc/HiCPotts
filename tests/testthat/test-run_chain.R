test_that("Per-dataset x_vars subsetting is wired correctly (regression test)", {
  set.seed(1L)
  N <- 5
  y_list <- list(
    matrix(rpois(N * N, 2), N, N),
    matrix(rpois(N * N, 8), N, N)
  )
  ## DIFFERENT covariate matrices for the two datasets
  x_vars <- list(
    distance = list(matrix(runif(N * N, 0, 1), N, N),
                    matrix(runif(N * N, 0, 1), N, N)),
    GC       = list(matrix(runif(N * N, 0, 1), N, N),
                    matrix(runif(N * N, 0, 1), N, N)),
    TES      = list(matrix(runif(N * N, 0, 1), N, N),
                    matrix(runif(N * N, 0, 1), N, N)),
    ACC      = list(matrix(runif(N * N, 0, 1), N, N),
                    matrix(runif(N * N, 0, 1), N, N))
  )
  res <- run_chain_betas(
    N = N, iterations = 3, x_vars = x_vars, y = y_list,
    theta_start = 0.5, size_start = c(1, 1, 1),
    use_data_priors = TRUE, dist = "Poisson",
    mc_cores = 1L
  )
  expect_length(res, 2L)
  expected_names <- c("chains", "gamma", "theta", "size", "z_final",
    "z_checkpoints", "z_probabilities", "z_probability_draws",
    "z_probability_burnin", "beta_mixing", "branch_mixing",
    "signal_block_mixing",
    "dispersion_mixing",
    "comp23_barrier", "noise_relationship", "regression_prior",
    "dispersion_prior", "pair_weighting",
    "z_probability_batches", "performance", "sampler_settings",
    "data_fingerprint", "component_definition", "provenance")
  expect_true(all(expected_names %in% names(res[[1]])))
  expect_true(all(expected_names %in% names(res[[2]])))
})

test_that("mc_cores > 1 is accepted on Windows without forced downgrade", {
  skip_on_os(c("mac", "linux", "solaris"))
  N <- 3
  y <- list(matrix(rpois(N * N, 2), N, N))
  x_vars <- list(
    distance = list(matrix(0, N, N)),
    GC       = list(matrix(0, N, N)),
    TES      = list(matrix(0, N, N)),
    ACC      = list(matrix(0, N, N))
  )
  result <- suppressWarnings(
    run_chain_betas(N = N, iterations = 2, x_vars = x_vars, y = y,
                    use_data_priors = TRUE, dist = "Poisson", mc_cores = 2L))
  expect_length(result, 1L)
})

test_that("Bad x_vars/ y structure errors early", {
  expect_error(run_chain_betas(N = 3, iterations = 2,
                               x_vars = list(), y = matrix(0, 3, 3),
                               use_data_priors = TRUE),
               "named list")
  expect_error(run_chain_betas(N = 3, iterations = 2,
                               x_vars = list(),
                               y = list(matrix(0, 3, 3)),
                               use_data_priors = TRUE),
               "named list")
})

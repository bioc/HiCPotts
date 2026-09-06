test_that("run_chain_betas controls the number of chains for one matrix", {
  set.seed(77)
  N <- 3L
  y <- matrix(rpois(N * N, 3), N, N)
  mk <- function() matrix(runif(N * N), N, N)
  x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())

  fits <- suppressWarnings(run_chain_betas(
    N = N, x_vars = x_vars, y = y, dist = "Poisson",
    iterations = 3L, seeds = c(101L, 202L), robust = FALSE,
    epsilon = 1))

  expect_length(fits, 2L)
  expect_named(fits, c("chain1", "chain2"))
  expect_identical(unname(vapply(fits, `[[`, integer(1), "seed")),
                   c(101L, 202L))
  expect_identical(
    unname(vapply(fits, `[[`, character(1), "initialization_method")),
    c("likelihood_informed", "likelihood_informed"))
})

test_that("run_chain_betas exposes and forwards public sampler controls", {
  public_controls <- c(
    "gamma_start", "initialization", "gamma_prior_shape", "verbose",
    "progress_interval", "abc_epsilon_quantile", "gamma_update_interval",
    "z_probability_burnin", "mcse_stop", "mcse_min_iterations",
    "mcse_check_interval", "mcse_relative_threshold", "diagnostic_control"
  )
  expect_true(all(public_controls %in% names(formals(run_chain_betas))))

  set.seed(79)
  N <- 3L
  y <- matrix(rpois(N * N, 3), N, N)
  mk <- function() matrix(runif(N * N), N, N)
  x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())

  output <- capture.output(fits <- run_chain_betas(
    N = N, x_vars = x_vars, y = y, dist = "Poisson",
    iterations = 2L, burnin = 0L, n_chains = 1L,
    initialization = "random", gamma_prior_shape = c(2, 3),
    epsilon = 1, gamma_update_interval = 1L,
    z_probability_burnin = 0L, mcse_stop = FALSE,
    verbose = TRUE, progress_interval = 1L
  ))
  fit <- fits[[1L]]

  expect_true(any(grepl("Iteration: 1", output, fixed = TRUE)))
  expect_identical(fit$initialization_method, "random")
  expect_equal(attr(fit$gamma, "gamma_prior_shape1"), 2)
  expect_equal(attr(fit$gamma, "gamma_prior_shape2"), 3)
  expect_identical(fit$sampler_settings$gamma_update_interval, 1L)
  expect_identical(fit$z_probability_burnin, 0L)
  expect_false(attr(fit$gamma, "mcse_stopping_enabled"))
  expect_true(fit$sampler_settings$verbose)
  expect_identical(fit$sampler_settings$progress_interval, 1L)
})

test_that("deprecated gamma_prior remains a checked compatibility alias", {
  N <- 3L
  y <- matrix(1, N, N)
  zero <- matrix(0, N, N)
  x_vars <- list(distance = zero, GC = zero, TES = zero, ACC = zero)

  expect_warning(
    fit <- run_chain_betas(
      N = N, x_vars = x_vars, y = y, dist = "Poisson",
      gamma_prior = 0.4, iterations = 1L, burnin = 0L,
      initialization = "random", epsilon = 1,
      z_probability_burnin = 0L, mcse_stop = FALSE
    ),
    "deprecated"
  )
  expect_equal(fit[[1L]]$gamma[[1L]], 0.4)
})

test_that("run_chain_betas exposes the robust workflow without exporting helpers", {
  set.seed(78)
  N <- 3L
  y <- matrix(rpois(N * N, 3), N, N)
  mk <- function() matrix(runif(N * N), N, N)
  x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())

  fit <- suppressWarnings(run_chain_betas(
    N = N, x_vars = x_vars, y = y, dist = "Poisson",
    iterations = 3L, burnin = 1L, n_chains = 1L, robust = TRUE,
    epsilon = 1,
    diagnostic_control = list(
      minimum_component_cells = 0L, minimum_ess = 1,
      maximum_rhat = 1.2, gamma_boundary_tolerance = 0.02,
      maximum_gamma_boundary_fraction = 0.9,
      minimum_gamma_unique = 2L, mode_agreement_threshold = 0.7,
      minimum_mode_chains = 1L
    )))

  expect_s3_class(fit, "hicpotts_robust_fit")
  expect_length(fit$all_fits, 1L)
  expect_identical(fit$settings$seeds, 1L)
  expect_identical(fit$settings$minimum_component_cells, 0L)
  expect_equal(fit$settings$minimum_ess, 1)
  expect_equal(fit$settings$maximum_rhat, 1.2)
  expect_equal(fit$settings$mode_agreement_threshold, 0.7)
  expect_identical(fit$settings$minimum_mode_chains, 1L)
  expect_true(fit$settings$shared_empirical_bayes_prior)
  expect_true(attr(fit$priors, "shared_across_production_chains"))
  expect_identical(
    fit$all_fits[[1L]]$regression_prior$method,
    "shared_soft_allocation_empirical_Bayes_pilot_frozen"
  )
  expect_true(fit$all_fits[[1L]]$regression_prior$shared_across_chains)
    expect_true(is.list(fit$empirical_bayes_pilot))
})

test_that("pilot progress distinguishes pilot and requested iterations", {
    expect_identical(
        .hicpotts_pilot_progress_message(4L, 5000L, 12000L),
        paste0(
            "Running 4 shared-prior pilot chains for 5000 pilot iterations; ",
            "production will then run up to 12000 requested iterations per ",
            "chain."
        )
    )
})

test_that("zero is the default and accepted component occupancy threshold", {
  expect_identical(
    formals(diagnose_hicpotts_fit)$minimum_component_cells,
    0L
  )
})

test_that("diagnostic_control rejects unknown settings before fitting", {
  expect_error(
    run_chain_betas(
      N = 2L, x_vars = list(), y = matrix(0, 2L, 2L),
      diagnostic_control = list(unknown_setting = 1)
    ),
    "Unknown 'diagnostic_control'"
  )
})

test_that("parameter summaries can diagnose covariates internally", {
  diagnostics <- list(
    parameters = data.frame(parameter = "gamma", estimate = 0.3),
    reliability_flags = data.frame(
      criterion = "covariate_conditioning", passed = TRUE,
      threshold = "correlation/condition thresholds"),
    reliable = TRUE, warnings = character())
  N <- 3L
  constant <- matrix(0, N, N)
  x_vars <- list(distance = constant, GC = constant,
                 TES = constant, ACC = constant)

  summary <- summarise_hicpotts_parameters(
    diagnostics, require_reliable = FALSE, x_vars = x_vars)

  # The overall verdict attribute was removed; the failed gate is named
  # explicitly instead.
  expect_true("covariate_conditioning" %in%
              attr(summary, "failed_reliability_gates"))
  expect_type(attr(summary, "covariate_diagnostics"), "list")
  expect_error(
    summarise_hicpotts_parameters(diagnostics, x_vars = x_vars),
    "covariate_conditioning")
})

test_that("component definitions are attached to computed probabilities", {
  n_draw <- 6L
  make_chain <- function(intercept) cbind(
    rep(intercept, n_draw), matrix(0, n_draw, 4L))
  fit <- list(
    chains = list(make_chain(-1), make_chain(1), make_chain(2)),
    theta = rep(0.2, n_draw),
    size = matrix(10, 3L, n_draw),
    gamma = rep(0.3, n_draw))
  data <- data.frame(
    start = 1:4, end = 2:5, interactions = c(0, 2, 4, 6),
    GC = c(0.1, 0.2, 0.3, 0.4), TES = c(0.2, 0.3, 0.4, 0.5),
    ACC = c(0.3, 0.4, 0.5, 0.6))
  definition <- data.frame(
    component = 1:3, label = c("background", "interaction", "artefact"),
    definition = c("low", "biological", "elevated background"))

  probabilities <- suppressWarnings(compute_HMRFHiC_probabilities(
    data, fit, iterations = n_draw - 1L, dist = "ZINB",
    method = "plugin", component_definition = definition))

  expect_identical(attr(probabilities, "component_definition"), definition)
  expect_identical(
    unname(attr(probabilities, "probability_components")), definition$label)
})

test_that("only the intended public interface is exported", {
  public <- getNamespaceExports("HiCPotts")
  expect_true(all(c("run_chain_betas", "summarise_hicpotts_parameters",
                    "compute_HMRFHiC_probabilities") %in% public))
  expect_false(any(c("fit_hicpotts_robust", "run_hicpotts_chains",
                     "diagnose_hicpotts_covariates",
                     "hicpotts_component_definition",
                     "run_metropolis_MCMC_betas") %in% public))
})

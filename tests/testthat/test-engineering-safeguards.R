test_that("cached native component posteriors equal the public reference", {
  set.seed(315)
  N <- 4L
  y <- matrix(rpois(N * N, 3), N, N)
  z <- matrix(rep(1:3, length.out = N * N), N, N)
  covariates <- lapply(c(8, 1, 2, 3), function(scale)
    matrix(runif(N * N, 0, scale), N, N))
  x_vars <- lapply(covariates, list)
  priors <- setNames(lapply(1:3, function(component) list(
    meany = component / 10, meanx1 = -0.1, meanx2 = 0.2,
    meanx3 = 0.05, meanx4 = -0.05,
    sdy = 1.2, sdx1 = 0.7, sdx2 = 1.4, sdx3 = 0.9, sdx4 = 1.1)),
    paste0("component", 1:3))
  betas <- list(
    c(0.2, -0.1, 0.1, 0.05, -0.02),
    c(0.9, -0.2, 0.3, -0.1, 0.15),
    c(0.6, -0.08, 0.12, 0.02, -0.03))

  for (dist in c("Poisson", "ZIP", "NB", "ZINB")) {
    for (component in 1:3) {
      size_value <- 2.5 + component
      reference <- posterior_combined(
        pred_combined, betas[[component]], z, y, x_vars,
        component = component, theta = 0.23, N = N,
        use_data_priors = FALSE, user_fixed_priors = priors,
        dist = dist,
        size = if (dist %in% c("NB", "ZINB")) size_value else NULL)
      cached <- HiCPotts:::.hicpotts_cached_posterior_cpp(
        betas[[component]], z, y, component, 0.23, size_value, dist,
        covariates[[1]], covariates[[2]], covariates[[3]],
        covariates[[4]], priors)
      expect_equal(cached, reference, tolerance = 1e-11,
        info = paste(dist, "component", component))
    }
  }
})

test_that("direct native allocation weights equal the reference conditional", {
  set.seed(316)
  N <- 4L
  y <- matrix(rpois(N * N, 5), N, N)
  z <- matrix(sample(1:3, N * N, replace = TRUE), N, N)
  covariates <- lapply(c(10, 1, 2, 3), function(scale)
    matrix(runif(N * N, 0, scale), N, N))
  x_vars <- setNames(lapply(covariates, list),
                     c("distance", "GC", "TES", "ACC"))
  betas <- rbind(
    c(0.1, -0.1, 0.2, 0.05, -0.02),
    c(1.0, -0.3, 0.4, -0.1, 0.2),
    c(0.7, -0.08, 0.18, 0.03, -0.01))
  sizes <- c(3, 7, 12)
  theta <- 0.27
  gamma <- 0.34
  row <- 2L
  column <- 3L

  for (dist in c("Poisson", "ZIP", "NB", "ZINB")) {
    native <- HiCPotts:::.hicpotts_cell_logweights_cpp(
      z, y, betas, sizes, theta, gamma, dist,
      covariates[[1]], covariates[[2]], covariates[[3]],
      covariates[[4]], row, column)
    reference <- vapply(1:3, function(component) {
      candidate <- z
      candidate[row, column] <- component
      neighbours <- Neighbours_combined(candidate, N)
      chains <- lapply(1:3, function(k)
        matrix(betas[k, ], nrow = 1L))
      pz_123(candidate, neighbours, y, pred_combined, chains, gamma,
             x_vars, theta, matrix(sizes, nrow = 3L), N, 0L,
             dist)[row, column]
    }, numeric(1))
    native_probability <- exp(native - max(native))
    native_probability <- native_probability / sum(native_probability)
    reference_probability <- exp(reference - max(reference))
    reference_probability <- reference_probability / sum(reference_probability)
    expect_equal(unname(native_probability), reference_probability,
                 tolerance = 1e-11, info = dist)
  }
})

test_that("public repeated-chain fits are reproducible and self-describing", {
  N <- 4L
  set.seed(317)
  y <- matrix(rpois(N * N, 3), N, N)
  mk <- function() list(matrix(runif(N * N), N, N))
  x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
  arguments <- list(
    N = N, gamma_prior = 0.3, iterations = 6L, x_vars = x_vars,
    y = y, dist = "Poisson", seeds = 91L, initialization = "random",
    mcse_stop = FALSE, gamma_update_interval = 2L, abc_sim_reps = 1L,
    abc_potts_sweeps = 1L, relabel = FALSE)
  first <- suppressWarnings(do.call(run_hicpotts_chains, arguments)[[1L]])
  second <- suppressWarnings(do.call(run_hicpotts_chains, arguments)[[1L]])

  expect_equal(first$gamma, second$gamma)
  expect_equal(first$chains, second$chains)
  expect_identical(first$seed, 91L)
  expect_identical(first$initialization_method, "random")
  expect_identical(first$component_definition$label,
                   c("noise", "signal", "false signal"))
  expect_s3_class(first$data_fingerprint, "hicpotts_fingerprint")
  expect_true(is.list(first$performance))
  expect_true(is.finite(first$performance$seconds_per_iteration))
  expect_true(is.list(first$sampler_settings))
})

test_that("validation rejects inputs that cannot define the fitted likelihood", {
  N <- 3L
  y <- matrix(1, N, N)
  mk <- function(value = 0) list(matrix(value, N, N))
  x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
  x_bad <- x_vars
  x_bad$GC[[1]][1, 1] <- -1
  expect_error(run_hicpotts_chains(
    N, iterations = 2L, x_vars = x_bad, y = y, seeds = 1L),
    "greater than -1|<= -1")
  expect_error(run_hicpotts_chains(
    N, iterations = 2L, x_vars = x_vars, y = y + 0.1, seeds = 1L),
    "non-integer|RAW integer")
  expect_error(run_hicpotts_chains(
    N, iterations = 2L, x_vars = x_vars, y = y, dist = "NB",
    seeds = 1L), "size_start")
})

test_that("component identities have one canonical package definition", {
  definition <- hicpotts_component_definition()
  expect_identical(definition$component, 1:3)
  expect_identical(definition$label, c("noise", "signal", "false signal"))
  expect_match(definition$definition[[3L]], "component 1")
})

test_that("oracle-validation controls hold every requested state block fixed", {
  truth <- simulate_hicpotts_potts_truth(
    N = 5L, gamma = 0.3, potts_sweeps = 20L,
    theta = 0.3, size = c(3, 8, 12), seed = 3151L)
  set.seed(3152L)
  fit <- suppressWarnings(run_metropolis_MCMC_betas(
    N = 5L, gamma_prior = 0.3, iterations = 4L,
    x_vars = truth$x_vars, y = truth$y, use_data_priors = TRUE,
    dist = "ZINB", theta_start = 0.3, size_start = c(3, 8, 12),
    z_start = truth$z_true, beta_start = truth$beta_true,
    validation_freeze_z = TRUE, validation_freeze_beta = TRUE,
    validation_freeze_gamma = TRUE, validation_freeze_size = TRUE,
    validation_freeze_theta = TRUE, mcse_stop = FALSE,
    z_probability_burnin_arg = 2L))

  expect_true(all(vapply(fit$chains, function(chain)
    all(apply(chain, 2L, function(value) length(unique(value)) == 1L)),
    logical(1))))
  expect_equal(unname(fit$chains[[1L]][1L, ]),
               unname(truth$beta_true[1L, ]))
  expect_true(all(fit$gamma == 0.3))
  expect_true(all(fit$theta == 0.3))
  expect_equal(fit$size,
               matrix(rep(c(3, 8, 12), 5L), nrow = 3L))
  expect_equal(fit$z_checkpoints$iter_1,
               matrix(as.numeric(truth$z_true), 5L, 5L))
  expect_true(all(unlist(fit$sampler_settings[
    grep("^validation_freeze_", names(fit$sampler_settings))])))
})

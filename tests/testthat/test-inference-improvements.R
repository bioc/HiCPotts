test_that("scaled priors preserve original coefficient structure", {
  N <- 4L
  x_vars <- list(
    distance = list(matrix(0:(N * N - 1), N)),
    GC = list(matrix(seq(0, 1, length.out = N * N), N)),
    TES = list(matrix(seq(0, 2, length.out = N * N), N)),
    ACC = list(matrix(seq(0, 3, length.out = N * N), N))
  )
  priors <- make_hicpotts_scaled_priors(x_vars, slope_sd_standardized = 0.5)
  expect_named(priors, paste0("component", 1:3))
  expect_named(priors[[1]], c("meany", "meanx1", "meanx2", "meanx3", "meanx4",
                              "sdy", "sdx1", "sdx2", "sdx3", "sdx4"))
  expect_equal(priors[[1]]$sdx1 * sd(log1p(as.numeric(x_vars$distance[[1]]))), 0.5)
})

test_that("z_final is the unforced allocation and proposal metadata is original scale", {
  set.seed(44)
  N <- 4L
  y <- matrix(c(0, 0, 2, 5, 0, 1, 3, 7,
                0, 2, 4, 8, 0, 1, 6, 9), N, N)
  x_vars <- list(
    distance = list(outer(1:N, 1:N, function(i, j) abs(i - j))),
    GC = list(matrix(seq(0, 1, length.out = N * N), N)),
    TES = list(matrix(seq(0, 2, length.out = N * N), N)),
    ACC = list(matrix(seq(0, 3, length.out = N * N), N))
  )
  fit <- suppressWarnings(run_metropolis_MCMC_betas(
    N, 0.3, 3L, x_vars, y, TRUE, NULL, "Poisson"))

  ## H7: z_final is now the UNFORCED posterior allocation. The forced mask is
  ## still available, but under an explicit name and never as the primary
  ## result. (This block previously asserted the opposite, encoding the defect:
  ## every observed zero was reported as component 1 regardless of the sampled
  ## latent state and regardless of distribution family.)
  forced <- attr(fit$z_final, "zero_forced_classification")
  expect_false(is.null(forced))
  expect_true(all(forced[y == 0] == 1))

  ## The primary allocation must agree with the internal parameter field.
  expect_identical(as.vector(fit$z_final),
                   as.vector(attr(fit$z_final, "z_parameter_final")))

  expect_equal(sum(attr(fit$z_final, "parameter_component_counts")), N * N)
  expect_equal(sum(attr(fit$z_final, "classification_component_counts")), N * N)
  expect_match(attr(fit$z_final, "zero_allocation"), "UNFORCED")
  expect_identical(attr(fit$chains[[1]], "coefficient_scale"),
                   "original log1p-covariate scale")
  expect_named(attr(fit$chains[[1]], "proposal_covariate_sds"),
               c("distance", "GC", "TES", "ACC"))
})

test_that("diagnostics report uncertainty and occupancy", {
  set.seed(1)
  make_fit <- function(offset) {
    chains <- lapply(1:3, function(k) {
      out <- matrix(rnorm(100, offset + k), 20, 5)
      colnames(out) <- paste0("b", 0:4)
      out
    })
    z <- matrix(rep(1:3, length.out = 16), 4)
    attr(z, "parameter_component_counts") <- c(6L, 5L, 5L)
    attr(z, "classification_component_counts") <- c(8L, 4L, 4L)
    list(chains = chains, gamma = runif(20), theta = runif(20),
         size = matrix(rexp(60), 3, 20), z_final = z,
         z_checkpoints = list())
  }
  out <- diagnose_hicpotts_fit(list(make_fit(0), make_fit(0.01)),
                               burnin = 5L, minimum_component_cells = 6L)
  expect_true(all(c("estimate", "CI_lower", "CI_upper", "CI_width",
                    "probability_positive", "ESS", "Rhat") %in%
                  names(out$parameters)))
  expect_equal(nrow(out$component_occupancy), 6L)
  expect_true(any(!out$component_occupancy$sufficient_parameter_cells))
  expect_true(all(is.finite(out$parameters$Rhat)))
})

test_that("initial allocation strategies are valid and do not force zeros", {
  N <- 6L
  y <- matrix(c(rep(0, 12), seq_len(24)), N)
  distance <- outer(seq_len(N), seq_len(N), function(i, j) abs(i - j))
  x_vars <- list(distance = list(distance))
  for (method in c("random", "count_quantile", "distance_adjusted",
                   "noise_anchored_random")) {
    z <- make_hicpotts_initial_z(y, x_vars, method = method, seed = 10)
    expect_identical(dim(z), dim(y))
    expect_true(all(z %in% 1:3))
    expect_equal(length(unique(as.integer(z))), 3L)
  }
  z_random <- make_hicpotts_initial_z(y, x_vars, method = "random", seed = 11)
  expect_true(any(z_random[y == 0] != 1L))
  low_count <- rank(log1p(as.numeric(y)), ties.method = "first") <= ceiling(length(y) / 3)
  for (method in c("distance_adjusted", "noise_anchored_random")) {
    z <- make_hicpotts_initial_z(y, x_vars, method = method, seed = 12)
    expect_true(all(as.integer(z)[low_count] == 1L))
    expect_true(all(as.integer(z)[!low_count] %in% 2:3))
  }
})

test_that("covariate diagnostics flag collinearity", {
  N <- 5L
  base <- matrix(seq(0, 1, length.out = N * N), N)
  x_vars <- list(
    distance = list(base), GC = list(base * 2),
    TES = list(matrix(seq(0, 2, length.out = N * N), N)),
    ACC = list(matrix(rev(seq(0, 3, length.out = N * N)), N)))
  out <- diagnose_hicpotts_covariates(x_vars)
  expect_false(out$reliable)
  expect_true(any(out$pairs$high_correlation))
  expect_true(length(out$warnings) > 0L)
})

test_that("relabel swaps occupancy metadata with components 2 and 3", {
  set.seed(30)
  n_draw <- 20L
  make_chain <- function(beta) matrix(rep(beta, each = n_draw), n_draw, 5) +
    matrix(rnorm(n_draw * 5, sd = 0.01), n_draw, 5)
  z <- matrix(rep(1:3, length.out = 25), 5)
  attr(z, "parameter_component_counts") <- c(component1 = 10L,
    component2 = 4L, component3 = 11L)
  attr(z, "classification_component_counts") <- c(component1 = 12L,
    component2 = 3L, component3 = 10L)
  fit <- list(chains = list(
      make_chain(c(0, 0.2, -0.1, 0.3, 0)),
      make_chain(c(3, 0.2, -0.1, 0.3, 0)),
      make_chain(c(2, -1, 1, 0.8, -0.7))),
    gamma = runif(n_draw), theta = runif(n_draw),
    size = matrix(1, 3, n_draw), z_final = z, z_checkpoints = list())
  relabelled <- relabel_hicpotts(fit)
  expect_equal(unname(attr(relabelled$z_final, "parameter_component_counts")),
               c(10L, 11L, 4L))
  expect_equal(unname(attr(relabelled$z_final, "classification_component_counts")),
               c(12L, 10L, 3L))
})

test_that("multi-chain runner records starts and returns relabelled fits", {
  set.seed(52)
  N <- 4L
  y <- matrix(rpois(N * N, 3), N)
  x_vars <- list(
    distance = list(outer(1:N, 1:N, function(i, j) abs(i - j))),
    GC = list(matrix(runif(N * N), N)),
    TES = list(matrix(runif(N * N), N)),
    ACC = list(matrix(runif(N * N), N)))
  fits <- suppressWarnings(run_hicpotts_chains(
    N = N, iterations = 2L, x_vars = x_vars, y = y,
    dist = "Poisson", seeds = 21:22,
    initialization = c("count_quantile", "distance_adjusted")))
  expect_length(fits, 2L)
  expect_identical(vapply(fits, `[[`, character(1), "initialization_method"),
                   c(chain1 = "count_quantile", chain2 = "distance_adjusted"))
  expect_true(all(vapply(fits, function(x) !is.null(x$relabel_rule), logical(1))))
})

test_that("diagnostics expose explicit reliability gates", {
  set.seed(72)
  make_fit <- function(seed) {
    set.seed(seed)
    chains <- lapply(1:3, function(k) matrix(rnorm(250, k), 50, 5))
    z <- matrix(rep(1:3, length.out = 400), 20)
    attr(z, "parameter_component_counts") <- c(134L, 133L, 133L)
    attr(z, "classification_component_counts") <- c(134L, 133L, 133L)
    for (k in 1:3) attr(chains[[k]], "coefficient_scale") <-
      "original log1p-covariate scale"
    relabel_hicpotts(list(chains = chains, gamma = rnorm(50), theta = runif(50),
      size = matrix(rexp(150), 3, 50), z_final = z, z_checkpoints = list()))
  }
  fits <- lapply(1:4, make_fit)
  out <- diagnose_hicpotts_fit(fits, burnin = 10L, minimum_ess = 1)
  expect_named(out$reliability_flags, c("criterion", "passed", "threshold"))
  expect_true(all(c("split_rhat", "independent_chains", "relabelled") %in%
                  out$reliability_flags$criterion))
  expect_true("beta_movement" %in% out$reliability_flags$criterion)
  expect_true("common_posterior_target" %in%
    out$reliability_flags$criterion)
  # Summary verdicts were removed: reliability_flags is the reported result.
  expect_false(any(c(
    "allocation_reliable", "parameter_reliable", "gamma_reliable",
    "reliable", "reliability_by_domain"
  ) %in% names(out)))
  expect_identical(out$coefficient_scale, "original log1p-covariate scale")
})

test_that("zero-cell components pass the default occupancy threshold", {
  make_fit <- function(seed) {
    set.seed(seed)
    chains <- lapply(1:3, function(k) matrix(rnorm(250, k), 50, 5))
    for (k in 1:3) {
      attr(chains[[k]], "coefficient_scale") <-
        "original log1p-covariate scale"
    }
    z <- matrix(rep(1:2, length.out = 400L), 20L)
    attr(z, "parameter_component_counts") <- c(200L, 200L, 0L)
    attr(z, "classification_component_counts") <- c(200L, 200L, 0L)
    list(
      chains = chains, gamma = runif(50), theta = runif(50),
      size = matrix(rexp(150), 3L, 50L), z_final = z,
      z_checkpoints = list(), relabel_rule = list(applied = TRUE)
    )
  }
  fits <- lapply(1:4, make_fit)
  out <- diagnose_hicpotts_fit(
    fits, burnin = 10L, minimum_ess = 1, relabel = FALSE
  )
  component3 <- out$component_occupancy$component == 3L
  expect_true(all(out$component_occupancy$parameter_cells[component3] == 0L))
  expect_true(all(
    out$component_occupancy$sufficient_parameter_cells[component3]
  ))
  occupancy_flag <- out$reliability_flags$criterion == "component_occupancy"
  expect_true(out$reliability_flags$passed[occupancy_flag])
})

test_that("likelihood-informed initialization is valid and reproducible", {
  set.seed(901)
  N <- 6L
  y <- matrix(rnbinom(N * N, mu = rep(c(1, 6, 15), length.out = N * N), size = 4), N)
  x_vars <- list(
    distance = list(outer(1:N, 1:N, function(i, j) abs(i - j))),
    GC = list(matrix(runif(N * N), N)),
    TES = list(matrix(runif(N * N), N)),
    ACC = list(matrix(runif(N * N), N)))
  z1 <- make_hicpotts_initial_z(y, x_vars, method = "likelihood_informed",
    seed = 17L, dist = "ZINB", size_start = c(3, 8, 20))
  z2 <- make_hicpotts_initial_z(y, x_vars, method = "likelihood_informed",
    seed = 17L, dist = "ZINB", size_start = c(3, 8, 20))
  expect_identical(z1, z2)
  expect_identical(dim(z1), c(N, N))
  expect_true(all(z1 %in% 1:3))
  expect_true(all(tabulate(z1, nbins = 3L) > 0L))
})

test_that("likelihood-informed initialization skips unstable candidate fits", {
  path <- system.file("extdata", "test_data2.csv", package = "HiCPotts")
  if (!nzchar(path)) {
    path <- testthat::test_path("..", "..", "inst", "extdata",
                               "test_data2.csv")
  }
  upper <- utils::read.csv(path, check.names = FALSE)
  bins <- sort(unique(c(upper$start, upper[["start.j."]])))
  bins <- bins[seq_len(20L)]
  upper <- upper[
    upper$start %in% bins & upper[["start.j."]] %in% bins,
    , drop = FALSE
  ]
  symmetric_matrix <- function(value) {
    out <- matrix(NA_real_, 20L, 20L)
    i <- match(upper$start, bins)
    j <- match(upper[["start.j."]], bins)
    out[cbind(i, j)] <- value
    out[cbind(j, i)] <- value
    out
  }
  y <- symmetric_matrix(upper$interactions)
  x_vars <- list(
    distance = list(abs(outer(bins, bins, "-"))),
    GC = list(symmetric_matrix(upper$GC)),
    TES = list(symmetric_matrix(upper$TES)),
    ACC = list(symmetric_matrix(upper$ACC))
  )
  expect_no_error({
    z <- make_hicpotts_initial_z(
      y, x_vars, method = "likelihood_informed", seed = 20260822L,
      dist = "ZINB", size_start = c(1, 1, 1)
    )
  })
  expect_identical(dim(z), c(20L, 20L))
  expect_true(all(z %in% 1:3))
  expect_true(all(tabulate(z, nbins = 3L) > 0L))
})

test_that("20,000 is the robust default and MCSE metadata is retained", {
  set.seed(902)
  N <- 3L
  y <- matrix(rpois(N * N, 2), N)
  x_vars <- list(
    distance = list(outer(1:N, 1:N, function(i, j) abs(i - j))),
    GC = list(matrix(runif(N * N), N)),
    TES = list(matrix(runif(N * N), N)),
    ACC = list(matrix(runif(N * N), N)))
  fit <- suppressWarnings(run_hicpotts_chains(N, iterations = 2L,
    x_vars = x_vars, y = y, dist = "Poisson", seeds = 91L,
    initialization = "count_quantile", relabel = FALSE)[[1L]])
  expect_identical(names(fit)[1:6],
    c("chains", "gamma", "theta", "size", "z_final", "z_checkpoints"))
  expect_identical(attr(fit$gamma, "iterations_requested"), 2L)
  expect_identical(formals(fit_hicpotts_robust)$iterations, 20000L)
  robust_starts <- eval(formals(fit_hicpotts_robust)$initialization)
  expect_identical(robust_starts, c("likelihood_informed", "likelihood_informed",
    "distance_adjusted", "noise_anchored_random"))
  expect_identical(attr(fit$gamma, "iterations_completed"), 2L)
  expect_true(isTRUE(attr(fit$gamma, "mcse_stopping_enabled")))

  stopped <- suppressWarnings(run_hicpotts_chains(N, iterations = 20001L,
    x_vars = x_vars, y = y, dist = "Poisson", seeds = 92L,
    initialization = "count_quantile", relabel = FALSE,
    mcse_min_iterations = 100L, mcse_check_interval = 50L,
    mcse_relative_threshold = 100)[[1L]])
  completed <- attr(stopped$gamma, "iterations_completed")
  expect_identical(attr(stopped$gamma, "iterations_requested"), 20001L)
  expect_lt(completed, 20001L)
  expect_identical(nrow(stopped$chains[[1L]]), completed + 1L)
  expect_true(isTRUE(attr(stopped$gamma, "mcse_converged")))
  expect_true(all(stopped$beta_mixing$retained_acceptances > 0))
})

test_that("legacy biological relabelling resolves internal 2/3 permutations", {
  set.seed(903)
  make_chain <- function(beta)
    matrix(rep(beta, each = 20L), 20L, 5L) +
      matrix(rnorm(100, sd = 0.01), 20L, 5L)
  make_fit <- function(z, swapped = FALSE) {
    b1 <- c(0, 0.2, -0.1, 0.3, 0)
    signal <- c(2, -1, 1, 0.8, -0.7)
    false_signal <- c(3, 0.2, -0.1, 0.3, 0)
    elevated <- if (swapped) list(false_signal, signal) else
      list(signal, false_signal)
    chains <- list(make_chain(b1), make_chain(elevated[[1L]]),
                   make_chain(elevated[[2L]]))
    attr(z, "z_parameter_final") <- z
    list(chains = chains, gamma = runif(20), theta = runif(20),
      size = matrix(2, 3, 20), z_final = z, z_checkpoints = list(),
      noise_relationship = list(enabled = FALSE, link_sd = 0.5,
        order_strength = 10, order_width = 0.5))
  }
  z1 <- matrix(rep(c(1L, 2L, 3L, 2L, 3L), length.out = 25), 5)
  z2 <- z1
  z2[z1 == 2L] <- 3L
  z2[z1 == 3L] <- 2L
  aligned <- relabel_hicpotts(list(make_fit(z1), make_fit(z2, swapped = TRUE)))
  expect_identical(as.integer(aligned[[1L]]$z_final),
                   as.integer(aligned[[2L]]$z_final))
  expect_true(all(vapply(aligned, function(x)
    is.finite(x$noise_relationship_probability), logical(1))))
  expect_true(all(vapply(aligned, function(x)
    identical(x$noise_relationship$threshold_policy,
      "none in the package; downstream users choose any probability threshold"),
    logical(1))))
  expect_true(all(vapply(aligned, function(x)
    is.null(x$noise_relationship_identified), logical(1))))
  expect_true(all(vapply(aligned, function(x)
    is.null(x$allocation_reference_chain), logical(1))))
})

test_that("posterior predictive checks use unforced internal Z", {
  set.seed(84)
  N <- 4L
  n_draw <- 20L
  make_chain <- function(intercept) {
    out <- cbind(rep(intercept, n_draw), matrix(0, n_draw, 4))
    attr(out, "coefficient_scale") <- "original log1p-covariate scale"
    out
  }
  z_parameter <- matrix(rep(1:3, length.out = N * N), N)
  z_reported <- z_parameter
  y <- matrix(rnbinom(N * N, mu = 4, size = 3), N)
  z_reported[y == 0] <- 1L
  attr(z_reported, "parameter_component_counts") <- tabulate(z_parameter, 3L)
  attr(z_reported, "classification_component_counts") <- tabulate(z_reported, 3L)
  fit <- list(chains = list(make_chain(log(2)), make_chain(log(4)),
                            make_chain(log(8))),
    gamma = rep(0.4, n_draw), theta = rep(0.2, n_draw),
    size = matrix(3, 3, n_draw), z_final = z_reported,
    z_parameter_final = z_parameter, z_checkpoints = list())
  x_vars <- list(
    distance = list(outer(1:N, 1:N, function(i, j) abs(i - j))),
    GC = list(matrix(runif(N * N), N)),
    TES = list(matrix(runif(N * N), N)),
    ACC = list(matrix(runif(N * N), N)))
  predictive <- posterior_predictive_hicpotts(
    list(fit, fit), x_vars, y, dist = "NB", burnin = 5L,
    n_rep = 5L, seed = 2L)
  expect_identical(predictive$conditioning, "final internal unforced Z")
  expect_equal(nrow(predictive$discrepancies), 5L)
  expect_true(all(is.finite(as.matrix(predictive$summary[, -1L]))))

  families <- compare_hicpotts_families(
    list(NB = list(fit, fit), ZINB = list(fit, fit)),
    x_vars, y, burnin = 5L, n_rep = 3L, seed = 3L)
  expect_setequal(unique(families$comparison$family), c("NB", "ZINB"))
  expect_true(all(c("zero_fraction", "near_diagonal_zero_fraction",
                    "variance_to_mean") %in% families$comparison$metric))
})

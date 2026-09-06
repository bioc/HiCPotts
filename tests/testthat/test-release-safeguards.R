make_1310_diagnostic_fit <- function(gamma, seed = 1L) {
  set.seed(seed)
  n <- length(gamma)
  chains <- lapply(1:3, function(k) {
    out <- matrix(rnorm(n * 5L, k, 0.2), n, 5L)
    attr(out, "coefficient_scale") <- "original log1p-covariate scale"
    out
  })
  z <- matrix(rep(1:3, length.out = 144), 12L)
  attr(z, "parameter_component_counts") <- c(48L, 48L, 48L)
  attr(z, "classification_component_counts") <- c(48L, 48L, 48L)
  list(
    chains = chains,
    gamma = gamma,
    theta = runif(n, 0.2, 0.5),
    size = matrix(rexp(3L * n, rate = 0.2), 3L, n),
    z_final = z,
    z_checkpoints = list())
}

test_that("gamma boundary and movement gates reject a stuck chain", {
  fits <- lapply(1:2, function(seed)
    make_1310_diagnostic_fit(rep(0.998461748985909, 100L), seed))
  diagnosed <- diagnose_hicpotts_fit(
    fits,
    burnin = 50L,
    minimum_component_cells = 1L,
    minimum_ess = 1,
    minimum_chains = 2L,
    relabel = FALSE)

  flags <- stats::setNames(
    diagnosed$reliability_flags$passed,
    diagnosed$reliability_flags$criterion)
  expect_false(flags[["gamma_not_boundary"]])
  expect_false(flags[["gamma_movement"]])
  expect_equal(diagnosed$gamma_diagnostics$boundary_fraction, c(1, 1))
  expect_equal(diagnosed$gamma_diagnostics$unique_draws, c(1L, 1L))
  expect_false(all(diagnosed$reliability_flags$passed))
})

test_that("safe parameter reporting rejects failed gates by default", {
  fits <- lapply(1:2, function(seed)
    make_1310_diagnostic_fit(rep(0.999, 100L), seed))
  diagnosed <- diagnose_hicpotts_fit(
    fits,
    burnin = 50L,
    minimum_component_cells = 1L,
    minimum_ess = 1,
    minimum_chains = 2L,
    relabel = FALSE)

  # gamma sitting at a boundary no longer blocks the parameter table: it does
  # not make a regression coefficient wrong. Rejection now cites either an
  # unmet precondition or an unresolved parameter.
  expect_error(
    summarise_hicpotts_parameters(diagnosed),
    "Parameter summary rejected")
  provisional <- summarise_hicpotts_parameters(
    diagnosed, require_reliable = FALSE)
  expect_s3_class(provisional, "hicpotts_parameter_summary")
  expect_false(attr(provisional, "resolved"))
  expect_true("gamma_not_boundary" %in%
    attr(provisional, "failed_reliability_gates"))
})

make_1310_classification_fit <- function() {
  probabilities <- list(
    component1 = matrix(c(.80, .10, .10, .10), 2L),
    component2 = matrix(c(.10, .75, .20, .15), 2L),
    component3 = matrix(c(.10, .15, .70, .75), 2L))
  z <- matrix(max.col(cbind(
    as.vector(probabilities[[1L]]),
    as.vector(probabilities[[2L]]),
    as.vector(probabilities[[3L]])), ties.method = "first"), 2L)
  list(
    chains = list(
      cbind(rep(-1, 40L), matrix(0, 40L, 4L)),
      cbind(rep(0.5, 40L), matrix(0, 40L, 4L)),
      cbind(rep(1.5, 40L), matrix(0, 40L, 4L))),
    gamma = rep(0.2, 40L),
    theta = rep(0.1, 40L),
    size = matrix(c(3, 8, 20), 3L, 40L),
    z_final = z,
    z_checkpoints = list(),
    z_probabilities = probabilities,
    z_probability_draws = 200L,
    pair_weighting = list(symmetric_input = FALSE))
}

test_that("classification comparison keeps sampled z official and three-way", {
  dat <- data.frame(
    start = c(0, 0, 0, 0),
    end = c(0, 1, 1, 2),
    interactions = c(0, 2, 5, 7),
    GC = c(.1, .2, .3, .4),
    TES = c(.2, .3, .4, .5),
    ACC = c(.3, .4, .5, .6))
  out <- compare_hicpotts_classifications(
    make_1310_classification_fit(),
    data = dat,
    N = 2L,
    dist = "Poisson",
    relabel = FALSE,
    min_draws = 10L,
    reflect = "never",
    n_draws = 10L,
    potts_iterations = 1L)

  expect_s3_class(out, "hicpotts_classification_sensitivity")
  expect_match(out$summary$official_source, "latent-state MCMC frequencies")
  expect_match(out$summary$sensitivity_source, "secondary parameter-integrated")
  expect_identical(levels(out$official$classification),
                   c("noise", "signal", "false signal"))
  expect_identical(levels(out$parameter_sensitivity$classification),
                   c("noise", "signal", "false signal"))
  expect_false(anyNA(out$official$classification))
  expect_false(anyNA(out$parameter_sensitivity$classification))
  expect_equal(rowSums(out$parameter_sensitivity[, c("prob1", "prob2", "prob3")]),
               rep(1, 4))
  expect_equal(out$summary$cells_compared, 4L)
  expect_identical(dim(out$confusion), c(3L, 3L))
})

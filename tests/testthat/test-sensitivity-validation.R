## Items 12 and 13: the structured sensitivity workflow and the known-truth
## validation suite. Kept small so the suite stays fast; the scientific runs
## are performed by the user at production settings.

test_that("the truth simulator matches the component-1/3 noise relationship", {
  truth <- simulate_hicpotts_truth(N = 16L, seed = 2L)
  expect_identical(dim(truth$y), c(16L, 16L))
  expect_true(all(truth$z_true %in% 1:3))
  expect_true(all(c(1L, 2L, 3L) %in% unique(as.integer(truth$z_true))))
  expect_named(truth$x_vars, c("distance", "GC", "TES", "ACC"))

  expect_equal(truth$beta_true[3L, -1L], truth$beta_true[1L, -1L])
  expect_gt(truth$beta_true[3L, 1L], truth$beta_true[1L, 1L])
  expect_gt(sqrt(sum((truth$beta_true[2L, -1L] -
                        truth$beta_true[1L, -1L])^2)), 1)
})

test_that("validation reports accuracy, calibration and stability", {
  truth <- simulate_hicpotts_truth(N = 12L, seed = 4L)
  v <- validate_hicpotts_simulation(truth, iterations = 300L, n_chains = 2L,
                                    min_draws = 10L, seed = 2L)
  expect_s3_class(v, "hicpotts_validation")
  expect_true(v$accuracy >= 0 && v$accuracy <= 1)
  expect_identical(dim(v$confusion), c(3L, 3L))
  expect_identical(colnames(v$per_class), c("noise", "signal", "false signal"))
  expect_true(is.finite(v$expected_calibration_error))
  expect_true(all(c("distance_signal_to_noise",
                    "distance_false_signal_to_noise",
                    "false_signal_intercept_gap", "correct_relationship",
                    "relationship_probability",
                    "threshold_policy") %in%
                  names(v$relationship_recovery)))
  expect_false(is.null(v$stability))
  expect_output(print(v), "known-truth validation")
})

test_that("sensitivity reports movement across setting families", {
  set.seed(6); N <- 10L
  y <- matrix(rpois(N * N, 5), N, N)
  x_vars <- list(distance = list(abs(row(y) - col(y))),
                 GC = list(matrix(runif(N * N), N, N)),
                 TES = list(matrix(runif(N * N), N, N)),
                 ACC = list(matrix(runif(N * N), N, N)))
  s <- hicpotts_sensitivity(N = N, x_vars = x_vars, y = y, dist = "Poisson",
                            iterations = 200L,
                            scenarios = c("relationship", "abc"),
                            seed = 4L, min_draws = 10L)
  expect_s3_class(s, "hicpotts_sensitivity")
  expect_true(all(c("family", "scenario", "label_agreement",
                    "mean_probability_shift") %in% names(s$summary)))
  expect_true(all(s$summary$label_agreement >= 0 & s$summary$label_agreement <= 1))
  expect_true(all(s$summary$mean_probability_shift >= 0))
  expect_true(any(s$summary$family == "relationship"))
  expect_true(any(s$summary$family == "abc"))
  expect_true(any(grepl("prior disabled", s$summary$scenario)))
  expect_output(print(s), "sensitivity analysis")
})

test_that("sensitivity does not erase biological 2/3 disagreement by swapping", {
  set.seed(7); N <- 8L
  y <- matrix(rpois(N * N, 6), N, N)
  x_vars <- list(distance = list(abs(row(y) - col(y))),
                 GC = list(matrix(runif(N * N), N, N)),
                 TES = list(matrix(runif(N * N), N, N)),
                 ACC = list(matrix(runif(N * N), N, N)))
  s <- hicpotts_sensitivity(N = N, x_vars = x_vars, y = y, dist = "Poisson",
                            iterations = 150L, scenarios = "relationship",
                            seed = 5L, min_draws = 10L)
  ## Biological labels have already been assigned by relabelling; an
  ## agreement-maximising second swap would hide signal/false-signal changes.
  expect_type(s$summary$permutation_aligned, "logical")
  expect_false(anyNA(s$summary$permutation_aligned))
  expect_false(any(s$summary$permutation_aligned))
})

## Draw counts are well above classify_hicpotts()'s min_draws floor. The
## classifier refuses too few post-burn-in latent-state draws and does not fall
## back to the single final state.
make_classification_fit <- function(probabilities, draws = 200L,
                                    intercepts = c(0, 1, 2)) {
  lattice_dim <- dim(probabilities[[1L]])
  z <- matrix(max.col(cbind(
    as.vector(probabilities[[1L]]),
    as.vector(probabilities[[2L]]),
    as.vector(probabilities[[3L]])), ties.method = "first"),
    nrow = lattice_dim[1L], ncol = lattice_dim[2L])
  list(
    chains = lapply(intercepts, function(value)
      cbind(rep(value, 6L), matrix(0, 6L, 4L))),
    gamma = rep(0.2, 6L), theta = rep(0.1, 6L),
    size = matrix(2, 3L, 6L), z_final = z, z_checkpoints = list(),
    z_probabilities = probabilities, z_probability_draws = draws)
}

test_that("classifier pools latent-state frequencies by their draw counts", {
  p1 <- list(
    component1 = matrix(c(.8, .2, .1, .1), 2L),
    component2 = matrix(c(.1, .7, .7, .2), 2L),
    component3 = matrix(c(.1, .1, .2, .7), 2L))
  p2 <- list(
    component1 = matrix(c(.6, .2, .1, .1), 2L),
    component2 = matrix(c(.2, .6, .8, .3), 2L),
    component3 = matrix(c(.2, .2, .1, .6), 2L))
  out <- classify_hicpotts(
    list(make_classification_fit(p1, 300L), make_classification_fit(p2, 100L)),
    relabel = FALSE)
  expect_equal(out$prob1, as.vector((3 * p1[[1L]] + p2[[1L]]) / 4))
  expect_equal(rowSums(out[, c("prob1", "prob2", "prob3")]), rep(1, 4))
  expect_identical(attr(out, "posterior_state_draws"), 400)
  expect_identical(attr(out, "chains_combined"), 2L)
})

test_that("classifier assigns every cell to its MAP component", {
  probs <- list(
    component1 = matrix(c(.70, .45, .20), 1L),
    component2 = matrix(c(.20, .40, .55), 1L),
    component3 = matrix(c(.10, .15, .25), 1L))
  out <- classify_hicpotts(make_classification_fit(probs), relabel = FALSE)
  ## The classification is a factor with exactly the three biological levels,
  ## so no fourth category can be introduced downstream.
  expect_s3_class(out$classification, "factor")
  expect_identical(levels(out$classification),
                   c("noise", "signal", "false signal"))
  expect_identical(as.character(out$classification),
                   c("noise", "noise", "signal"))
  expect_identical(out$map_component, c(1L, 1L, 2L))
  expect_true(all(out$normalized_entropy >= 0 & out$normalized_entropy <= 1))
  expect_false(anyNA(out$classification))
})

test_that("classification refuses too few post-burn-in draws and never falls back", {
  probs <- list(
    component1 = matrix(c(.70, .45, .20), 1L),
    component2 = matrix(c(.20, .40, .55), 1L),
    component3 = matrix(c(.10, .15, .25), 1L))
  expect_error(
    classify_hicpotts(make_classification_fit(probs, 10L), relabel = FALSE),
    "below the required minimum")
  ## ...but it is configurable.
  expect_s3_class(
    classify_hicpotts(make_classification_fit(probs, 10L), relabel = FALSE,
                      min_draws = 5L),
    "hicpotts_classification")
  ## A fit that cannot report its draw count is rejected outright rather than
  ## classified from the single final latent state.
  bad <- make_classification_fit(probs, 200L)
  bad$z_probability_draws <- NULL
  expect_error(classify_hicpotts(bad, relabel = FALSE),
               "z_probability_draws")
})

test_that("reflection averaging gives (i,j) and (j,i) the same classification", {
  set.seed(9)
  N <- 6L
  a <- matrix(runif(N * N), N); b <- matrix(runif(N * N), N)
  c3 <- matrix(runif(N * N), N)
  tot <- a + b + c3
  probs <- list(component1 = a / tot, component2 = b / tot, component3 = c3 / tot)
  f <- make_classification_fit(probs, 500L)
  f$pair_weighting <- list(symmetric_input = TRUE, offdiagonal_weight = 0.5)

  out <- classify_hicpotts(f, relabel = FALSE, reflect = "always")
  expect_true(attr(out, "reflection_averaged"))
  m <- matrix(as.integer(out$map_component), N, N)
  expect_identical(m, t(m))
  pm <- matrix(out$prob2, N, N)
  expect_equal(pm, t(pm))

  ## Disabled explicitly -> asymmetric probabilities are left alone.
  off <- classify_hicpotts(f, relabel = FALSE, reflect = "never")
  expect_false(attr(off, "reflection_averaged"))
})

test_that("relabel keeps sampled membership probabilities aligned", {
  probs <- list(
    component1 = matrix(.1, 2L, 2L),
    component2 = matrix(.2, 2L, 2L),
    component3 = matrix(.7, 2L, 2L))
  fit <- make_classification_fit(probs, intercepts = c(0, 3, 1))
  relabelled <- relabel_hicpotts(fit)
  expect_equal(relabelled$z_probabilities[[2L]], probs[[3L]])
  expect_equal(relabelled$z_probabilities[[3L]], probs[[2L]])
  expect_named(relabelled$z_probabilities, paste0("component", 1:3))
})

test_that("sampler records normalized post-burn-in membership frequencies", {
  set.seed(2401)
  N <- 3L
  y <- matrix(rpois(N * N, 3), N)
  x_vars <- list(
    distance = list(outer(1:N, 1:N, function(i, j) abs(i - j))),
    GC = list(matrix(runif(N * N), N)),
    TES = list(matrix(runif(N * N), N)),
    ACC = list(matrix(runif(N * N), N)))
  fit <- suppressWarnings(run_metropolis_MCMC_betas(
    N, 0.3, 4L, x_vars, y, TRUE, NULL, "Poisson",
    mcse_stop = FALSE, z_probability_burnin_arg = 1L))
  expect_length(fit$z_probabilities, 3L)
  expect_true(all(vapply(fit$z_probabilities,
    function(x) identical(dim(x), c(N, N)), logical(1))))
  expect_equal(Reduce(`+`, fit$z_probabilities), matrix(1, N, N))
  expect_identical(fit$z_probability_burnin, 1L)
  expect_identical(fit$z_probability_draws, 3L)
})

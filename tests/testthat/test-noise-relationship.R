make_relationship_fit <- function(b2, b3, n = 400L,
                                  b1 = c(0, 0.2, -0.1, 0.3, 0),
                                  enabled = TRUE) {
  make_chain <- function(beta) {
    matrix(rep(beta, each = n), nrow = n) +
      matrix(stats::rnorm(n * length(beta), sd = 0.005), nrow = n)
  }
  chains <- list(make_chain(b1), make_chain(b2), make_chain(b3))
  for (k in seq_len(3L))
    attr(chains[[k]], "proposal_covariate_sds") <-
      c(distance = 0.8, GC = 0.5, TES = 0.6, ACC = 0.7)
  z <- matrix(c(1L, 2L, 3L, 2L), 2L)
  attr(z, "parameter_component_counts") <- c(1L, 2L, 1L)
  list(
    chains = chains, size = matrix(2, 3L, n), z_final = z,
    z_checkpoints = list(), noise_relationship = list(
      enabled = enabled, link_sd = 0.5, order_strength = 10,
      order_width = 0.5))
}

test_that("legacy fits receive one biological component-2/3 swap", {
  set.seed(11)
  false_signal <- c(3.0, 0.2, -0.1, 0.3, 0)
  true_signal <- c(2.0, -1.2, 1.0, 0.9, -0.8)
  out <- relabel_hicpotts(make_relationship_fit(
    false_signal, true_signal, enabled = FALSE))

  expect_identical(
    out$relabel_basis,
    paste0("all four standardised covariate slopes relative to component 1, ",
           "plus elevated-noise intercept preference"))
  expect_equal(unname(colMeans(out$chains[[3]])), false_signal,
               tolerance = 0.02)
  expect_equal(unname(colMeans(out$chains[[2]])), true_signal,
               tolerance = 0.02)
  expect_null(out$noise_relationship_identified)
  expect_gte(out$noise_relationship_probability, 0.95)
  expect_match(out$noise_relationship$threshold_policy, "none in the package")
})

test_that("an already biologically correct component order is retained", {
  set.seed(12)
  true_signal <- c(2.0, -1.2, 1.0, 0.9, -0.8)
  false_signal <- c(3.0, 0.2, -0.1, 0.3, 0)
  out <- relabel_hicpotts(make_relationship_fit(true_signal, false_signal))
  expect_equal(unname(colMeans(out$chains[[3]])), false_signal,
               tolerance = 0.02)
  expect_identical(unname(out$label_permutation[1L, ]), 1:3)
})

test_that("component 1 is never reordered and population cannot override biology", {
  set.seed(13)
  b1 <- c(-2, 0.4, -0.2, 0.1, 0.2)
  false_signal <- c(2.5, 0.4, -0.2, 0.1, 0.2)
  true_signal <- c(1.5, -1, 1, 1, -1)
  fit <- make_relationship_fit(
    false_signal, true_signal, b1 = b1, enabled = FALSE)
  fit$z_final[,] <- 3L
  out <- relabel_hicpotts(fit)
  expect_equal(unname(colMeans(out$chains[[1]])), b1, tolerance = 0.02)
  expect_equal(unname(colMeans(out$chains[[3]])), false_signal,
               tolerance = 0.02)
})

test_that("component 2 is absent from the coupled prior", {
  b1 <- c(0, 0.2, -0.1, 0.3, 0)
  b2 <- c(2, -1, 1, 0.8, -0.7)
  another_b2 <- c(-5, 4, -3, 2, -1)
  b3 <- c(3, 0.2, -0.1, 0.3, 0)
  sds <- c(0.8, 0.5, 0.6, 0.7)
  first <- HiCPotts:::.hicpotts_noise_relationship_logprior_cpp(
    b1, b2, b3, sds)
  second <- HiCPotts:::.hicpotts_noise_relationship_logprior_cpp(
    b1, another_b2, b3, sds)
  expect_equal(first, second, tolerance = 1e-12)
})

test_that("the relationship prior favours slope matching plus elevation", {
  b1 <- c(0, 0.2, -0.1, 0.3, 0)
  unrelated2 <- c(2, -1, 1, 0.8, -0.7)
  unrelated3 <- c(3, 1.1, -1.0, -0.9, 0.8)
  matched3 <- c(3, 0.2, -0.1, 0.3, 0)
  sds <- c(0.8, 0.5, 0.6, 0.7)
  without_match <- HiCPotts:::.hicpotts_noise_relationship_logprior_cpp(
    b1, unrelated2, unrelated3, sds)
  with_match <- HiCPotts:::.hicpotts_noise_relationship_logprior_cpp(
    b1, unrelated2, matched3, sds)
  below_noise <- matched3
  below_noise[1L] <- -2
  below <- HiCPotts:::.hicpotts_noise_relationship_logprior_cpp(
    b1, unrelated2, below_noise, sds)
  expect_gt(with_match, without_match)
  expect_gt(with_match, below)
})

test_that("a current fit in the opposite branch is globally relabelled", {
  set.seed(14)
  false_signal <- c(3.0, 0.2, -0.1, 0.3, 0)
  true_signal <- c(2.0, -1.2, 1.0, 0.9, -0.8)
  out <- relabel_hicpotts(make_relationship_fit(false_signal, true_signal))
  expect_identical(unname(out$label_permutation[1L, ]), c(1L, 3L, 2L))
  expect_null(out$noise_relationship_identified)
  expect_gt(out$noise_relationship_probability, 0.95)
  expect_match(out$noise_relationship$threshold_policy, "none in the package")
  expect_equal(unname(colMeans(out$chains[[3]])), false_signal,
               tolerance = 0.02)
})

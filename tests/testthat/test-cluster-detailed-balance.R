test_that("component-2/3 cluster move balances a zero-count boundary", {
  ## A is a uniform component-2 field. B differs only at one corner. With
  ## identical component-2/3 emissions their posterior ratio is exp(-2*gamma),
  ## because B breaks two equal-label Potts edges. The production transition
  ## is exercised through the internal hook, not reimplemented in this test.
  set.seed(1401)
  N <- 2L
  y <- matrix(0, N, N)
  y[1L, 1L] <- 5
  zero <- matrix(0, N, N)
  beta <- matrix(0, 3L, 5L)
  beta[, 1L] <- log(3)
  sizes <- rep(7, 3L)
  gamma <- 0.7
  repetitions <- 100000L
  state_a <- matrix(2, N, N)
  state_b <- state_a
  state_b[1L, 1L] <- 3

  from_a <- HiCPotts:::.hicpotts_cluster_transition_counts(
    state_a, y, zero, zero, zero, zero, beta, sizes,
    gamma, "NB", repetitions)
  from_b <- HiCPotts:::.hicpotts_cluster_transition_counts(
    state_b, y, zero, zero, zero, zero, beta, sizes,
    gamma, "NB", repetitions)
  probability_a_to_b <- from_a[2L] / repetitions
  probability_b_to_a <- from_b[1L] / repetitions
  flow_ratio <- probability_a_to_b /
    (exp(-2 * gamma) * probability_b_to_a)

  expect_equal(flow_ratio, 1, tolerance = 0.04)
})

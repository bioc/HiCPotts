## Item 7: mode detection uses soft distances between z_probabilities and
## preserves the biological component 2/3 distinction, rather than comparing
## single final allocation matrices.

soft_fit <- function(p1, p2, p3) {
  list(z_probabilities = list(component1 = p1, component2 = p2, component3 = p3))
}

test_that("soft distances are used when z_probabilities are present", {
  set.seed(31)
  N <- 8L
  mk <- function(seed) {
    set.seed(seed)
    a <- matrix(runif(N * N), N); b <- matrix(runif(N * N), N); c3 <- matrix(runif(N * N), N)
    s <- a + b + c3
    soft_fit(a / s, b / s, c3 / s)
  }
  res <- .hicpotts_mode_consensus(list(mk(1), mk(1), mk(2)),
                                  agreement_threshold = 0.8)
  expect_match(res$distance_basis, "soft z_probabilities")
  ## Identical fits must agree exactly.
  expect_equal(res$agreement[1, 2], 1)
  expect_true(all(res$agreement >= 0 & res$agreement <= 1))
})

test_that("a component 2/3 exchange is a biological mode difference", {
  N <- 6L
  low <- matrix(0.01, N, N)
  high <- matrix(0.98, N, N)
  base <- soft_fit(low, high, low)
  flipped <- soft_fit(low, low, high)

  res <- .hicpotts_mode_consensus(list(base, flipped), agreement_threshold = 0.8)
  expect_lt(res$agreement[1, 2], 0.1)
  expect_false(res$swap23[1, 2])
  expect_false(res$component[1] == res$component[2])
  expect_identical(res$selected, 1:2)
  expect_false(res$consensus_available)
  expect_false(res$multimodal)
})

test_that("genuinely different modes stay separate", {
  N <- 6L
  n <- N * N
  confident <- function(k) {
    m <- matrix(0.01, N, N); mm <- list(m, m, m)
    mm[[k]] <- matrix(0.98, N, N)
    soft_fit(mm[[1]], mm[[2]], mm[[3]])
  }
  res <- .hicpotts_mode_consensus(
    list(confident(1), confident(1), confident(2), confident(2)),
    agreement_threshold = 0.9, minimum_mode_chains = 2L)
  ## Two replicated, mutually distinct modes -> multimodal, not "consensus".
  expect_true(res$multimodal)
  expect_gte(res$competing_modes, 1L)
})

test_that("soft agreement degrades smoothly, unlike hard label matching", {
  N <- 4L
  near <- soft_fit(matrix(0.51, N, N), matrix(0.49, N, N), matrix(0, N, N))
  flip <- soft_fit(matrix(0.49, N, N), matrix(0.51, N, N), matrix(0, N, N))
  res <- .hicpotts_mode_consensus(list(near, flip), agreement_threshold = 0.8)
  ## Every argmax differs, so hard matching would score 0; the soft distance
  ## correctly reports these as almost identical.
  expect_gt(res$agreement[1, 2], 0.95)
})

test_that("fits without z_probabilities fall back to hard allocations", {
  z <- matrix(rep(1:3, length.out = 100), 10)
  res <- .hicpotts_mode_consensus(list(list(z_final = z), list(z_final = z)),
                                  agreement_threshold = 0.8)
  expect_match(res$distance_basis, "hard final allocations")
  expect_equal(res$agreement[1, 2], 1)
})

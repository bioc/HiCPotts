test_that("mode consensus selects a coherent allocation group", {
  set.seed(9901)
  z1 <- matrix(rep(1:3, length.out = 100), 10)
  z2 <- z1
  z3 <- z1
  z3[sample.int(length(z3), 10)] <- sample.int(3L, 10, replace = TRUE)
  z4 <- matrix(sample.int(3L, 100, replace = TRUE), 10)
  fits <- lapply(list(z1, z2, z3, z4), function(z) list(z_final = z))

  result <- .hicpotts_mode_consensus(
    fits, agreement_threshold = 0.8, minimum_mode_chains = 2L)

  expect_identical(result$selected, 1:3)
  expect_identical(result$excluded, 4L)
  expect_true(result$consensus_available)
  expect_true(all(result$agreement[1:3, 1:3] >= 0.8))
})

test_that("mode consensus does not discard chains without a replicated mode", {
  set.seed(9902)
  fits <- lapply(1:3, function(i)
    list(z_final = matrix(sample.int(3L, 100, replace = TRUE), 10)))
  result <- .hicpotts_mode_consensus(
    fits, agreement_threshold = 0.9, minimum_mode_chains = 2L)
  expect_identical(result$selected, 1:3)
  expect_length(result$excluded, 0L)
  expect_false(result$consensus_available)
})

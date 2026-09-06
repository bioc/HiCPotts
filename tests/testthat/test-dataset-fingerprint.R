## Item 8: fingerprints and coordinate keys stop unrelated datasets being
## pooled and shuffled classification data being silently misaligned.

test_that("a fingerprint is order-sensitive, not just content-sensitive", {
  set.seed(3); N <- 6L
  y <- matrix(rpois(N * N, 5), N, N)
  shuffled <- matrix(sample(as.vector(y)), N, N)
  a <- hicpotts_fingerprint(y)
  b <- hicpotts_fingerprint(shuffled)
  expect_true(isTRUE(hicpotts_fingerprints_agree(a, a)))
  expect_false(isTRUE(hicpotts_fingerprints_agree(a, b)))
})

test_that("different datasets and different lattice sizes are detected", {
  set.seed(4)
  y1 <- matrix(rpois(36, 5), 6L)
  y2 <- matrix(rpois(36, 9), 6L)
  y3 <- matrix(rpois(25, 5), 5L)
  expect_false(isTRUE(hicpotts_fingerprints_agree(
    hicpotts_fingerprint(y1), hicpotts_fingerprint(y2))))
  v <- hicpotts_fingerprints_agree(hicpotts_fingerprint(y1),
                                   hicpotts_fingerprint(y3))
  expect_type(v, "character")
  expect_match(paste(v, collapse = " "), "lattice dimensions differ")
})

test_that("differing covariates are detected", {
  set.seed(5); N <- 5L
  y <- matrix(rpois(N * N, 4), N)
  xa <- list(GC = list(matrix(0.1, N, N)))
  xb <- list(GC = list(matrix(0.9, N, N)))
  v <- hicpotts_fingerprints_agree(hicpotts_fingerprint(y, xa),
                                   hicpotts_fingerprint(y, xb))
  expect_match(paste(v, collapse = " "), "covariate 'GC' differs")
})

test_that("pooling chains from different data is refused", {
  set.seed(6); N <- 4L
  probs <- list(component1 = matrix(0.6, N, N),
                component2 = matrix(0.3, N, N),
                component3 = matrix(0.1, N, N))
  mk <- function(y) list(
    chains = lapply(1:3, function(k) matrix(k, 6L, 5L)),
    size = matrix(2, 3L, 6L), z_final = matrix(1L, N, N),
    z_checkpoints = list(), z_probabilities = probs,
    z_probability_draws = 300L,
    data_fingerprint = hicpotts_fingerprint(y))
  f1 <- mk(matrix(rpois(N * N, 5), N))
  f2 <- mk(matrix(rpois(N * N, 50), N))
  expect_error(classify_hicpotts(list(f1, f2), relabel = FALSE),
               "Refusing classification pooling")
  ## Same data pools without complaint.
  y <- matrix(rpois(N * N, 5), N)
  expect_s3_class(classify_hicpotts(list(mk(y), mk(y)), relabel = FALSE),
                  "hicpotts_classification")
})

test_that("pooling different gamma transitions is refused", {
  N <- 3L
  probs <- list(component1 = matrix(0.6, N, N),
                component2 = matrix(0.3, N, N),
                component3 = matrix(0.1, N, N))
  mk <- function(method) {
    gamma <- structure(rep(0.3, 6L), gamma_method = method)
    list(chains = lapply(1:3, function(k) matrix(k, 6L, 5L)),
         gamma = gamma, size = matrix(2, 3L, 6L),
         z_final = matrix(1L, N, N), z_checkpoints = list(),
         z_probabilities = probs, z_probability_draws = 300L,
         data_fingerprint = hicpotts_fingerprint(matrix(1, N, N)))
  }
  expect_error(classify_hicpotts(list(mk("exchange"), mk("abc")),
                                 relabel = FALSE),
               "different gamma transitions")
})

test_that("shuffled classification data is rejected", {
  N <- 4L
  probs <- list(component1 = matrix(0.6, N, N),
                component2 = matrix(0.3, N, N),
                component3 = matrix(0.1, N, N))
  fit <- list(chains = lapply(1:3, function(k) matrix(k, 6L, 5L)),
              size = matrix(2, 3L, 6L), z_final = matrix(1L, N, N),
              z_checkpoints = list(), z_probabilities = probs,
              z_probability_draws = 300L)
  idx <- arrayInd(seq_len(N * N), .dim = c(N, N))
  ok <- data.frame(row = idx[, 1L], column = idx[, 2L])
  expect_s3_class(classify_hicpotts(fit, data = ok, relabel = FALSE),
                  "hicpotts_classification")
  set.seed(8)
  bad <- ok[sample.int(nrow(ok)), ]
  expect_error(classify_hicpotts(fit, data = bad, relabel = FALSE),
               "not in the stored column-major lattice order")
})

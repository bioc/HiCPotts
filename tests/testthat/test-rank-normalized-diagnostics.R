## Item 10: rank-normalized split R-hat, bulk ESS and tail ESS replace the
## normality-assuming diagnostics; posterior expected occupancy is reported.

test_that("rank-normalized R-hat is essentially invariant to monotone warping", {
  set.seed(21)
  ch <- lapply(1:4, function(i) rnorm(600))
  plain <- .hicpotts_rhat(ch)
  warped <- .hicpotts_rhat(lapply(ch, exp))
  ## The rank-normalized BULK term is exactly invariant. The folded tail term
  ## is not -- folding about the median is a function of the values, not their
  ## ranks -- so the reported maximum of the two agrees only to within
  ## Monte-Carlo noise. That is the intended behaviour: exact invariance where
  ## it is available, without giving up tail sensitivity.
  expect_equal(plain, warped, tolerance = 1e-3)
})

test_that("converged chains pass and disagreeing chains fail", {
  set.seed(22)
  good <- lapply(1:4, function(i) rnorm(800))
  expect_lt(.hicpotts_rhat(good), 1.05)
  bad <- lapply(1:4, function(i) rnorm(800, mean = i * 5))
  expect_gt(.hicpotts_rhat(bad), 1.5)
})

test_that("a frozen chain is non-diagnostic rather than perfect", {
  frozen <- lapply(1:4, function(i) rep(2.5, 500))
  expect_true(is.na(.hicpotts_rhat(frozen)))
  expect_true(is.na(.hicpotts_bulk_ess(frozen)))
  expect_true(is.na(.hicpotts_tail_ess(frozen)))
})

test_that("bulk and tail ESS are finite, positive and bounded by the draws", {
  set.seed(23)
  n_per_chain <- 800L; n_chains <- 4L
  ch <- lapply(seq_len(n_chains), function(i)
    as.numeric(stats::filter(rnorm(n_per_chain), 0.6, "recursive")))
  b <- .hicpotts_bulk_ess(ch); tl <- .hicpotts_tail_ess(ch)
  ## Splitting each chain in half preserves the total number of draws.
  total <- n_per_chain * n_chains
  expect_true(is.finite(b) && b > 0 && b <= total)
  expect_true(is.finite(tl) && tl > 0 && tl <= total)
})

test_that("autocorrelated chains have lower bulk ESS than independent ones", {
  set.seed(24)
  indep <- lapply(1:4, function(i) rnorm(1000))
  corr  <- lapply(1:4, function(i) as.numeric(stats::filter(rnorm(1000), 0.9, "recursive")))
  expect_gt(.hicpotts_bulk_ess(indep), .hicpotts_bulk_ess(corr))
})

test_that("tail ESS also degrades under autocorrelation", {
  set.seed(25)
  indep <- lapply(1:4, function(i) rnorm(1000))
  corr  <- lapply(1:4, function(i) as.numeric(stats::filter(rnorm(1000), 0.9, "recursive")))
  expect_gt(.hicpotts_tail_ess(indep), .hicpotts_tail_ess(corr))
})

test_that("bulk and tail ESS are reported separately and can disagree", {
  set.seed(26)
  ch <- lapply(1:4, function(i) as.numeric(stats::filter(rnorm(1200), 0.7, "recursive")))
  b <- .hicpotts_bulk_ess(ch); tl <- .hicpotts_tail_ess(ch)
  ## Both are computed independently; the point of reporting both is that one
  ## can be adequate while the other is not.
  expect_true(is.finite(b) && is.finite(tl))
  expect_false(isTRUE(all.equal(b, tl, tolerance = 1e-6)))
})

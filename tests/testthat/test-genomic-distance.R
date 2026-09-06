test_that("genomic distance is zero on the diagonal and symmetric (H2)", {
  res <- 100000
  bs <- c(1, 100001, 200001)
  be <- bs + res - 1
  g <- expand.grid(i = seq_along(bs), j = seq_along(bs))
  d <- data.frame(
    start      = bs[g$i],
    `end.i.`   = be[g$i],
    `start.j.` = bs[g$j],
    end        = be[g$j],
    check.names = FALSE
  )

  dist <- suppressMessages(HiCPotts:::.hicpotts_genomic_distance(d))
  m <- matrix(dist, nrow = length(bs), ncol = length(bs))

  # A bin pair with itself is at distance zero.
  expect_equal(diag(m), rep(0, length(bs)))

  # The two mirrored copies of one contact must have the same distance.
  expect_equal(m, t(m))

  # Adjacent bins are exactly one bin width apart.
  expect_equal(m[1, 2], res)

  # Guard the specific defect: abs(end - start) gives resolution - 1 on the
  # diagonal and is asymmetric off it.
  legacy <- abs(d$end - d$start)
  expect_false(isTRUE(all.equal(legacy[1], 0)))
  expect_false(isTRUE(all.equal(matrix(legacy, 3, 3), t(matrix(legacy, 3, 3)))))
})

test_that("legacy start/end-only designs still work and are reported", {
  d <- data.frame(start = c(1, 1, 100001), end = c(100000, 200000, 200000))

  expect_message(
    HiCPotts:::.hicpotts_genomic_distance(d),
    "abs\\(end - start\\)",
    fixed = FALSE
  )
  expect_equal(
    suppressMessages(HiCPotts:::.hicpotts_genomic_distance(d)),
    abs(d$end - d$start)
  )
})

test_that("start.j. is preferred when present and the choice is reported", {
  d <- data.frame(
    start      = c(1, 1),
    `start.j.` = c(1, 100001),
    end        = c(100000, 200000),
    check.names = FALSE
  )

  expect_message(
    HiCPotts:::.hicpotts_genomic_distance(d),
    "abs\\(start.j. - start\\)",
    fixed = FALSE
  )
  expect_equal(
    suppressMessages(HiCPotts:::.hicpotts_genomic_distance(d)),
    c(0, 100000)
  )
})

test_that("an unusable start.j. falls back rather than propagating NA", {
  d <- data.frame(
    start      = c(1, 1),
    `start.j.` = c(NA_real_, 100001),
    end        = c(100000, 200000),
    check.names = FALSE
  )

  expect_equal(
    suppressMessages(HiCPotts:::.hicpotts_genomic_distance(d)),
    abs(d$end - d$start)
  )
})

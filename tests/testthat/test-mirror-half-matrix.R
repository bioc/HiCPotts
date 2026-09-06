half_matrix_fixture <- function(n_bins = 6L, res = 1000L, diag = TRUE) {
  bins <- seq(1L, by = res, length.out = n_bins)
  ut <- which(upper.tri(matrix(0, n_bins, n_bins), diag = diag),
              arr.ind = TRUE)
  set.seed(99)
  data.frame(
    start = bins[ut[, "row"]], start.j. = bins[ut[, "col"]],
    end.i. = bins[ut[, "row"]] + res - 1L,
    end = bins[ut[, "col"]] + res - 1L,
    chrom = "chr1",
    GC = runif(nrow(ut), 0.3, 0.7), ACC = runif(nrow(ut), 0, 1),
    TES = rpois(nrow(ut), 2), interactions = rpois(nrow(ut), 5)
  )
}

test_that("mirror turns a half matrix into the full lattice", {
  n <- 6L
  half <- half_matrix_fixture(n)
  expect_equal(nrow(half), n * (n + 1L) / 2L)

  processed <- process_data(half, N = n, standardization_y = FALSE,
                            mirror = TRUE)
  expect_equal(dim(processed$y[[1]]), c(n, n))

  # reflection must produce a symmetric count matrix and invent nothing
  ym <- processed$y[[1]]
  expect_equal(ym, t(ym))
  expect_equal(sort(as.vector(ym)[as.vector(upper.tri(ym, diag = TRUE))]),
               sort(half$interactions))
})

test_that("mirror rejects a full matrix and duplicate pairs", {
  n <- 4L
  half <- half_matrix_fixture(n)
  full <- process_data(half, N = n, standardization_y = FALSE,
                       mirror = TRUE)
  # feeding an already-complete matrix back through mirror must error
  grid <- expand.grid(start = sort(unique(c(half$start, half$start.j.))),
                      start.j. = sort(unique(c(half$start, half$start.j.))))
  rebuilt <- merge(grid, half, all.x = TRUE)
  rebuilt <- rebuilt[stats::complete.cases(rebuilt), ]
  expect_error(
    process_data(rbind(half, half), N = n, mirror = TRUE),
    "duplicate bin pairs")
})

test_that("a triangular row count is diagnosed in the error message", {
  half <- half_matrix_fixture(6L)
  # 21 rows over 6 bins is not a multiple of any N^2
  expect_error(process_data(half, N = 4L),
               "looks like a half matrix")
  expect_error(process_data(half, N = 4L), "mirror = TRUE")
})

test_that("mirror = FALSE leaves existing behaviour untouched", {
  n <- 4L
  bins <- seq(1L, by = 1000L, length.out = n)
  grid <- expand.grid(start = bins, start.j. = bins)
  set.seed(5)
  full <- data.frame(
    start = grid$start, start.j. = grid$start.j.,
    end.i. = grid$start + 999L, end = grid$start.j. + 999L,
    chrom = "chr1", GC = runif(n * n, 0.3, 0.7), ACC = runif(n * n),
    TES = rpois(n * n, 2), interactions = rpois(n * n, 4))
  a <- process_data(full, N = n, standardization_y = FALSE)
  b <- process_data(full, N = n, standardization_y = FALSE, mirror = FALSE)
  expect_equal(a, b)
})

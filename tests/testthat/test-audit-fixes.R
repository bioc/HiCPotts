cool_file <- function() {
  system.file("extdata", "BG3_WT_merged_hic_matrix_chr4_100Kb.cool",
              package = "HiCPotts")
}

test_that("C1: cooler counts survive at native resolution under any scipen", {
  skip_if_not(nzchar(cool_file()) && file.exists(cool_file()))
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("rtracklayer")

  totals <- vapply(c(0L, 999L, -5L), function(sp) {
    old <- options(scipen = sp)
    on.exit(options(old), add = TRUE)
    sum(get_data(cool_file(), "chr4", 1, 400000, 100000)$interactions)
  }, numeric(1))

  # Identical across print settings, and emphatically not zero: the defect
  # returned an all-zero matrix at default scipen.
  expect_equal(length(unique(totals)), 1L)
  expect_gt(totals[1], 0)
})

test_that("C1: native and rebinned resolutions agree on the total", {
  skip_if_not(nzchar(cool_file()) && file.exists(cool_file()))
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("rtracklayer")

  native <- sum(get_data(cool_file(), "chr4", 1, 400000, 100000)$interactions)
  rebin  <- sum(get_data(cool_file(), "chr4", 1, 400000, 200000)$interactions)
  expect_equal(native, rebin)
})

test_that("H1: a non-bin-aligned window returns aligned bins with data", {
  skip_if_not(nzchar(cool_file()) && file.exists(cool_file()))
  skip_if_not_installed("rhdf5")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("rtracklayer")

  d <- get_data(cool_file(), "chr4", 50001, 450000, 200000)
  # Previously produced user-start-anchored bins holding zero contacts.
  expect_gt(sum(d$interactions), 0)
  # Bins sit on the shared grid: start = k * resolution + 1.
  expect_true(all((d$start - 1) %% 200000 == 0))
})

test_that("H6: a frozen chain is non-diagnostic, not perfectly converged", {
  constant <- rep(1, 1000)
  expect_true(is.na(HiCPotts:::.hicpotts_ess(constant)))
  expect_true(is.na(HiCPotts:::.hicpotts_rhat(replicate(4, constant, simplify = FALSE))))

  # A genuinely moving chain still gets a finite, positive ESS.
  set.seed(1)
  moving <- rnorm(1000)
  expect_true(is.finite(HiCPotts:::.hicpotts_ess(moving)))
  expect_gt(HiCPotts:::.hicpotts_ess(moving), 0)
})

test_that("H5: components 2 and 3 default to the fitted (non-inflated) family", {
  expect_false(formals(compute_HMRFHiC_probabilities)$consistent_dist)
})

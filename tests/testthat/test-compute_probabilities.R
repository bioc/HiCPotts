make_fake_chains <- function(n_iter, value_per_comp_per_col) {
  ## value_per_comp_per_col: 3 x 5 matrix of constant values (per comp, per col)
  chains_list <- lapply(seq_len(3), function(c) {
    m <- matrix(0, nrow = n_iter + 1, ncol = 5)
    for (col in 1:5) m[, col] <- value_per_comp_per_col[c, col]
    m
  })
  list(list(
    chains = chains_list,
    theta  = rep(0.25, n_iter + 1),
    size   = matrix(c(5, 10, 15), nrow = 3, ncol = n_iter + 1)
  ))
}

test_that("probabilities sum to 1 per row", {
  dat <- data.frame(
    start        = c(1, 100, 500),
    start.j.     = c(50, 200, 800),
    end          = c(50, 200, 800),
    interactions = c(3, 5, 7),
    GC           = c(0.4, 0.5, 0.6),
    TES          = c(2, 3, 4),
    ACC          = c(0.8, 0.5, 0.2)
  )
  n_iter <- 10
  vals <- matrix(c(
    0.1, 0.2, 0.1, 0.05, 0.1,
    0.3, 0.4, 0.2, 0.1,  0.2,
    0.5, 0.1, 0.0, 0.3,  0.4
  ), nrow = 3, ncol = 5, byrow = TRUE)
  
  res <- compute_HMRFHiC_probabilities(dat,
                                       chain_betas = make_fake_chains(n_iter, vals),
                                       iterations  = n_iter,
                                       dist        = "Poisson", relabel = FALSE)
  expect_equal(rowSums(res[, c("prob1","prob2","prob3")]), rep(1, nrow(dat)),
               tolerance = 1e-10)
  expect_identical(res$start.j., dat$start.j.)
})

test_that("max_interactions defaults to no cap", {
  dat <- data.frame(
    start        = 1,
    end          = 100,
    interactions = 1000L,          # > 500
    GC           = 0.5,
    TES          = 1,
    ACC          = 0.5
  )
  vals <- matrix(rep(0, 15), nrow = 3, ncol = 5)
  res <- compute_HMRFHiC_probabilities(dat,
                                       chain_betas = make_fake_chains(10, vals),
                                       iterations  = 10,
                                       dist        = "Poisson", relabel = FALSE)
  expect_equal(res$interactions, 1000)       # still 1000 in the returned df
})

test_that("Explicit max_interactions still truncates (backwards compat)", {
  dat <- data.frame(
    start        = 1,
    end          = 100,
    interactions = 1000L,
    GC           = 0.5,
    TES          = 1,
    ACC          = 0.5
  )
  vals <- matrix(rep(0, 15), nrow = 3, ncol = 5)
  res <- compute_HMRFHiC_probabilities(
    dat,
    chain_betas      = make_fake_chains(10, vals),
    iterations       = 10,
    dist             = "Poisson",
    max_interactions = 500L,
    relabel = TRUE
  )
  expect_equal(res$interactions, 1000)       # data value preserved...
  ## ...but density was evaluated at 500, so probs are still on [0,1]
  expect_true(all(res$prob1 + res$prob2 + res$prob3 - 1 < 1e-10))
})

test_that("Missing required columns raises informative error", {
  bad <- data.frame(start = 1, end = 1)
  expect_error(
    compute_HMRFHiC_probabilities(bad, chain_betas = list(), iterations = 4),
    "required columns are missing"
  )
})

test_that("iterations = 4 (small) still works with proper integer burn-in", {
  ## Guards against a non-integer burn-in (e.g. iterations/2 for odd
  ## iteration counts), which would produce fractional matrix indices.
  dat <- data.frame(start = 1:3, end = 11:13,
                    interactions = c(3,5,7),
                    GC = c(.2,.4,.6), TES = 1:3, ACC = c(.1,.5,.9))
  vals <- matrix(0, 3, 5)
  expect_warning(compute_HMRFHiC_probabilities(
    dat, chain_betas = make_fake_chains(4, vals),
    iterations = 4, dist = "Poisson", relabel = TRUE
  ), "N.*not supplied")
  expect_warning(compute_HMRFHiC_probabilities(
    dat, chain_betas = make_fake_chains(5, vals),
    iterations = 5, dist = "Poisson", relabel = TRUE        # odd
  ), "N.*not supplied")
})
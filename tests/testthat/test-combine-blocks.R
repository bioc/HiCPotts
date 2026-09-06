make_block_chain <- function(component, N = 2L) {
    probabilities <- lapply(seq_len(3L), function(k) {
        matrix(if (k == component) 0.9 else 0.05, N, N)
    })
    batches <- array(0, dim = c(N, N, 3L, 2L))
    for (k in seq_len(3L)) {
        batches[, , k, 1L] <- probabilities[[k]]
        batches[, , k, 2L] <- probabilities[[k]]
    }
    list(
        chains = lapply(seq_len(3L), function(k) {
            matrix(k, nrow = 4L, ncol = 5L)
        }),
        z_final = matrix(component, N, N),
        z_probabilities = probabilities,
        z_probability_draws = 10L,
        z_probability_batches = list(
            batch_means = batches, batch_draws = c(5, 5)
        )
    )
}

make_block_robust <- function(selected, N = 2L) {
    all_fits <- lapply(seq_len(3L), make_block_chain, N = N)
    structure(list(
        fits = all_fits[selected], all_fits = all_fits,
        mode_selection = list(selected = selected, multimodal = FALSE),
        diagnostics = list(
            parameters = data.frame(
                parameter = "gamma", estimate = selected[1L] / 10,
                bulk_ESS = 500, tail_ESS = 500, Rhat = 1
            ),
            reliability_flags = data.frame(
                criterion = c("effective_sample_size", "split_rhat"),
                passed = TRUE, threshold = c(">= 200", "<= 1.01")
            ),
            warnings = character()
        ),
        settings = list(N = N, iterations = 4L)
    ), class = "hicpotts_robust_fit")
}

make_processed_blocks <- function(n_blocks = 2L, N = 2L) {
    make_values <- function(offset) matrix(seq_len(N * N) + offset, N, N)
    list(
        x_vars = setNames(lapply(seq_len(4L), function(k) {
            lapply(seq_len(n_blocks), function(j) make_values(k + j))
        }), c("distance", "GC", "TES", "ACC")),
        y = lapply(seq_len(n_blocks), function(j) make_values(10L + j))
    )
}

test_that("one block needs no combined container", {
    fit <- make_block_robust(c(1L, 3L))
    expect_identical(combine_hicpotts_blocks(fit), fit)
    expect_identical(combine_hicpotts_blocks(list(fit)), fit)
})

test_that("each block keeps its own mode-selected chains", {
    fits <- list(left = make_block_robust(c(1L, 3L)),
                 right = make_block_robust(2L))
    combined <- combine_hicpotts_blocks(fits, make_processed_blocks())

    expect_s3_class(combined, "hicpotts_block_fit")
    expect_equal(combined$block_table$selected_chains, c("1,3", "2"))
    expect_length(combined$blocks[[1L]]$fits, 2L)
    expect_length(combined$blocks[[2L]]$fits, 1L)

    result <- classify_hicpotts(
        combined, relabel = FALSE, min_draws = 1L, reflect = "never"
    )
    expect_s3_class(result, "hicpotts_classification")
    expect_equal(as.integer(table(result$block)), c(4L, 4L))
    expect_equal(names(table(result$block)), c("left", "right"))
    expect_true(all(result$map_component[result$block == "left"] == 1L))
    expect_true(all(result$map_component[result$block == "right"] == 2L))
    expect_equal(attr(result, "chains_by_block")$selected_chains,
                 c("1,3", "2"))
    expect_identical(result$global_cell_index, seq_len(8L))
})

test_that("whole-map data are split using fitted block sizes", {
    fits <- list(make_block_robust(1L), make_block_robust(2L))
    original <- data.frame(id = seq_len(8L))
    combined <- combine_hicpotts_blocks(fits, data = original)
    expect_equal(combined$data[[1L]]$id, 1:4)
    expect_equal(combined$data[[2L]]$id, 5:8)
})

test_that("processed matrices cannot replace original probability data", {
    fits <- list(make_block_robust(1L), make_block_robust(2L))
    combined <- combine_hicpotts_blocks(fits, make_processed_blocks())
    expect_error(
        compute_HMRFHiC_probabilities(chain_betas = combined),
        "Original per-cell data are required"
    )
})

test_that("block allocation diagnostics are not pooled", {
    fits <- list(A = make_block_robust(c(1L, 3L)),
                 B = make_block_robust(2L))
    combined <- combine_hicpotts_blocks(fits, make_processed_blocks())
    diagnostics <- allocation_diagnostics(combined)
    expect_named(diagnostics$blocks, c("A", "B"))
    expect_equal(diagnostics$summary$selected_chains, c("1,3", "2"))
    expect_equal(nrow(diagnostics$summary), 2L)
})

test_that("parameter summaries remain block-specific", {
    fits <- list(A = make_block_robust(1L), B = make_block_robust(2L))
    combined <- combine_hicpotts_blocks(fits, make_processed_blocks())
    summary <- summarise_hicpotts_parameters(
        combined, require_reliable = FALSE, diagnose_covariates = FALSE
    )
    expect_s3_class(summary, "hicpotts_block_parameter_summary")
    expect_equal(summary$block, c("A", "B"))
    expect_equal(summary$estimate, c(0.1, 0.2))
    expect_equal(summary$selected_chains, c("1", "2"))
})

make_pooling_chain <- function(
    N, betas, gamma, theta, size, expected_component_cells
) {
    cells <- N * N
    stopifnot(
        identical(dim(betas), c(3L, 5L)), length(size) == 3L,
        length(expected_component_cells) == 3L,
        isTRUE(all.equal(sum(expected_component_cells), cells))
    )
    iterations <- 10L
    probabilities <- lapply(expected_component_cells, function(value) {
        matrix(value / cells, N, N)
    })
    size_draws <- do.call(rbind, lapply(size, function(value) {
        rep(value, iterations)
    }))
    list(
        chains = lapply(seq_len(3L), function(component) {
            matrix(
                rep(betas[component, ], each = iterations),
                nrow = iterations, ncol = 5L
            )
        }),
        gamma = rep(gamma, iterations),
        theta = rep(theta, iterations),
        size = size_draws,
        z_final = matrix(which.max(expected_component_cells), N, N),
        z_probabilities = probabilities,
        z_probability_draws = iterations,
        sampler_settings = list(distribution = "ZINB")
    )
}

make_pooling_robust <- function(chain) {
    sets <- .hicpotts_draw_sets(list(chain), burnin = 0L)
    parameters <- do.call(rbind, lapply(names(sets), function(parameter) {
        draws <- unlist(sets[[parameter]], use.names = FALSE)
        data.frame(
            parameter = parameter, estimate = mean(draws),
            posterior_sd = stats::sd(draws),
            CI_lower = min(draws), CI_upper = max(draws),
            CI_width = max(draws) - min(draws),
            probability_positive = mean(draws > 0),
            probability_negative = mean(draws < 0),
            ESS = 500, bulk_ESS = 500, tail_ESS = 500, Rhat = 1,
            stringsAsFactors = FALSE
        )
    }))
    structure(list(
        fits = list(chain), all_fits = list(chain),
        mode_selection = list(selected = 1L, multimodal = FALSE),
        diagnostics = list(
            parameters = parameters,
            reliability_flags = data.frame(
                criterion = c("effective_sample_size", "split_rhat"),
                passed = TRUE, threshold = c(">= 200", "<= 1.01"),
                stringsAsFactors = FALSE
            ),
            warnings = character()
        ),
        settings = list(N = nrow(chain$z_final), burnin = 0L)
    ), class = "hicpotts_robust_fit")
}

test_that("posterior block pooling uses the parameter-appropriate weights", {
    beta_a <- matrix(1, 3L, 5L)
    beta_b <- matrix(5, 3L, 5L)
    block_a <- make_pooling_robust(make_pooling_chain(
        N = 2L, betas = beta_a, gamma = 0.2, theta = 0.1,
        size = c(2, 4, 6), expected_component_cells = c(1, 3, 0)
    ))
    block_b <- make_pooling_robust(make_pooling_chain(
        N = 3L, betas = beta_b, gamma = 0.8, theta = 0.4,
        size = c(8, 10, 12), expected_component_cells = c(4, 1, 4)
    ))
    combined <- combine_hicpotts_blocks(list(A = block_a, B = block_b))

    result <- summarise_hicpotts_parameters(
        combined, require_reliable = FALSE,
        diagnose_covariates = FALSE,
        pool_blocks = "posterior_weighted", pooling_draws = 2000L
    )

    expect_s3_class(result, "hicpotts_pooled_parameter_summary")
    component2 <- result[result$parameter == "component2:distance", ]
    expect_equal(component2$estimate, (3 * 1 + 1 * 5) / 4)
    expect_equal(
        component2$weighting_method,
        "posterior_expected_component_cells"
    )
    expect_equal(component2$effective_cells, 4)

    gamma <- result[result$parameter == "gamma", ]
    expect_equal(gamma$estimate, (4 * 0.2 + 9 * 0.8) / 13)
    expect_equal(gamma$weighting_method, "analysed_cells")
    expect_equal(gamma$effective_cells, 13)

    component3 <- result[result$parameter == "component3:intercept", ]
    expect_equal(component3$estimate, 5)
    expect_equal(component3$n_blocks_used, 1L)
    expect_equal(component3$blocks_excluded, "A")
    expect_equal(component3$CI_lower, 5)
    expect_equal(component3$CI_upper, 5)

    expect_equal(
        attr(result, "posterior_expected_component_cells"),
        rbind(A = c(1, 3, 0), B = c(4, 1, 4)),
        ignore_attr = TRUE
    )
    expect_s3_class(
        attr(result, "block_summaries"),
        "hicpotts_block_parameter_summary"
    )
})

test_that("minimum expected occupancy excludes weak component blocks", {
    beta_a <- matrix(1, 3L, 5L)
    beta_b <- matrix(5, 3L, 5L)
    combined <- combine_hicpotts_blocks(list(
        A = make_pooling_robust(make_pooling_chain(
            2L, beta_a, 0.2, 0.1, c(2, 4, 6), c(1, 3, 0)
        )),
        B = make_pooling_robust(make_pooling_chain(
            3L, beta_b, 0.8, 0.4, c(8, 10, 12), c(4, 1, 4)
        ))
    ))
    result <- summarise_hicpotts_parameters(
        combined, require_reliable = FALSE,
        diagnose_covariates = FALSE,
        pool_blocks = "posterior_weighted", pooling_draws = 1000L,
        pooling_min_expected_cells = 2
    )
    component2 <- result[result$parameter == "component2:distance", ]
    expect_equal(component2$estimate, 1)
    expect_equal(component2$n_blocks_used, 1L)
    expect_equal(component2$blocks_excluded, "B")
})

test_that("one fitted matrix retains the ordinary summary behavior", {
    beta <- matrix(2, 3L, 5L)
    fit <- make_pooling_robust(make_pooling_chain(
        2L, beta, 0.3, 0.1, c(2, 4, 6), c(1, 2, 1)
    ))
    one <- combine_hicpotts_blocks(list(only = fit))
    expect_identical(one, fit)
    result <- summarise_hicpotts_parameters(
        one, require_reliable = FALSE, diagnose_covariates = FALSE,
        pool_blocks = "posterior_weighted"
    )
    expect_s3_class(result, "hicpotts_parameter_summary")
    expect_false(inherits(result, "hicpotts_pooled_parameter_summary"))
})

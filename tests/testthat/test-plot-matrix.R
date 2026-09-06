test_that("dual-triangle heatmap keeps one reversed y scale", {
    skip_if_not_installed("ggnewscale")

    results <- data.frame(
        start = c(1e6, 1e6, 2e6, 2e6, 3e6),
        end = c(1e6, 2e6, 2e6, 3e6, 3e6),
        prob2 = c(0.10, 0.80, 0.20, 0.70, 0.30),
        interactions = c(5, 20, 8, 15, 6)
    )

    expect_silent(
        plot <- plot_upper_prob_lower_count(results)
    )
    expect_s3_class(plot, "ggplot")
    expect_identical(plot$scales$get_scales("y")$trans$name, "reverse")
})

test_that("symmetric plotting rejects inconsistent mirrored counts", {
    skip_if_not_installed("ggnewscale")

    results <- data.frame(
        start = c(1, 2, 1, 2),
        end = c(1, 1, 2, 2),
        prob2 = c(0.1, 0.2, 0.8, 0.3),
        interactions = c(5, 20, 10, 6)
    )

    expect_error(
        plot_upper_prob_lower_count(results, symmetric_matrix = TRUE),
        "mirrored count values differ"
    )
})

test_that("ordered plotting preserves different upper and lower cells", {
    skip_if_not_installed("ggnewscale")

    results <- data.frame(
        start = c(1, 2, 1, 2),
        end = c(1, 1, 2, 2),
        prob2 = c(0.1, 0.2, 0.8, 0.3),
        interactions = c(5, 20, 10, 6)
    )
    plot <- plot_upper_prob_lower_count(
        results, use_log_count = FALSE, symmetric_matrix = FALSE
    )
    upper <- plot$layers[[1L]]$data
    lower <- plot$layers[[2L]]$data

    expect_equal(
        upper$prob[upper$bin1 == 2 & upper$bin2 == 1],
        0.2
    )
    expect_equal(
        lower$count_plot[lower$bin1 == 1 & lower$bin2 == 2],
        10
    )
})

test_that("symmetric_matrix must be one logical value", {
    skip_if_not_installed("ggnewscale")

    results <- data.frame(
        start = 1, end = 1, prob2 = 0.5, interactions = 2
    )
    expect_error(
        plot_upper_prob_lower_count(results, symmetric_matrix = "yes"),
        "must be TRUE or FALSE"
    )
})

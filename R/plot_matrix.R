#' @title Plot a dual-triangle Hi-C heatmap
#'
#' @description Creates a Hi-C style heatmap where the visually upper triangle
#' shows posterior probabilities and the visually lower triangle shows observed
#' interaction counts. Symmetric and ordered input matrices are both supported.
#' This is designed for HiCPotts-style long-format output where each row
#' represents an interaction pair. By default it uses columns commonly found in
#' package outputs: \code{start}, \code{end}, \code{prob2}, and
#' \code{interactions}.
#'
#' @param results A data.frame containing interaction results.
#' @param bin1_col Column name for the first bin coordinate.
#' @param bin2_col Column name for the second bin coordinate.
#' @param prob_col Column name for the posterior probability to display in the
#'   upper triangle.
#' @param count_col Column name for the observed interaction counts to display
#'   in the lower triangle.
#' @param chr_label Chromosome label for the axes.
#' @param title Plot title.
#' @param prob_agg Aggregation function for duplicate probability rows. When
#'   \code{symmetric_matrix = TRUE}, it also combines the two orientations of
#'   an unordered pair. Defaults to \code{max}.
#' @param count_agg Aggregation function for duplicate count rows. When
#'   \code{symmetric_matrix = TRUE}, it also combines the two orientations of
#'   an unordered pair. Defaults to \code{max}.
#' @param use_log_count Logical; if \code{TRUE}, plot \code{log1p(count)} in
#'   the lower triangle. If \code{FALSE}, plot raw counts.
#' @param symmetric_matrix Logical; if \code{TRUE} (the default), treat
#'   \eqn{(i,j)} and \eqn{(j,i)} as two orientations of the same Hi-C contact.
#'   A single complete triangle is reflected, while a full matrix is checked
#'   for mirrored counts before its duplicate orientations are aggregated. If
#'   \code{FALSE}, preserve ordered pairs and use the supplied upper and lower
#'   triangles without reflection or cross-triangle aggregation.
#'
#' @return A ggplot object.
#'
#' @examples
#' sim_res <- data.frame(
#'     start = c(1e6, 1e6, 2e6, 2e6, 3e6),
#'     end = c(1e6, 2e6, 2e6, 3e6, 3e6),
#'     prob2 = c(0.10, 0.80, 0.20, 0.70, 0.30),
#'     interactions = c(5, 20, 8, 15, 6)
#' )
#'
#' p <- plot_upper_prob_lower_count(sim_res)
#' print(p)
#'
#' @importFrom rlang .data
#'
#' @export
#'


plot_upper_prob_lower_count <- function(
    results,
    bin1_col = "start",
    bin2_col = "end",
    prob_col = "prob2",
    count_col = "interactions",
    chr_label = "2L",
    title = "Significant interactions detected by HiCPotts",
    prob_agg = max,
    count_agg = max,
    use_log_count = TRUE,
    symmetric_matrix = TRUE
) {
    if (!is.data.frame(results)) {
        stop("'results' must be a data.frame.")
    }

    if (!requireNamespace("ggnewscale", quietly = TRUE)) {
        stop(
            "Package 'ggnewscale' is required for ",
            "plot_upper_prob_lower_count(). ",
            "Please install it with install.packages('ggnewscale')."
        )
    }


    needed <- c(bin1_col, bin2_col, prob_col, count_col)
    missing_cols <- setdiff(needed, names(results))
    if (length(missing_cols) > 0L) {
        stop(
            "Missing required columns in 'results': ",
            paste(missing_cols, collapse = ", ")
        )
    }

    df <- results[, needed, drop = FALSE]
    names(df) <- c("bin1", "bin2", "prob", "count")

    df <- df[
        is.finite(df$bin1) &
            is.finite(df$bin2) &
            is.finite(df$prob) &
            is.finite(df$count), ,
        drop = FALSE
    ]

    if (nrow(df) == 0L) {
        stop("No finite rows found in 'results' after filtering.")
    }

    if (!is.function(prob_agg)) {
        stop("'prob_agg' must be a function.")
    }
    if (!is.function(count_agg)) {
        stop("'count_agg' must be a function.")
    }
    if (!is.logical(symmetric_matrix) || length(symmetric_matrix) != 1L ||
        is.na(symmetric_matrix)) {
        stop("'symmetric_matrix' must be TRUE or FALSE.", call. = FALSE)
    }

    bins <- sort(unique(c(df$bin1, df$bin2)))
    grid <- expand.grid(
        bin1 = bins,
        bin2 = bins,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    if (isTRUE(symmetric_matrix)) {
        # A full matrix declared symmetric must actually contain mirrored
        # counts. This catches column shifts and malformed reflection before
        # they produce a misleading heatmap. A single complete triangle is
        # valid and is reflected below.
        ordered_counts <- stats::aggregate(
            count ~ bin1 + bin2, data = df, FUN = count_agg
        )
        pair_key <- function(a, b) paste(a, b, sep = "\r")
        reverse_index <- match(
            pair_key(ordered_counts$bin2, ordered_counts$bin1),
            pair_key(ordered_counts$bin1, ordered_counts$bin2)
        )
        if (all(!is.na(reverse_index))) {
            reflected_count <- ordered_counts$count[reverse_index]
            tolerance <- sqrt(.Machine$double.eps) * pmax(
                1, abs(ordered_counts$count), abs(reflected_count)
            )
            mismatch <- abs(ordered_counts$count - reflected_count) > tolerance
            if (any(mismatch)) {
                stop(
                    "'symmetric_matrix = TRUE' but mirrored count values ",
                    "differ. Correct the input matrix or set ",
                    "'symmetric_matrix = FALSE' for ordered data.",
                    call. = FALSE
                )
            }
        }

        # Treat interaction pairs as unordered and build a reflected grid.
        df$key1 <- pmin(df$bin1, df$bin2)
        df$key2 <- pmax(df$bin1, df$bin2)
        prob_df <- stats::aggregate(
            prob ~ key1 + key2, data = df, FUN = prob_agg
        )
        count_df <- stats::aggregate(
            count ~ key1 + key2, data = df, FUN = count_agg
        )
        pair_df <- merge(
            prob_df, count_df, by = c("key1", "key2"), all = TRUE
        )
        grid$key1 <- pmin(grid$bin1, grid$bin2)
        grid$key2 <- pmax(grid$bin1, grid$bin2)
        plot_df <- merge(
            grid, pair_df, by = c("key1", "key2"), all.x = TRUE,
            sort = FALSE
        )
    } else {
        # Preserve direction: do not combine (i,j) with (j,i).
        prob_df <- stats::aggregate(
            prob ~ bin1 + bin2, data = df, FUN = prob_agg
        )
        count_df <- stats::aggregate(
            count ~ bin1 + bin2, data = df, FUN = count_agg
        )
        pair_df <- merge(
            prob_df, count_df, by = c("bin1", "bin2"), all = TRUE
        )
        plot_df <- merge(
            grid, pair_df, by = c("bin1", "bin2"), all.x = TRUE,
            sort = FALSE
        )
    }

    plot_df$prob[is.na(plot_df$prob)] <- 0
    plot_df$count[is.na(plot_df$count)] <- 0

    plot_df$x_mb <- plot_df$bin1 / 1e6
    plot_df$y_mb <- plot_df$bin2 / 1e6

    if (use_log_count) {
        plot_df$count_plot <- log1p(plot_df$count)
        count_legend <- "log1p(count)"
    } else {
        plot_df$count_plot <- plot_df$count
        count_legend <- "Count"
    }

    # With scale_y_reverse(), the visually upper triangle corresponds to bin1 >
    # bin2
    upper_df <- plot_df[plot_df$bin1 > plot_df$bin2, , drop = FALSE]
    lower_df <- plot_df[plot_df$bin1 <= plot_df$bin2, , drop = FALSE]

    axis_breaks <- pretty(range(c(plot_df$x_mb, plot_df$y_mb)), n = 6)

    ggplot2::ggplot() +
        ggplot2::geom_raster(
            data = upper_df,
            ggplot2::aes(
                x = .data[["x_mb"]],
                y = .data[["y_mb"]],
                fill = .data[["prob"]]
            )
        ) +
        ggplot2::scale_fill_gradientn(
            colours = c("#F7F7F7", "#DCE6F2", "#9FC0DE", "#4F92C6", "#0B3C78"),
            limits = c(0, 1),
            name = "Probability"
        ) +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_raster(
            data = lower_df,
            ggplot2::aes(
                x = .data[["x_mb"]],
                y = .data[["y_mb"]],
                fill = .data[["count_plot"]]
            )
        ) +
        ggplot2::scale_fill_gradientn(
            colours = c("#F7F7F7", "#E5E5E5", "#BDBDBD", "#7F7F7F", "#252525"),
            name = count_legend
        ) +
        ggplot2::coord_fixed() +
        ggplot2::scale_x_continuous(
            breaks = axis_breaks,
            labels = axis_breaks
        ) +
        ggplot2::scale_y_reverse(
            breaks = axis_breaks,
            labels = axis_breaks
        ) +
        ggplot2::labs(
            title = title,
            x = paste0(chr_label, "\nMb"),
            y = paste0(chr_label, "\nMb")
        ) +
        ggplot2::theme_bw(base_size = 16) +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
            axis.title = ggplot2::element_text(face = "bold"),
            axis.text = ggplot2::element_text(colour = "black")
        )
}

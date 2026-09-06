#' @title Extract one HiCPotts fit object
#'
#' @description
#' Internal helper that returns a single HiCPotts fit object whether the input
#' is already one fit or a list of fits.
#'
#' @param x A HiCPotts fit object, or a list of such objects.
#' @param index Integer index used when \code{x} is a list.
#'
#' @return A single HiCPotts fit object.
#'
#' @noRd
.get_hicpotts_result <- function(x, index = 1L) {
    if (!is.numeric(index) || length(index) != 1L || index < 1L) {
        stop("'index' must be a positive scalar.")
    }
    index <- as.integer(index)
    if (inherits(x, "hicpotts_robust_fit")) x <- x$fits

    is_single_result <- is.list(x) &&
        !is.null(x$chains) &&
        is.list(x$chains) &&
        length(x$chains) == 3L

    if (is_single_result) {
        return(x)
    }

    if (is.list(x) &&
        length(x) >= index &&
        is.list(x[[index]]) &&
        !is.null(x[[index]]$chains) &&
        is.list(x[[index]]$chains) &&
        length(x[[index]]$chains) == 3L) {
        return(x[[index]])
    }

    stop("Input does not look like a HiCPotts fit or list of fits.")
}


#' @title Plot HiCPotts MCMC traces by component
#'
#' @description
#' Produces trace plots for the component-specific regression coefficients
#' and, when available, the size parameters. Optionally also plots the global
#' parameters \code{gamma} and \code{theta}.
#'
#' This function is intended for raw fitted HiCPotts objects returned by the
#' MCMC fitting functions, not for the data frame returned by
#' \code{compute_HMRFHiC_probabilities()}.
#'
#' @param fit A single HiCPotts fit object or a list of such objects.
#' @param index Integer index used when \code{fit} is a list.
#' @param burnin Optional burn-in iteration. Defaults to half the iterations.
#' @param component_names Character vector of length 3 giving component labels.
#' @param beta_names Character vector giving the regression coefficient names.
#' @param plot_globals Logical; if \code{TRUE}, also plot \code{gamma} and
#'   \code{theta} when present.
#' @param plot_size Logical; if \code{TRUE}, plot the component-wise size
#'   parameters when present.
#' @param ask Logical; if \code{TRUE} in an interactive session, pause between
#'   pages.
#'
#' @return Invisibly returns the extracted HiCPotts fit object used for
#' plotting.
#'
#'
#'
#' @examples
#' set.seed(4921)
#' n_draw <- 12L
#'
#' make_chain <- function(intercept) {
#'     cbind(
#'         rnorm(n_draw, intercept, 0.05),
#'         rnorm(n_draw),
#'         rnorm(n_draw),
#'         rnorm(n_draw),
#'         rnorm(n_draw)
#'     )
#' }
#'
#' fit <- list(
#'     chains = list(
#'         make_chain(-1),
#'         make_chain(0),
#'         make_chain(1)
#'     ),
#'     size = matrix(5, nrow = 3, ncol = n_draw),
#'     theta = runif(n_draw, 0.1, 0.2),
#'     gamma = runif(n_draw, 0.5, 0.8)
#' )
#'
#' plot_hicpotts_mcmc_by_component(
#'     fit,
#'     burnin = 5,
#'     ask = FALSE
#' )
#'
#' @export
plot_hicpotts_mcmc_by_component <- function(
    fit,
    index = 1L,
    burnin = NULL,
    # Component 2 is centred on the low quantile and component 3 on the high
    # quantile of pooled signal counts (see prior_modelling.R), so component 2
    # is the true signal component and component 3 is the false-signal
    # component.
    component_names = c("Noise", "Signal", "False signal"),
    beta_names = c("Intercept", "Distance", "GC", "TES", "ACC"),
    plot_globals = TRUE,
    plot_size = TRUE,
    ask = interactive()
) {
    res <- .get_hicpotts_result(fit, index = index)

    if (length(component_names) != 3L) {
        stop("'component_names' must have length 3.")
    }
    if (length(beta_names) < 1L) {
        stop("'beta_names' must contain at least one name.")
    }

    if (!all(vapply(res$chains, is.matrix, logical(1)))) {
        stop("Each element of res$chains must be a matrix.")
    }

    n_iter <- nrow(res$chains[[1]])
    if (!all(vapply(res$chains, nrow, integer(1)) == n_iter)) {
        stop("All component chains must have the same number of rows.")
    }

    n_beta <- length(beta_names)
    if (any(vapply(res$chains, ncol, integer(1)) < n_beta)) {
        stop("Each chain matrix must have at least length(beta_names) columns.")
    }

    if (is.null(burnin)) {
        burnin <- floor(n_iter / 2)
    }
    if (!is.numeric(burnin) || length(burnin) != 1L ||
        burnin < 0 || burnin >= n_iter) {
        stop("'burnin' must be between 0 and n_iter - 1.")
    }
    burnin <- as.integer(burnin)

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    has_size <- isTRUE(plot_size) &&
        !is.null(res$size) &&
        is.matrix(res$size) &&
        nrow(res$size) == 3L &&
        ncol(res$size) == n_iter

    for (k in seq_len(3L)) {
        graphics::par(
            mfrow = c(3, 2),
            mar = c(3, 3, 3, 1),
            oma = c(0, 0, 3, 0)
        )

        for (j in seq_len(n_beta)) {
            y <- res$chains[[k]][, j]
            graphics::plot(
                y,
                type = "l",
                main = paste(component_names[k], "-", beta_names[j]),
                xlab = "Iteration",
                ylab = "Value"
            )
            graphics::abline(v = burnin, col = "red", lwd = 2, lty = 2)
        }

        if (has_size) {
            graphics::plot(
                res$size[k, ],
                type = "l",
                main = paste(component_names[k], "- Size"),
                xlab = "Iteration",
                ylab = "Value"
            )
            graphics::abline(v = burnin, col = "red", lwd = 2, lty = 2)
        } else {
            graphics::plot.new()
            graphics::text(0.5, 0.5, "No size parameter", cex = 1)
        }

        graphics::mtext(
            paste("MCMC traces for", component_names[k]),
            outer = TRUE,
            cex = 1.2
        )

        if (isTRUE(ask) && interactive() && (k < 3L || isTRUE(plot_globals))) {
            invisible(base::readline("Press <Enter> for next plot page..."))
        }
    }

    if (isTRUE(plot_globals)) {
        n_panels <- 0L
        if (!is.null(res$gamma)) n_panels <- n_panels + 1L
        if (!is.null(res$theta)) n_panels <- n_panels + 1L

        if (n_panels > 0L) {
            graphics::par(
                mfrow = c(1, max(1L, n_panels)),
                mar = c(3, 3, 3, 1),
                oma = c(0, 0, 3, 0)
            )

            if (!is.null(res$gamma)) {
                graphics::plot(
                    res$gamma,
                    type = "l",
                    main = "Gamma",
                    xlab = "Iteration",
                    ylab = "Value"
                )
                graphics::abline(v = burnin, col = "red", lwd = 2, lty = 2)
            }

            if (!is.null(res$theta)) {
                graphics::plot(
                    res$theta,
                    type = "l",
                    main = "Theta",
                    xlab = "Iteration",
                    ylab = "Value"
                )
                graphics::abline(v = burnin, col = "red", lwd = 2, lty = 2)
            }

            graphics::mtext(
                "MCMC traces for global parameters",
                outer = TRUE,
                cex = 1.2
            )
        }
    }

    invisible(res)
}


#' @title Summarise HiCPotts posterior probabilities
#'
#' @description
#' Summarises the data frame returned by \code{compute_HMRFHiC_probabilities()}.
#' For each component, this function reports the mean probability, median
#' probability, quantile-based interval across interactions, and optional
#' hard-assignment counts based on the maximum posterior probability.
#'
#' This replaces the older chain-based posterior summary table when the main
#' downstream object is the probability output from
#' \code{compute_HMRFHiC_probabilities()}.
#'
#' @param prob_result A data.frame returned by
#' \code{compute_HMRFHiC_probabilities()}.
#' @param ci_level Numeric confidence level in \code{(0, 1)} used for the
#' quantile interval across interactions.
#' @param component_names Character vector of length 3 naming the components.
#' @param include_hard_calls Logical; if \code{TRUE}, include
#' maximum-probability assignment counts and proportions.
#' @param include_interaction_summary Logical; if \code{TRUE} and the
#' \code{interactions} column is present, include mean and median interaction
#' counts among rows hard-assigned to each component.
#'
#' @return A data.frame with one row per component.
#'
#' @usage
#' summarise_hicpotts_probabilities(prob_result, ci_level = 0.95,
#' component_names = c("noise", "signal", "false signal"),
#' include_hard_calls = TRUE, include_interaction_summary = TRUE)
#'
#' @examples
#' prob_res <- data.frame(
#'     start = c(1e6, 2e6, 3e6),
#'     end = c(2e6, 3e6, 4e6),
#'     interactions = c(0, 5, 12),
#'     prob1 = c(0.80, 0.20, 0.10),
#'     prob2 = c(0.15, 0.60, 0.25),
#'     prob3 = c(0.05, 0.20, 0.65)
#' )
#'
#' summarise_hicpotts_probabilities(
#'     prob_result = prob_res, ci_level = 0.95,
#'     component_names = c("noise", "signal", "false signal"),
#'     include_hard_calls = TRUE, include_interaction_summary = TRUE
#' )
#'
#' @export
summarise_hicpotts_probabilities <- function(
    prob_result,
    ci_level = 0.95,
    # Component 2 is the true signal component and component 3 is the
    # false-signal component (see plot_hicpotts_mcmc_by_component() above).
    component_names = c("noise", "signal", "false signal"),
    include_hard_calls = TRUE,
    include_interaction_summary = TRUE
) {
    if (!is.data.frame(prob_result)) {
        stop("'prob_result' must be a data.frame.")
    }

    needed <- c("prob1", "prob2", "prob3")
    missing_cols <- setdiff(needed, names(prob_result))
    if (length(missing_cols) > 0L) {
        stop(
            "'prob_result' must contain columns: ",
            paste(needed, collapse = ", ")
        )
    }

    if (!is.numeric(ci_level) || length(ci_level) != 1L ||
        ci_level <= 0 || ci_level >= 1) {
        stop("'ci_level' must be a single number strictly between 0 and 1.")
    }

    if (length(component_names) != 3L) {
        stop("'component_names' must have length 3.")
    }

    alpha <- 1 - ci_level
    q_probs <- c(alpha / 2, 1 - alpha / 2)

    probs <- as.matrix(prob_result[, needed, drop = FALSE])

    bad_rows <- !apply(probs, 1L, function(z) all(is.finite(z)))
    if (any(bad_rows)) {
        probs <- probs[!bad_rows, , drop = FALSE]
        prob_result <- prob_result[!bad_rows, , drop = FALSE]
    }

    if (nrow(probs) == 0L) {
        stop("No finite probability rows remain after filtering.")
    }

    hard_class <- if (isTRUE(include_hard_calls)) {
        if ("map_component" %in% names(prob_result) &&
            is.numeric(prob_result$map_component)) {
            as.integer(prob_result$map_component)
        } else {
            max.col(probs, ties.method = "first")
        }
    } else {
        NULL
    }

    has_interactions <- "interactions" %in% names(prob_result) &&
        isTRUE(include_interaction_summary)

    out_list <- vector("list", 3L)

    for (k in seq_len(3L)) {
        vals <- probs[, k]

        row_k <- data.frame(
            component = component_names[k],
            mean_probability = mean(vals, na.rm = TRUE),
            median_probability = stats::median(vals, na.rm = TRUE),
            ## These are quantiles of the per-interaction probabilities ACROSS
            ## interactions -- they describe the spread of the fitted
            ## probability field, not uncertainty about any parameter or mean.
            ## They were previously named CI_lower/CI_upper, which invited
            ## reading them as credible intervals; renamed so the output cannot
            ## be misread.
            quantile_lower = stats::quantile(vals,
                probs = q_probs[1],
                na.rm = TRUE, names = FALSE
            ),
            quantile_upper = stats::quantile(vals,
                probs = q_probs[2],
                na.rm = TRUE, names = FALSE
            ),
            stringsAsFactors = FALSE
        )

        if (isTRUE(include_hard_calls)) {
            assigned <- hard_class == k
            row_k$assigned_n <- sum(assigned, na.rm = TRUE)
            row_k$assigned_prop <- mean(assigned, na.rm = TRUE)

            if (has_interactions) {
                assigned_interactions <- prob_result$interactions[assigned]
                if (length(assigned_interactions) > 0L) {
                    row_k$mean_interactions_assigned <- mean(
                        assigned_interactions,
                        na.rm = TRUE
                    )
                    row_k$median_interactions_assigned <- stats::median(
                        assigned_interactions,
                        na.rm = TRUE
                    )
                } else {
                    row_k$mean_interactions_assigned <- NA_real_
                    row_k$median_interactions_assigned <- NA_real_
                }
            }
        }

        out_list[[k]] <- row_k
    }

    out <- do.call(rbind, out_list)
    rownames(out) <- NULL
    out
}

#' Relabel one HiCPotts result by the biological noise relationship
#'
#' Component 1 is the structurally distinct zero-inflated, low-mean noise
#' component and is never reordered. Components 2 and 3 use the same elevated
#' emission family, but the component-1/3 relationship prior gives them fixed
#' biological roles. The sampler first attempts a reversible whole-state 2/3
#' orientation move internally. This function then makes one global 2/3
#' decision for the complete retained chain:
#'
#' \itemize{
#'   \item component 3 is the elevated component whose standardised
#'     covariate-response slopes most closely resemble component 1 and whose
#'     intercept is above component 1;
#'   \item component 2 is the remaining elevated component, whose slopes are
#'     unrestricted and represent true interaction.
#' }
#'
#' The decision uses the same branch score as the sampler's noise-relationship
#' prior: all four standardised covariate slopes and the component-1-to-3
#' elevation are evaluated jointly. It never falls back to component size or
#' baseline intensity alone. A posterior branch probability is stored as a
#' continuous diagnostic; the package does not impose a hard probability
#' threshold.
#'
#' @param res A single HiCPotts result containing three component-specific
#'   coefficient chains.
#' @return A relabelled HiCPotts result. Component allocations, probabilities,
#'   checkpoints, dispersion chains and batch summaries are permuted
#'   consistently. The result also contains \code{label_permutation},
#'   \code{relabel_basis}, and \code{noise_relationship_probability}.
#' @noRd
.hicpotts_log_sigmoid <- function(x) {
    pmin(x, 0) - log1p(exp(-abs(x)))
}

.hicpotts_relationship_settings <- function(res) {
    model <- res$noise_relationship
    value <- function(name, default) {
        x <- if (is.list(model)) model[[name]] else NULL
        if (is.null(x) || length(x) != 1L || !is.finite(
            x
        )) {
            default
        } else {
            as.numeric(x)
        }
    }
    list(
        model_enabled = is.list(model) && isTRUE(model$enabled),
        link_sd = value("link_sd", 0.5),
        order_strength = value("order_strength", 10),
        order_width = value("order_width", 0.5)
    )
}

.hicpotts_noise_branch_score <- function(beta1, candidate, covariate_sds,
                                        link_sd = 0.5,
                                        order_strength = 10,
                                        order_width = 0.5) {
    beta1 <- as.matrix(beta1)
    candidate <- as.matrix(candidate)
    if (!identical(dim(beta1), dim(candidate)) || ncol(beta1) < 1L) {
        stop("The two coefficient inputs must have identical dimensions.")
    }
    n_slopes <- min(ncol(beta1) - 1L, length(covariate_sds))
    slope_score <- rep(0, nrow(beta1))
    if (n_slopes > 0L) {
        columns <- seq_len(n_slopes) + 1L
        differences <- sweep(
            candidate[, columns, drop = FALSE] -
                beta1[, columns, drop = FALSE],
            2L, covariate_sds[seq_len(n_slopes)], "*"
        )
        slope_score <- -0.5 * rowSums((differences / link_sd)^2)
    }
    elevation <- .hicpotts_log_sigmoid(
        (candidate[, 1L] - beta1[, 1L]) / order_width
    )
    slope_score + order_strength * elevation
}

.hicpotts_noise_relationship_scores <- function(res, idx = NULL) {
    ch <- res$chains
    if (!is.list(ch) || length(ch) != 3L ||
        !all(vapply(ch, is.matrix, logical(1)))) {
        stop("res$chains must contain three coefficient-chain matrices.")
    }
    n_iter <- nrow(ch[[1L]])
    if (!all(vapply(ch, nrow, integer(1)) == n_iter)) {
        stop("All component chains must contain the same number of draws.")
    }
    if (is.null(idx)) idx <- (floor(n_iter / 2) + 1L):n_iter
    idx <- as.integer(idx)
    if (!length(idx) || anyNA(idx) || any(idx < 1L | idx > n_iter)) {
        stop("'idx' must select valid coefficient-chain rows.")
    }

    sds <- attr(ch[[1L]], "proposal_covariate_sds", exact = TRUE)
    n_slopes <- max(0L, min(4L, ncol(ch[[1L]]) - 1L))
    if (is.null(sds) || length(sds) < n_slopes ||
        any(!is.finite(sds[seq_len(n_slopes)])) ||
        any(sds[seq_len(n_slopes)] <= 0)) {
        sds <- rep(1, n_slopes)
        scale_source <- "unit-scale fallback (legacy fit)"
    } else {
        sds <- as.numeric(sds[seq_len(n_slopes)])
        scale_source <- "recorded log1p-covariate standard deviations"
    }
    settings <- .hicpotts_relationship_settings(res)
    score2 <- .hicpotts_noise_branch_score(
        ch[[1L]][idx, , drop = FALSE], ch[[2L]][idx, , drop = FALSE], sds,
        settings$link_sd, settings$order_strength, settings$order_width
    )
    score3 <- .hicpotts_noise_branch_score(
        ch[[1L]][idx, , drop = FALSE], ch[[3L]][idx, , drop = FALSE], sds,
        settings$link_sd, settings$order_strength, settings$order_width
    )
    list(
        score2 = score2, score3 = score3, delta = score3 - score2,
        covariate_sds = sds, scale_source = scale_source, settings = settings
    )
}

.hicpotts_swap_labels23 <- function(z) {
    if (!is.matrix(z)) {
        return(z)
    }
    out <- z
    out[z == 2L] <- 3L
    out[z == 3L] <- 2L
    out
}

relabel_hicpotts_result <- function(res) {
    scores <- .hicpotts_noise_relationship_scores(res)
    ## Direct posterior probability that component 3 has the stronger
    ## component-1-similarity-plus-elevation relationship. With labelled
    ## components this event probability is the relevant identification
    ## diagnostic; the logistic branch weight used by the earlier exchangeable
    ## prototype is not.
    branch3_probability <- mean(scores$delta > 0, na.rm = TRUE)
    if (!is.finite(branch3_probability)) {
        stop("The component-2/3 biological relationship score is not finite.")
    }
    ## The relationship prior and the per-iteration reversible branch move do
    ## most of the orientation work. A final global decision is nevertheless
    ## retained so a whole chain that remains in the opposite 2/3 branch is
    ## reported with the intended biological labels. Every stored state and
    ## allocation summary is permuted together below.
    swap23 <- branch3_probability < 0.5
    selected_probability <- if (swap23) {
        1 - branch3_probability
    } else {
        branch3_probability
    }
    oriented_delta <- if (swap23) -scores$delta else scores$delta
    delta_ci <- stats::quantile(oriented_delta, c(0.025, 0.975),
        names = FALSE, na.rm = TRUE
    )

    out <- res
    n_iter <- nrow(res$chains[[1L]])
    if (swap23) {
        out$chains[c(2L, 3L)] <- res$chains[c(3L, 2L)]
        if (is.matrix(res$size) && nrow(res$size) == 3L) {
            out$size[c(2L, 3L), ] <- res$size[c(3L, 2L), , drop = FALSE]
        }

        if (is.matrix(res$z_final)) {
            out$z_final <- .hicpotts_swap_labels23(res$z_final)
            for (attribute_name in c(
                "parameter_component_counts",
                "classification_component_counts"
            )) {
                counts <- attr(res$z_final, attribute_name, exact = TRUE)
                if (!is.null(counts) && length(counts) == 3L) {
                    attr(out$z_final, attribute_name) <- counts[c(1L, 3L, 2L)]
                }
            }
            parameter_z <- attr(res$z_final, "z_parameter_final", exact = TRUE)
            if (is.matrix(parameter_z)) {
                attr(out$z_final, "z_parameter_final") <-
                    .hicpotts_swap_labels23(parameter_z)
            }
            forced_z <- attr(res$z_final, "zero_forced_classification",
                exact = TRUE
            )
            if (is.matrix(forced_z)) {
                attr(out$z_final, "zero_forced_classification") <-
                    .hicpotts_swap_labels23(forced_z)
            }
        }
        if (is.matrix(res$z_parameter_final)) {
            out$z_parameter_final <- .hicpotts_swap_labels23(
                res$z_parameter_final
            )
        }
        if (is.list(res$z_checkpoints)) {
            out$z_checkpoints <- lapply(
                res$z_checkpoints,
                .hicpotts_swap_labels23
            )
        }
        if (is.list(res$z_probabilities) && length(res$z_probabilities) == 3L &&
            all(vapply(res$z_probabilities, is.matrix, logical(1)))) {
            out$z_probabilities <- res$z_probabilities[c(1L, 3L, 2L)]
            names(out$z_probabilities) <- paste0("component", seq_len(3L))
        }
        if (is.list(res$z_probability_batches) &&
            length(dim(res$z_probability_batches$batch_means)) == 4L) {
            out$z_probability_batches <- res$z_probability_batches
            out$z_probability_batches$batch_means <-
                res$z_probability_batches$batch_means[, , c(1L, 3L, 2L), ,
                    drop = FALSE
                ]
        }
        if (is.list(res$beta_mixing)) {
            out$beta_mixing <- res$beta_mixing
            for (field in c(
                "attempts", "acceptances", "acceptance_rate",
                "warmup_attempts", "warmup_acceptances",
                "retained_attempts", "retained_acceptances",
                "retained_acceptance_rate"
            )) {
                value <- res$beta_mixing[[field]]
                if (!is.null(value) && length(value) == 3L) {
                    out$beta_mixing[[field]] <- value[c(1L, 3L, 2L)]
                }
            }
            step <- res$beta_mixing$final_qr_step
            if (is.matrix(step) && nrow(step) == 3L) {
                out$beta_mixing$final_qr_step <-
                    step[c(1L, 3L, 2L), , drop = FALSE]
            }
        }
        if (is.list(res$regression_prior) &&
            is.list(res$regression_prior$frozen_hyperparameters) &&
            length(res$regression_prior$frozen_hyperparameters) == 3L) {
            out$regression_prior <- res$regression_prior
            out$regression_prior$frozen_hyperparameters <-
                res$regression_prior$frozen_hyperparameters[c(1L, 3L, 2L)]
            names(out$regression_prior$frozen_hyperparameters) <-
                paste0("component", seq_len(3L))
        }
    }

    permutation <- if (swap23) c(1L, 3L, 2L) else c(1L, 2L, 3L)
    out$label_permutation <- matrix(
        rep(permutation, each = n_iter),
        nrow = n_iter, ncol = 3L,
        dimnames = list(NULL, c(
            "new_comp1_from", "new_comp2_from",
            "new_comp3_from"
        ))
    )
    out$relabel_basis <- paste0(
        "all four standardised covariate slopes relative to component 1, ",
        "plus elevated-noise intercept preference"
    )
    out$noise_relationship_probability <- selected_probability
    out$noise_relationship_score_difference_ci <- unname(delta_ci)
    relationship <- if (is.list(res$noise_relationship)) {
        res$noise_relationship
    } else {
        list(enabled = FALSE)
    }
    relationship$selected_probability <- selected_probability
    relationship$score_difference_ci <- unname(delta_ci)
    relationship$threshold_policy <- paste0(
        "none in the package; downstream users choose any probability threshold"
    )
    relationship$covariate_sds <- scores$covariate_sds
    relationship$scale_source <- scores$scale_source
    out$noise_relationship <- relationship
    out$relabel_rule <- paste0(
        "Component 1 retained as low-mean noise. Components 2/3 ",
        if (swap23) {
            "swapped globally after the internal branch moves"
        } else {
            "left in the internally selected orientation"
        },
        "; component 3 is tested as the elevated component most consistent ",
        "with ",
        "component 1's standardised slopes and component 2 remains ",
        "unrestricted. ",
        "Component-3 relationship probability = ",
        formatC(selected_probability, digits = 4L, format = "f"),
        "; no hard probability threshold is imposed by the package."
    )
    out$component_definition <- hicpotts_component_definition()
    out
}

#' Relabel HiCPotts output
#'
#' Applies the component-1/3 biological noise relationship to either one
#' HiCPotts result or a list of independent results. Every chain is labelled by
#' the same scientific rule; no allocation-agreement or population-size
#' tie-break can override it.
#'
#' @param x A single HiCPotts result or a list of HiCPotts results.
#' @return A relabelled object with the same outer structure as \code{x}.
#' @examples
#' make_chain <- function(intercept, slopes) {
#'     cbind(
#'         intercept = rep(intercept, 4L),
#'         matrix(rep(slopes, each = 4L), nrow = 4L)
#'     )
#' }
#' fit <- list(
#'     chains = list(
#'         make_chain(0, rep(0, 4L)),
#'         make_chain(2, rep(1, 4L)),
#'         make_chain(3, rep(0, 4L))
#'     ),
#'     size = matrix(1, 3, 4),
#'     z_final = matrix(c(1, 2, 3, 1), 2, 2),
#'     z_checkpoints = list(),
#'     z_probabilities = lapply(seq_len(3L), function(k) matrix(1 / 3, 2, 2))
#' )
#' relabelled <- relabel_hicpotts(fit)
#' relabelled$component_definition
#' @export
relabel_hicpotts <- function(x) {
    is_single_result <- is.list(x) && !is.null(x$chains) &&
        is.list(x$chains) && length(x$chains) == 3L
    if (is_single_result) {
        return(relabel_hicpotts_result(x))
    }

    if (is.list(x) && length(x) > 0L) {
        valid <- vapply(x, function(element) {
            is.list(element) && is.list(element$chains) &&
                length(element$chains) == 3L
        }, logical(1))
        if (all(valid)) {
            return(lapply(x, relabel_hicpotts_result))
        }
    }
    stop("x must be either one HiCPotts result or a list of HiCPotts results.")
}

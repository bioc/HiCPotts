#' Officially classify interactions from sampled HiCPotts latent states
#'
#' Combines the post-burn-in latent-state membership frequencies stored by the
#' HiCPotts sampler across one or more relabelled chains. This is the package's
#' official three-state classification and a
#' classification-layer summary of the fitted model: it does not refit the
#' model, change an MCMC transition, or apply additional spatial smoothing.
#'
#' Unlike a classification based on the single final state or on posterior-mean
#' parameters, the returned probabilities average the latent allocations
#' actually sampled under the full count and Potts model. Each cell is assigned
#' to its maximum-posterior component; probability margin and entropy are
#' returned as diagnostics.
#'
#' @param fit A raw chain, list of chains, \code{hicpotts_robust_fit}, or
#'   block-aware result from \code{combine_hicpotts_blocks()}.
#' @param data Optional data frame in the same column-major lattice order as
#'   \code{process_data()}. It must contain one row per lattice cell. When
#'   omitted, row and column indices are returned.
#' @param component_names Character vector of length three giving the
#'   semantic names of components 1, 2 and 3, in that order. The default
#'   is the package's canonical naming,
#'   \code{c("noise", "signal", "false signal")}, where component 1 is
#'   low-mean noise, component 2 is a true interaction with an
#'   unrestricted covariate-response pattern, and component 3 is elevated
#'   noise whose covariate-response pattern resembles component 1. Supply
#'   your own wording to rename the reported classes; this relabels the
#'   output only and does not change which component a cell is assigned
#'   to. The default is produced internally, so it is shown in the usage
#'   above as a function call rather than as a literal vector.
#' @param use For a robust fit, use the mode-consistent chains selected for
#'   parameter summaries (\code{"selected"}, the default) or every retained
#'   chain (\code{"all"}). The selected result is conditional on that posterior
#'   mode; a warning is emitted when another replicated mode exists.
#' @param relabel Whether to apply the package's component-2/3 relabelling and
#'   cross-chain allocation alignment before pooling. The default is
#'   \code{TRUE}.
#' @param min_draws Minimum number of post-burn-in latent-state draws, pooled
#'   across chains, required before a classification is returned. Classifying
#'   from very few draws produces near-degenerate probabilities whose margin and
#'   entropy diagnostics are not meaningful, so this is an error rather than a
#'   warning. There is no fallback to the single final latent state. Defaults to
#'   \code{100}; lower it deliberately if a coarser classification is
#'   acceptable.
#' @param reflect Reflection averaging for mirrored contact matrices. On a
#'   symmetric input, cells \code{(i,j)} and \code{(j,i)} are two stored copies
#'   of one measured contact and must receive the same label; averaging their
#'   pooled probabilities guarantees this. \code{"auto"} (the default) applies
#'   it when every fit records a confirmed symmetric input and the lattice is
#'   square, \code{"always"} forces it, \code{"never"} disables it.
#'
#' @return A data frame containing \code{prob1}, \code{prob2}, \code{prob3},
#'   the maximum-posterior component and label, the winning probability,
#'   probability margin, normalized entropy, and the classification. Attributes
#'   record the number of chains and latent-state draws pooled and the
#'   classification scope.
#'
#'
#' @section Working with sub-blocks of a large map:
#' This function summarises one contact map. It cannot span sub-blocks:
#' chains fitted to different blocks have different data and different
#' frozen priors, so they are refused rather than pooled, and passing a
#' whole-map \code{data} frame to a single block's fit fails the
#' one-row-per-lattice-cell check.
#'
#' Run it per block and stitch the results, which lets every block keep
#' the chains its own mode screen accepted -- block 1 may be summarised
#' from chains 1 and 3 while block 2 uses chains 1 and 2:
#' \preformatted{parts <- list()
#' for (j in seq_along(fits)) {
#'     parts[[j]] <- classify_hicpotts(fits[[j]], data = blocks[[j]])
#'     parts[[j]]$block <- names(fits)[j]
#' }
#' combined <- do.call(rbind, parts)}
#' \code{\link{combine_hicpotts_blocks}()} does this with the checks that
#' matter (matching lengths, one row per cell, non-overlapping blocks) and
#' records the chains each block used. Note that component labels are only
#' aligned within a block, not across blocks.
#' @examples
#' N <- 3L
#' probs <- list(
#'     component1 = matrix(0.8, N, N),
#'     component2 = matrix(0.15, N, N),
#'     component3 = matrix(0.05, N, N)
#' )
#' fake_fit <- list(
#'     chains = lapply(1:3, function(k) matrix(k, 4, 5)),
#'     size = matrix(1, 3, 4),
#'     z_final = matrix(1, N, N),
#'     z_checkpoints = list(),
#'     z_probabilities = probs,
#'     z_probability_draws = 200L
#' )
#' classify_hicpotts(fake_fit, relabel = FALSE)
#'
#' @seealso \code{\link{compute_HMRFHiC_probabilities}},
#'   \code{\link{relabel_hicpotts}}
#' @export
classify_hicpotts <- function(
    fit,
    data = NULL,
    component_names = hicpotts_component_definition()$label,
    use = c("selected", "all"),
    relabel = TRUE,
    min_draws = 100L,
    reflect = c("auto", "always", "never")
) {
    use <- match.arg(use)
    reflect <- match.arg(reflect)
    if (!is.character(component_names) || length(component_names) != 3L ||
        anyNA(component_names) || any(!nzchar(component_names))) {
        stop("'component_names' must contain three non-empty character labels.")
    }
    if (anyDuplicated(component_names)) {
        stop(
            "'component_names' must be three DISTINCT labels; duplicates ",
            "would collapse two biological states into one class."
        )
    }
    if (!is.logical(relabel) || length(relabel) != 1L || is.na(relabel)) {
        stop("'relabel' must be TRUE or FALSE.")
    }
    if (!is.numeric(min_draws) || length(min_draws) != 1L ||
        !is.finite(min_draws) || min_draws < 1) {
        stop("'min_draws' must be a single finite number >= 1.")
    }
    if (inherits(fit, "hicpotts_block_fit")) {
        return(.hicpotts_classify_blocks(
            fit = fit, data = data, component_names = component_names,
            use = use, relabel = relabel, min_draws = min_draws,
            reflect = reflect
        ))
    }

    robust <- inherits(fit, "hicpotts_robust_fit")
    scope <- "supplied fits"
    if (robust) {
        fits <- if (identical(use, "all")) fit$all_fits else fit$fits
        scope <- if (identical(use, "all")) {
            "all robust-fit chains"
        } else {
            "selected allocation mode"
        }
        if (identical(use, "selected") && isTRUE(
            fit$mode_selection$multimodal
        )) {
            warning(
                "Classification is conditional on the selected allocation ",
                "mode; ",
                "the robust fit contains at least one replicated competing ",
                "mode."
            )
        }
        if (!all(fit$diagnostics$reliability_flags$passed)) {
            warning(
                "Review the recorded diagnostic criteria in ",
                "fit$diagnostics$reliability_flags. Criteria below their ",
                "configured thresholds: ",
                paste(fit$diagnostics$reliability_flags$criterion[
                    !fit$diagnostics$reliability_flags$passed
                ], collapse = ", "), "."
            )
        }
    } else {
        is_single <- is.list(fit) && is.list(fit$chains) && length(
            fit$chains
        ) == 3L
        fits <- if (is_single) list(fit) else fit
    }
    if (!is.list(fits) || !length(fits) ||
        !all(vapply(
            fits, function(x) {
                is.list(x) && is.list(x$chains) && length(x$chains) == 3L
            },
            logical(1)
        ))) {
        stop(
            "'fit' must be a HiCPotts result, a non-empty fit list, or a ",
            "robust fit."
        )
    }

    ## Refuse to pool chains fitted to different data. Shapes alone are not
    ## enough: two unrelated datasets on the same lattice size would pool
    ## silently and return a meaningless answer.
    .hicpotts_check_poolable(fits, context = "classification pooling")

    if (isTRUE(relabel)) fits <- relabel_hicpotts(fits)

    valid_probability_list <- function(x) {
        is.list(x) && length(x) == 3L &&
            all(vapply(x, function(m) is.matrix(m) && is.numeric(m), logical(
                1
            )))
    }
    have_probabilities <- vapply(fits, function(x) {
        valid_probability_list(x$z_probabilities)
    }, logical(1))
    if (!all(have_probabilities)) {
        stop(
            "Every fit must contain $z_probabilities. Refit the data to ",
            "record post-burn-in latent-state membership frequencies."
        )
    }

    reference_dim <- dim(fits[[1L]]$z_probabilities[[1L]])
    if (length(reference_dim) != 2L || any(reference_dim < 1L)) {
        stop(
            "Stored z-probability matrices must be non-empty two-dimensional ",
            "matrices."
        )
    }
    for (i in seq_along(fits)) {
        zp <- fits[[i]]$z_probabilities
        if (!all(vapply(
            zp, function(m) identical(dim(m), reference_dim),
            logical(1)
        ))) {
            stop(
                "All stored z-probability matrices must have identical ",
                "dimensions."
            )
        }
        values <- unlist(zp, use.names = FALSE)
        if (any(!is.finite(values)) || any(values < -1e-10) || any(
            values > 1 + 1e-10
        )) {
            stop("Stored z probabilities must be finite and lie in [0, 1].")
        }
        total <- Reduce(`+`, zp)
        if (any(abs(total - 1) > 1e-6)) {
            stop(
                "Stored component probabilities must sum to one at every ",
                "lattice cell."
            )
        }
    }

    ## Post-burn-in draw counts. There is deliberately NO fallback here: a fit
    ## that cannot report how many valid latent-state draws it accumulated is
    ## rejected rather than classified from a single retained state. Classifying
    ## from one final allocation yields degenerate 0/1 "probabilities" that
    ## carry no posterior uncertainty, and the margin and entropy diagnostics
    ## computed from them are meaningless.
    draws <- vapply(fits, function(x) {
        value <- x$z_probability_draws
        if (is.null(value) || length(value) != 1L || !is.finite(
            value
        ) || value < 1) {
            stop(
                "Each fit must report a positive $z_probability_draws value. ",
                "Refit the data; classification from a single final latent ",
                "state is not supported."
            )
        }
        as.numeric(value)
    }, numeric(1))
    total_draws <- sum(draws)
    if (total_draws < min_draws) {
        template <- paste0(
            "Only %g post-burn-in latent-state draws are available across %d ",
            "chain(s), below the required minimum of %g. Run a longer chain, ",
            "lower the burn-in, or lower 'min_draws' if a coarser ",
            "classification is acceptable."
        )
        stop(
            sprintf(template, total_draws, length(fits), as.numeric(
                min_draws
            )),
            call. = FALSE
        )
    }
    pooled <- lapply(seq_len(3L), function(k) {
        numerator <- Reduce(`+`, Map(function(x, weight) {
            x$z_probabilities[[k]] * weight
        }, fits, draws))
        numerator / sum(draws)
    })

    ## ---- reflection averaging for symmetric inputs
    ## --------------------------- On a mirrored contact matrix (i,j) and (j,i)
    ## are two stored copies of ONE measured contact, so they must receive the
    ## same classification. The sampler treats z as a full lattice and updates
    ## the two cells separately, so their pooled probabilities differ by
    ## Monte-Carlo noise and can even straddle a decision boundary, producing
    ## two different labels for a single contact. Averaging the two orientations
    ## removes that, and because the average is taken component-wise the three
    ## probabilities still sum to one.
    ##
    ## Applied only when the fit records a confirmed symmetric input (or when
    ## explicitly requested), so genuinely directional data is never folded.
    symmetric_fit <- vapply(fits, function(x) {
        isTRUE(x$pair_weighting$symmetric_input)
    }, logical(1))
    reflection_applied <- switch(reflect,
        always = TRUE,
        never = FALSE,
        auto = all(symmetric_fit) && length(symmetric_fit) > 0L &&
            identical(reference_dim[1L], reference_dim[2L])
    )
    if (reflection_applied) {
        if (!identical(reference_dim[1L], reference_dim[2L])) {
            stop("Reflection averaging requires a square lattice.")
        }
        pooled <- lapply(pooled, function(m) (m + t(m)) / 2)
    }

    probability_matrix <- cbind(
        prob1 = as.vector(pooled[[1L]]),
        prob2 = as.vector(pooled[[2L]]),
        prob3 = as.vector(pooled[[3L]])
    )
    probability_matrix <- probability_matrix / rowSums(probability_matrix)
    ## Every contact is assigned to exactly one of the three biological states
    ## by maximum posterior probability. ties.method = "first" makes the
    ## assignment total and deterministic: there is no tie case that could yield
    ## NA, no "uncertain"/"ambiguous" bucket, and no binary collapse of the two
    ## signal states. Confidence is reported separately (class_probability,
    ## probability_margin, normalized_entropy) and never withheld as a class.
    map_component <- max.col(probability_matrix, ties.method = "first")
    if (anyNA(map_component) || any(map_component < 1L) || any(
        map_component > 3L
    )) {
        stop("MAP assignment produced a component outside 1:3.", call. = FALSE)
    }
    sorted <- t(apply(probability_matrix, 1L, sort, decreasing = TRUE))
    winning_probability <- sorted[, 1L]
    probability_margin <- sorted[, 1L] - sorted[, 2L]
    entropy_terms <- ifelse(probability_matrix > 0,
        probability_matrix * log(probability_matrix), 0
    )
    normalized_entropy <- -rowSums(entropy_terms) / log(3)
    map_label <- component_names[map_component]

    n_cells <- prod(reference_dim)
    if (is.null(data)) {
        index <- arrayInd(seq_len(n_cells), .dim = reference_dim)
        out <- data.frame(row = index[, 1L], column = index[, 2L])
    } else {
        if (!is.data.frame(data)) data <- as.data.frame(data)
        if (nrow(data) != n_cells) {
            stop("'data' must contain exactly one row per stored lattice cell.")
        }
        ## The probabilities are attached to 'data' positionally, so 'data' must
        ## be in the SAME column-major lattice order the sampler stored. A
        ## row-shuffled frame has the right number of rows and would silently
        ## attach every probability to the wrong genomic coordinate. When the
        ## frame carries lattice indices, verify the order rather than trusting
        ## it.
        idx_cols <- intersect(c("row", "column"), names(data))
        if (length(idx_cols) == 2L) {
            expected <- arrayInd(seq_len(n_cells), .dim = reference_dim)
            if (!identical(as.integer(data$row), as.integer(expected[, 1L])) ||
                !identical(as.integer(data$column), as.integer(expected[
                    ,
                    2L
                ]))) {
                stop(
                    "'data' is not in the stored column-major lattice order ",
                    "(its ",
                    "row/column indices do not match). Reorder it to match ",
                    "process_data() output; attaching probabilities to ",
                    "shuffled rows ",
                    "would mislabel every contact."
                )
            }
        }
        out <- data
    }
    out$prob1 <- probability_matrix[, 1L]
    out$prob2 <- probability_matrix[, 2L]
    out$prob3 <- probability_matrix[, 3L]
    out$map_component <- map_component
    out$map_label <- map_label
    out$class_probability <- winning_probability
    out$probability_margin <- probability_margin
    out$normalized_entropy <- normalized_entropy
    ## A factor with exactly the three biological levels, so downstream code
    ## (tables, plots, joins) cannot silently introduce a fourth category and an
    ## empty class stays visible rather than disappearing.
    out$classification <- factor(map_label, levels = component_names)
    if (anyNA(out$classification)) {
        stop("A contact was not assigned one of the three labels.",
            call. = FALSE
        )
    }

    attr(out, "classification_source") <-
        "post-burn-in unrestricted latent-state MCMC frequencies"
    attr(out, "classification_scope") <- scope
    attr(out, "chains_combined") <- length(fits)
    attr(out, "posterior_state_draws") <- sum(draws)
    attr(out, "min_draws_required") <- as.numeric(min_draws)
    attr(out, "component_names") <- component_names
    attr(out, "reflection_averaged") <- reflection_applied
    class(out) <- c("hicpotts_classification", "data.frame")
    out
}

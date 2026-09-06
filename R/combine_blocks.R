#' Combine independently fitted HiCPotts matrix blocks
#'
#' @description
#' Create one block-aware HiCPotts object from fits to submatrices of a larger
#' contact map. Each block retains the chains selected by its own mode screen.
#' Downstream functions analyse every block independently and combine only
#' compatible cell-level results.
#'
#' @details
#' MCMC draws from different matrix blocks target different posteriors and are
#' retained separately. For a robust fit, \code{mode_selection$selected}
#' is reapplied
#' to \code{all_fits}; consequently block 1 can use chains 1 and 3 while block
#' 2 uses chain 2. The selected chains are recorded in \code{block_table}.
#' \code{summarise_hicpotts_parameters(pool_blocks = "posterior_weighted")}
#' can additionally calculate a descriptive, occupancy-aware average from
#' independent posterior draws. It is explicitly labelled as an aggregate and
#' is not a substitute for a joint shared-parameter fit.
#'
#' A single fit is returned unchanged when \code{what = "fit"}. Thus callers
#' do not need a special branch for an analysis containing only one matrix.
#'
#' The structured output of \code{\link{process_data}} may be supplied through
#' \code{processed}. It is sufficient for classification, allocation
#' diagnostics, parameter summaries and posterior predictive checks. The
#' parameter-based probability sensitivity function additionally needs the
#' original data frame(s), because genomic coordinates and untransformed
#' covariates cannot be reconstructed from processed matrices.
#'
#' @param fits A fitted result, or a list containing one fitted result per
#'   matrix block. Names, when present, label the blocks.
#' @param data Optional original data frame, or list of one data frame per
#'   block. For multiple blocks, a whole-map data frame is split sequentially
#'   using the fitted lattice sizes.
#' @param processed Optional structured result from \code{\link{process_data}}.
#'   For convenience, it may be supplied as the second positional argument.
#' @param what Return a reusable block-aware \code{"fit"} (the default), or
#'   immediately return combined \code{"classification"} or
#'   \code{"probabilities"} output.
#' @param ... Arguments passed to \code{\link{classify_hicpotts}} or
#'   \code{\link{compute_HMRFHiC_probabilities}} when requested by \code{what}.
#'
#' @return For multiple blocks and \code{what = "fit"}, an object of class
#'   \code{hicpotts_block_fit}. For one block, the original fit is returned.
#'   The other values of \code{what} return a combined data frame with block
#'   provenance.
#'
#' @seealso \code{\link{classify_hicpotts}},
#'   \code{\link{compute_HMRFHiC_probabilities}},
#'   \code{\link{summarise_hicpotts_parameters}}
#'
#' @examples
#' make_example_fit <- function(component, N = 2L) {
#'     probabilities <- lapply(1:3, function(k) {
#'         matrix(if (k == component) 0.9 else 0.05, N, N)
#'     })
#'     list(
#'         chains = lapply(1:3, function(k) matrix(k, 4L, 5L)),
#'         z_final = matrix(component, N, N),
#'         z_probabilities = probabilities,
#'         z_probability_draws = 20L
#'     )
#' }
#' block_fit <- combine_hicpotts_blocks(list(
#'     left = make_example_fit(1L), right = make_example_fit(2L)
#' ))
#' classified <- classify_hicpotts(
#'     block_fit, relabel = FALSE, min_draws = 1L
#' )
#' table(classified$block, classified$map_component)
#'
#' @export
combine_hicpotts_blocks <- function(
    fits, data = NULL, processed = NULL,
    what = c("fit", "classification", "probabilities"), ...
) {
    what <- match.arg(what)
    if (.hicpotts_is_processed(data) && is.null(processed)) {
        processed <- data
        data <- NULL
    }

    single_direct <- .hicpotts_is_fit(fits)
    if (single_direct) {
        block_fits <- list(fits)
    } else {
        if (!is.list(fits) || !length(fits) ||
            !all(vapply(fits, .hicpotts_is_fit, logical(1)))) {
            stop(
                "'fits' must be a fitted result or a non-empty list of ",
                "per-block fitted results.",
                call. = FALSE
            )
        }
        block_fits <- fits
    }
    labels <- names(block_fits)
    if (is.null(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
        labels <- paste0("block", seq_along(block_fits))
    }
    names(block_fits) <- labels
    block_fits <- lapply(block_fits, .hicpotts_prepare_block_fit)

    dimensions <- lapply(block_fits, .hicpotts_fit_dimension)
    cells <- as.integer(vapply(dimensions, prod, numeric(1)))
    block_data <- .hicpotts_normalise_block_data(data, cells, labels)
    processed <- .hicpotts_validate_processed(
        processed, length(block_fits), dimensions
    )

    ## A one-element wrapper should behave exactly like the ordinary fit. This
    ## is important for pipelines that discover the number of blocks at run
    ## time and always call this function.
    if (length(block_fits) == 1L && identical(what, "fit")) {
        return(block_fits[[1L]])
    }

    provenance <- do.call(rbind, lapply(seq_along(block_fits), function(j) {
        info <- .hicpotts_block_chains(block_fits[[j]])
        reliable <- if (inherits(block_fits[[j]], "hicpotts_robust_fit")) {
            flags <- block_fits[[j]]$diagnostics$reliability_flags
            is.data.frame(flags) && nrow(flags) > 0L &&
                all(flags$passed %in% TRUE)
        } else {
            NA
        }
        data.frame(
            block = labels[j], block_index = j,
            rows = dimensions[[j]][1L],
            columns = dimensions[[j]][2L], cells = cells[j],
            selected_chains = paste(info$selected, collapse = ","),
            chains_run = info$total,
            competing_mode = info$multimodal,
            reliable = reliable,
            stringsAsFactors = FALSE
        )
    }))
    combined <- structure(
        list(
            blocks = block_fits,
            block_table = provenance,
            mode_selection = lapply(block_fits, function(x) {
                if (inherits(x, "hicpotts_robust_fit")) {
                    x$mode_selection
                } else {
                    NULL
                }
            }),
            diagnostics = lapply(block_fits, function(x) x$diagnostics),
            data = block_data,
            processed = processed
        ),
        class = "hicpotts_block_fit"
    )
    if (identical(what, "classification")) {
        return(classify_hicpotts(combined, ...))
    }
    if (identical(what, "probabilities")) {
        return(compute_HMRFHiC_probabilities(
            data = data, chain_betas = combined, ...
        ))
    }
    combined
}

#' @export
print.hicpotts_block_fit <- function(x, ...) {
    cat("HiCPotts block fit\n")
    cat("  blocks:", length(x$blocks), "\n")
    print(x$block_table, row.names = FALSE)
    invisible(x)
}

.hicpotts_is_fit <- function(x) {
    inherits(x, "hicpotts_robust_fit") ||
        (is.list(x) && is.list(x$chains) && length(x$chains) == 3L)
}

.hicpotts_is_processed <- function(x) {
    is.list(x) && is.list(x$x_vars) && is.list(x$y) &&
        all(c("distance", "GC", "TES", "ACC") %in% names(x$x_vars))
}

## Make mode_selection authoritative even if a caller has retained an older or
## manually modified robust object whose $fits field is stale.
.hicpotts_prepare_block_fit <- function(fit) {
    if (!inherits(fit, "hicpotts_robust_fit")) return(fit)
    selected <- fit$mode_selection$selected
    if (is.null(selected)) {
        selected <- seq_along(fit$fits)
        fit$mode_selection$selected <- selected
        return(fit)
    }
    selected <- as.integer(selected)
    if (!length(selected) || anyNA(selected) || any(selected < 1L) ||
        any(selected > length(fit$all_fits)) || anyDuplicated(selected)) {
        stop(
            "A block has invalid mode_selection$selected chain indices.",
            call. = FALSE
        )
    }
    fit$fits <- fit$all_fits[selected]
    fit
}

.hicpotts_fit_dimension <- function(fit) {
    chains <- if (inherits(fit, "hicpotts_robust_fit")) fit$fits else list(fit)
    first <- chains[[1L]]
    z <- first$z_probabilities
    d <- if (is.list(z) && length(z) && is.matrix(z[[1L]])) {
        dim(z[[1L]])
    } else if (is.matrix(first$z_final)) {
        dim(first$z_final)
    } else {
        NULL
    }
    if (is.null(d) || length(d) != 2L || any(!is.finite(d)) || any(d < 1L)) {
        stop("Cannot determine a block's fitted lattice dimensions.",
            call. = FALSE
        )
    }
    as.integer(d)
}

.hicpotts_validate_processed <- function(processed, n_blocks, dimensions) {
    if (is.null(processed)) return(NULL)
    if (!.hicpotts_is_processed(processed)) {
        stop("'processed' must be an object returned by process_data().",
            call. = FALSE
        )
    }
    if (length(processed$y) != n_blocks) {
        stop(
            sprintf(
                "'processed' contains %d matrices but 'fits' contains %d.",
                length(processed$y), n_blocks
            ),
            call. = FALSE
        )
    }
    for (j in seq_len(n_blocks)) {
        y <- processed$y[[j]]
        if (!is.matrix(y) || !identical(dim(y), dimensions[[j]])) {
            stop(sprintf(
                "processed$y[[%d]] does not match fitted block dimensions.", j
            ), call. = FALSE)
        }
        for (nm in c("distance", "GC", "TES", "ACC")) {
            value <- processed$x_vars[[nm]][[j]]
            if (!is.matrix(value) || !identical(dim(value), dimensions[[j]])) {
                stop(sprintf(
                    "processed$x_vars$%s[[%d]] differs from its fitted block.",
                    nm, j
                ), call. = FALSE)
            }
        }
    }
    processed
}

.hicpotts_normalise_block_data <- function(data, cells, labels) {
    if (is.null(data)) return(NULL)
    if (is.data.frame(data)) {
        if (nrow(data) != sum(cells)) {
            stop(sprintf(
                "'data' has %d rows; the fitted blocks require %d cells.",
                nrow(data), sum(cells)
            ), call. = FALSE)
        }
        ends <- cumsum(cells)
        starts <- c(1L, utils::head(ends, -1L) + 1L)
        data <- Map(function(a, b) data[a:b, , drop = FALSE], starts, ends)
    }
    if (!is.list(data) || length(data) != length(cells) ||
        !all(vapply(data, is.data.frame, logical(1)))) {
        stop("'data' must provide one data frame per fitted block.",
            call. = FALSE
        )
    }
    observed <- vapply(data, nrow, integer(1))
    if (!identical(unname(observed), as.integer(cells))) {
        stop(
            "Each data block must contain one row per fitted lattice cell.",
            call. = FALSE
        )
    }
    names(data) <- labels
    data
}

.hicpotts_processed_xvars <- function(processed, j) {
    lapply(processed$x_vars, function(value) list(value[[j]]))
}

.hicpotts_block_chains <- function(fit) {
    if (inherits(fit, "hicpotts_robust_fit")) {
        selected <- fit$mode_selection$selected
        if (is.null(selected)) selected <- seq_along(fit$fits)
        return(list(
            selected = as.integer(selected), total = length(fit$all_fits),
            multimodal = isTRUE(fit$mode_selection$multimodal)
        ))
    }
    list(selected = 1L, total = 1L, multimodal = FALSE)
}

.hicpotts_block_data <- function(fit, data = NULL, require = FALSE) {
    if (!is.null(data)) {
        return(.hicpotts_normalise_block_data(
            data, fit$block_table$cells, fit$block_table$block
        ))
    }
    if (!is.null(fit$data)) return(fit$data)
    if (isTRUE(require)) {
        stop(
            "Original per-cell data are required for this calculation. ",
            "Supply them to combine_hicpotts_blocks(data = ...) or to this ",
            "function; process_data() matrices alone do not retain genomic ",
            "coordinates.",
            call. = FALSE
        )
    }
    NULL
}

.hicpotts_bind_block_tables <- function(parts, fit, output_class = NULL) {
    for (j in seq_along(parts)) {
        parts[[j]]$block <- fit$block_table$block[j]
        parts[[j]]$block_index <- j
        parts[[j]]$block_cell_index <- seq_len(nrow(parts[[j]]))
    }
    combined <- do.call(rbind, lapply(parts, as.data.frame))
    rownames(combined) <- NULL
    combined$global_cell_index <- seq_len(nrow(combined))
    .hicpotts_check_block_overlap(combined)
    attr(combined, "blocks") <- fit$block_table
    if (!is.null(output_class)) class(combined) <- c(output_class, "data.frame")
    combined
}

.hicpotts_classify_blocks <- function(fit, data = NULL, ...) {
    block_data <- .hicpotts_block_data(fit, data, require = FALSE)
    parts <- lapply(seq_along(fit$blocks), function(j) {
        classify_hicpotts(
            fit$blocks[[j]],
            data = if (is.null(block_data)) NULL else block_data[[j]], ...
        )
    })
    out <- .hicpotts_bind_block_tables(
        parts, fit, "hicpotts_classification"
    )
    attr(out, "classification_source") <-
        "block-specific post-burn-in latent-state MCMC frequencies"
    attr(out, "classification_scope") <-
        "each block's independently selected allocation mode"
    attr(out, "chains_by_block") <- fit$block_table[, c(
        "block", "selected_chains", "chains_run"
    )]
    out
}

.hicpotts_probability_blocks <- function(
    fit, data = NULL, iterations = NULL, ...
) {
    block_data <- .hicpotts_block_data(fit, data, require = TRUE)
    if (is.null(iterations)) {
        iterations <- vapply(fit$blocks, function(block) {
            chains <- if (inherits(block, "hicpotts_robust_fit")) {
                block$fits
            } else {
                list(block)
            }
            max(vapply(chains, function(x) nrow(x$chains[[1L]]), integer(1)))
        }, integer(1))
    }
    if (length(iterations) == 1L) {
        iterations <- rep(iterations, length(fit$blocks))
    }
    if (length(iterations) != length(fit$blocks)) {
        stop("'iterations' must be scalar or have one value per block.",
            call. = FALSE
        )
    }
    dots <- list(...)
    parts <- lapply(seq_along(fit$blocks), function(j) {
        args <- c(list(
            data = block_data[[j]], chain_betas = fit$blocks[[j]],
            iterations = iterations[j]
        ), dots)
        if (is.null(args$N)) args$N <- fit$block_table$rows[j]
        do.call(compute_HMRFHiC_probabilities, args)
    })
    .hicpotts_bind_block_tables(parts, fit)
}

.hicpotts_block_xvars <- function(fit, x_vars = NULL) {
    n_blocks <- length(fit$blocks)
    if (is.null(x_vars)) {
        if (is.null(fit$processed)) return(NULL)
        return(lapply(seq_len(n_blocks), function(j) {
            .hicpotts_processed_xvars(fit$processed, j)
        }))
    }
    if (.hicpotts_is_processed(x_vars)) {
        checked <- .hicpotts_validate_processed(
            x_vars, n_blocks,
            lapply(fit$blocks, .hicpotts_fit_dimension)
        )
        return(lapply(seq_len(n_blocks), function(j) {
            .hicpotts_processed_xvars(checked, j)
        }))
    }
    required <- c("distance", "GC", "TES", "ACC")
    is_one <- is.list(x_vars) && all(required %in% names(x_vars))
    if (is_one && n_blocks == 1L) return(list(x_vars))
    is_many <- is.list(x_vars) && length(x_vars) == n_blocks &&
        all(vapply(x_vars, function(value) {
            is.list(value) && all(required %in% names(value))
        }, logical(1)))
    if (!is_many) {
        stop(
            "'x_vars' must be structured process_data() output or contain ",
            "one covariate list per block.",
            call. = FALSE
        )
    }
    x_vars
}

.hicpotts_split_block_result <- function(value, fit, name) {
    if (is.null(value)) return(rep(list(NULL), length(fit$blocks)))
    if (is.data.frame(value) && "block" %in% names(value)) {
        return(lapply(fit$block_table$block, function(label) {
            value[value$block == label, , drop = FALSE]
        }))
    }
    if (is.list(value) && length(value) == length(fit$blocks)) return(value)
    stop(sprintf(
        "'%s' must contain one value per block or a data frame with 'block'.",
        name
    ), call. = FALSE)
}

.hicpotts_allocation_diagnostics_blocks <- function(
    fit, component_names
) {
    blocks <- lapply(fit$blocks, allocation_diagnostics,
        component_names = component_names
    )
    names(blocks) <- fit$block_table$block
    summary <- do.call(rbind, lapply(seq_along(blocks), function(j) {
        value <- blocks[[j]]
        data.frame(
            block = fit$block_table$block[j],
            selected_chains = fit$block_table$selected_chains[j],
            worst_cell_ess = value$worst_cell_ess,
            max_mcse = value$max_mcse,
            mean_between_chain_disagreement = if (is.null(
                value$between_chain_disagreement
            )) {
                NA_real_
            } else {
                value$between_chain_disagreement$mean
            },
            stringsAsFactors = FALSE
        )
    }))
    structure(list(blocks = blocks, summary = summary),
        class = "hicpotts_block_diagnostics"
    )
}

.hicpotts_diagnose_blocks <- function(
    fit, burnin, ci_level, prob_result, minimum_component_cells,
    minimum_ess, maximum_rhat, minimum_chains,
    gamma_boundary_tolerance, maximum_gamma_boundary_fraction,
    minimum_gamma_unique, covariate_diagnostics, relabel
) {
    probabilities <- .hicpotts_split_block_result(
        prob_result, fit, "prob_result"
    )
    covariates <- .hicpotts_split_block_result(
        covariate_diagnostics, fit, "covariate_diagnostics"
    )
    blocks <- lapply(seq_along(fit$blocks), function(j) {
        diagnose_hicpotts_fit(
            fit$blocks[[j]], burnin = burnin, ci_level = ci_level,
            prob_result = probabilities[[j]],
            minimum_component_cells = minimum_component_cells,
            minimum_ess = minimum_ess, maximum_rhat = maximum_rhat,
            minimum_chains = minimum_chains,
            gamma_boundary_tolerance = gamma_boundary_tolerance,
            maximum_gamma_boundary_fraction =
                maximum_gamma_boundary_fraction,
            minimum_gamma_unique = minimum_gamma_unique,
            covariate_diagnostics = covariates[[j]], relabel = relabel
        )
    })
    names(blocks) <- fit$block_table$block
    summary <- do.call(rbind, lapply(seq_along(blocks), function(j) {
        flags <- blocks[[j]]$reliability_flags
        failed <- flags$criterion[!flags$passed]
        data.frame(
            block = fit$block_table$block[j],
            selected_chains = fit$block_table$selected_chains[j],
            reliable = !length(failed),
            failed_gates = paste(failed, collapse = ", "),
            stringsAsFactors = FALSE
        )
    }))
    structure(list(blocks = blocks, summary = summary),
        class = "hicpotts_block_fit_diagnostics"
    )
}

.hicpotts_summarise_parameter_blocks <- function(
    fit, require_reliable, x_vars, diagnose_covariates,
    covariate_diagnostic_args, pool_blocks, pooling_draws, pooling_seed,
    pooling_min_expected_cells, ...
) {
    block_x <- .hicpotts_block_xvars(fit, x_vars)
    parts <- lapply(seq_along(fit$blocks), function(j) {
        value <- summarise_hicpotts_parameters(
            fit$blocks[[j]], require_reliable = require_reliable,
            x_vars = if (is.null(block_x)) NULL else block_x[[j]],
            diagnose_covariates = diagnose_covariates,
            covariate_diagnostic_args = covariate_diagnostic_args, ...
        )
        value <- as.data.frame(value)
        value$block <- fit$block_table$block[j]
        value$selected_chains <- fit$block_table$selected_chains[j]
        value
    })
    out <- do.call(rbind, parts)
    rownames(out) <- NULL
    attr(out, "blocks") <- fit$block_table
    class(out) <- c("hicpotts_block_parameter_summary", "data.frame")
    if (identical(pool_blocks, "none")) return(out)
    .hicpotts_pool_parameter_blocks(
        fit = fit, block_summaries = out, pooling_draws = pooling_draws,
        pooling_seed = pooling_seed,
        pooling_min_expected_cells = pooling_min_expected_cells, ...
    )
}

.hicpotts_selected_block_fits <- function(fit) {
    if (inherits(fit, "hicpotts_robust_fit")) fit$fits else list(fit)
}

.hicpotts_block_pool_burnin <- function(fit, burnin) {
    if (!is.null(burnin)) return(burnin)
    if (inherits(fit, "hicpotts_robust_fit")) {
        value <- fit$settings$burnin
        if (length(value) == 1L && is.finite(value)) return(as.integer(value))
    }
    NULL
}

.hicpotts_expected_component_cells <- function(fit) {
    chains <- .hicpotts_selected_block_fits(fit)
    values <- lapply(chains, function(chain) {
        probabilities <- chain$z_probabilities
        if (!is.list(probabilities) || length(probabilities) != 3L ||
            !all(vapply(probabilities, is.matrix, logical(1L)))) {
            return(NULL)
        }
        occupancy <- vapply(probabilities, function(x) {
            sum(as.numeric(x), na.rm = TRUE)
        }, numeric(1L))
        draws <- chain$z_probability_draws
        if (length(draws) != 1L || !is.finite(draws) || draws <= 0) draws <- 1
        list(occupancy = occupancy, draws = as.numeric(draws))
    })
    available <- !vapply(values, is.null, logical(1L))
    values <- values[available]
    if (!length(values)) return(rep(NA_real_, 3L))
    denominator <- sum(vapply(values, `[[`, numeric(1L), "draws"))
    numerator <- Reduce(`+`, lapply(values, function(x) {
        x$occupancy * x$draws
    }))
    numerator / denominator
}

.hicpotts_pool_parameter_blocks <- function(
    fit, block_summaries, pooling_draws, pooling_seed,
    pooling_min_expected_cells, burnin = NULL, ci_level = 0.95, ...
) {
    if (length(ci_level) != 1L || !is.finite(ci_level) ||
        ci_level <= 0 || ci_level >= 1) {
        stop("'ci_level' must be strictly between zero and one.")
    }
    n_blocks <- length(fit$blocks)
    block_draw_sets <- lapply(seq_len(n_blocks), function(j) {
        current_burnin <- .hicpotts_block_pool_burnin(
            fit$blocks[[j]], burnin
        )
        .hicpotts_draw_sets(
            .hicpotts_selected_block_fits(fit$blocks[[j]]),
            burnin = current_burnin
        )
    })
    occupancy <- do.call(rbind, lapply(
        fit$blocks, .hicpotts_expected_component_cells
    ))
    colnames(occupancy) <- paste0("component", seq_len(3L))
    rownames(occupancy) <- fit$block_table$block

    parameters <- unique(block_summaries$parameter)
    alpha <- (1 - ci_level) / 2

    pooled <- withr::with_seed(pooling_seed, do.call(rbind, lapply(
        parameters, function(parameter) {
            component_text <- sub(
                "^component([123]):.*$", "\\1", parameter
            )
            component <- if (grepl("^component[123]:", parameter)) {
                as.integer(component_text)
            } else {
                NA_integer_
            }
            raw_weights <- if (is.na(component)) {
                as.numeric(fit$block_table$cells)
            } else {
                occupancy[, component]
            }
            weighting_method <- if (is.na(component)) {
                "analysed_cells"
            } else {
                "posterior_expected_component_cells"
            }

            draws <- lapply(block_draw_sets, function(sets) {
                value <- sets[[parameter]]
                if (is.null(value)) return(numeric())
                value <- unlist(value, use.names = FALSE)
                value[is.finite(value)]
            })
            available <- vapply(draws, length, integer(1L)) > 0L &
                is.finite(raw_weights) & raw_weights > 0
            if (!is.na(component)) {
                available <- available &
                    raw_weights >= pooling_min_expected_cells
            }
            excluded <- fit$block_table$block[!available]
            if (!any(available)) {
                return(data.frame(
                    parameter = parameter, estimate = NA_real_,
                    posterior_sd = NA_real_, CI_lower = NA_real_,
                    CI_upper = NA_real_, CI_width = NA_real_,
                    probability_positive = NA_real_,
                    probability_negative = NA_real_, ESS = NA_real_,
                    bulk_ESS = NA_real_, tail_ESS = NA_real_, Rhat = NA_real_,
                    resolved = FALSE, weighting_method = weighting_method,
                    effective_cells = 0, n_blocks_used = 0L,
                    blocks_used = "", blocks_excluded = paste(
                        excluded, collapse = ","
                    ), block_estimate_min = NA_real_,
                    block_estimate_max = NA_real_,
                    between_block_sd = NA_real_, stringsAsFactors = FALSE
                ))
            }

            used <- which(available)
            weights <- raw_weights[used] / sum(raw_weights[used])
            used_draws <- draws[used]
            sampled <- vapply(seq_along(used_draws), function(j) {
                index <- sample.int(
                    length(used_draws[[j]]), pooling_draws, replace = TRUE
                )
                used_draws[[j]][index] * weights[j]
            }, numeric(pooling_draws))
            if (is.null(dim(sampled))) sampled <- matrix(sampled, ncol = 1L)
            pooled_draws <- rowSums(sampled)
            interval <- stats::quantile(
                pooled_draws, c(alpha, 1 - alpha), names = FALSE,
                na.rm = TRUE
            )
            block_means <- vapply(used_draws, mean, numeric(1L))
            block_variances <- vapply(used_draws, stats::var, numeric(1L))
            block_variances[!is.finite(block_variances)] <- 0
            estimate <- sum(weights * block_means)
            between_sd <- sqrt(sum(weights * (block_means - estimate)^2))

            parameter_rows <- lapply(used, function(j) {
                block_summaries[
                    block_summaries$block == fit$block_table$block[j] &
                        block_summaries$parameter == parameter,
                    , drop = FALSE
                ]
            })
            resolution <- vapply(parameter_rows, function(x) {
                if (!nrow(x) || !"resolved" %in% names(x)) return(NA)
                isTRUE(x$resolved[1L])
            }, logical(1L))
            resolved <- if (anyNA(resolution)) NA else all(resolution)

            data.frame(
                parameter = parameter, estimate = estimate,
                posterior_sd = sqrt(sum(
                    weights^2 * block_variances
                )),
                CI_lower = interval[1L], CI_upper = interval[2L],
                CI_width = interval[2L] - interval[1L],
                probability_positive = mean(pooled_draws > 0),
                probability_negative = mean(pooled_draws < 0),
                ESS = NA_real_, bulk_ESS = NA_real_, tail_ESS = NA_real_,
                Rhat = NA_real_, resolved = resolved,
                weighting_method = weighting_method,
                effective_cells = sum(raw_weights[used]),
                n_blocks_used = length(used),
                blocks_used = paste(
                    fit$block_table$block[used], collapse = ","
                ),
                blocks_excluded = paste(excluded, collapse = ","),
                block_estimate_min = min(block_means),
                block_estimate_max = max(block_means),
                between_block_sd = between_sd,
                stringsAsFactors = FALSE
            )
        }
    )))
    rownames(pooled) <- NULL
    attr(pooled, "block_summaries") <- block_summaries
    attr(pooled, "block_table") <- fit$block_table
    attr(pooled, "posterior_expected_component_cells") <- occupancy
    attr(pooled, "pooling") <- list(
        estimand = paste0(
            "descriptive weighted average of independently fitted block ",
            "parameters; not a joint shared-parameter posterior"
        ),
        component_parameter_weights =
            "posterior expected component occupancy",
        global_parameter_weights = "analysed block cells",
        credible_interval = ci_level,
        monte_carlo_draws = pooling_draws,
        seed = pooling_seed,
        minimum_expected_component_cells = pooling_min_expected_cells,
        exclusions = stats::setNames(
            strsplit(pooled$blocks_excluded, ",", fixed = TRUE),
            pooled$parameter
        )
    )
    class(pooled) <- c(
        "hicpotts_pooled_parameter_summary", "hicpotts_parameter_summary",
        "data.frame"
    )
    pooled
}

.hicpotts_posterior_predictive_blocks <- function(
    fit, x_vars, y, dist, burnin, n_rep, seed
) {
    block_x <- .hicpotts_block_xvars(fit, x_vars)
    if (is.null(y)) {
        if (is.null(fit$processed)) {
            stop(
                "Supply 'x_vars' and 'y', or store process_data() output ",
                "when combining the block fits.",
                call. = FALSE
            )
        }
        y <- fit$processed$y
    } else if (is.matrix(y) && length(fit$blocks) == 1L) {
        y <- list(y)
    }
    if (is.null(block_x) || !is.list(y) ||
        length(y) != length(fit$blocks)) {
        stop("'x_vars' and 'y' must contain one entry per fitted block.",
            call. = FALSE
        )
    }
    blocks <- lapply(seq_along(fit$blocks), function(j) {
        posterior_predictive_hicpotts(
            fit$blocks[[j]], x_vars = block_x[[j]], y = y[[j]],
            dist = dist, burnin = burnin, n_rep = n_rep,
            seed = as.integer(seed) + j - 1L
        )
    })
    names(blocks) <- fit$block_table$block
    summary <- do.call(rbind, lapply(seq_along(blocks), function(j) {
        value <- as.data.frame(blocks[[j]]$summary)
        value$block <- fit$block_table$block[j]
        value
    }))
    discrepancies <- do.call(rbind, lapply(seq_along(blocks), function(j) {
        value <- as.data.frame(blocks[[j]]$discrepancies)
        value$block <- fit$block_table$block[j]
        value
    }))
    rownames(summary) <- rownames(discrepancies) <- NULL
    structure(
        list(
            blocks = blocks, summary = summary,
            discrepancies = discrepancies,
            block_table = fit$block_table
        ),
        class = "hicpotts_block_posterior_predictive"
    )
}

.hicpotts_check_block_overlap <- function(combined) {
    coords <- c("start", "start.j.")
    if (!all(coords %in% names(combined))) return(invisible(NULL))
    key <- paste(combined$start, combined$start.j., sep = "\r")
    duplicate <- duplicated(key)
    if (any(duplicate)) {
        stop(sprintf(
            "%d bin pair(s) occur in more than one block; blocks overlap.",
            sum(duplicate)
        ), call. = FALSE)
    }
    invisible(NULL)
}

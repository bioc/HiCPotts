#' Construct covariate-scale-aware fixed priors
#'
#' Builds the existing `user_fixed_priors` structure so each slope has the same
#' prior SD on a standardised `log1p` covariate scale. The returned priors and
#' all fitted coefficients remain on the original manuscript scale.
#'
#' @param x_vars Named HiCPotts covariate list.
#' @param dataset_index Dataset index for lists containing several matrices.
#' @param slope_sd_standardized Positive prior SD for a one-SD covariate effect.
#' @param intercept_mean Scalar or length-three intercept prior mean.
#' @param intercept_sd Scalar or length-three positive intercept prior SD.
#' @param minimum_covariate_sd Lower bound for nearly constant covariates.
#' @return A list accepted by `user_fixed_priors`.
#' @examples
#' N <- 6
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' priors <- make_hicpotts_scaled_priors(x_vars, slope_sd_standardized = 0.5)
#' str(priors$component1)
#' @noRd
make_hicpotts_scaled_priors <- function(
    x_vars, dataset_index = 1L, slope_sd_standardized = 0.5,
    intercept_mean = 0, intercept_sd = 5, minimum_covariate_sd = 1e-8
) {
    required <- c("distance", "GC", "TES", "ACC")
    if (!is.list(x_vars) || !all(required %in% names(x_vars))) {
        stop("'x_vars' must contain distance, GC, TES and ACC.")
    }
    dataset_index <- as.integer(dataset_index)
    if (length(dataset_index) != 1L || is.na(
        dataset_index
    ) || dataset_index < 1L) {
        stop("'dataset_index' must be a positive integer.")
    }
    if (length(slope_sd_standardized) != 1L ||
        !is.finite(slope_sd_standardized) || slope_sd_standardized <= 0) {
        stop("'slope_sd_standardized' must be finite and positive.")
    }
    if (length(minimum_covariate_sd) != 1L ||
        !is.finite(minimum_covariate_sd) || minimum_covariate_sd <= 0) {
        stop("'minimum_covariate_sd' must be finite and positive.")
    }
    intercept_mean <- rep(intercept_mean, length.out = 3L)
    intercept_sd <- rep(intercept_sd, length.out = 3L)
    if (any(!is.finite(intercept_mean)) || any(!is.finite(intercept_sd)) || any(
        intercept_sd <= 0
    )) {
        stop("Intercept prior means must be finite and SDs must be positive.")
    }

    cov_sd <- vapply(required, function(nm) {
        item <- x_vars[[nm]]
        mat <- if (is.list(item)) item[[dataset_index]] else item
        if (!is.matrix(mat) || !is.numeric(mat) || any(!is.finite(mat)) || any(
            mat <= -1
        )) {
            stop(
                "Each selected covariate must be a finite numeric matrix ",
                "greater than -1."
            )
        }
        value <- stats::sd(log1p(as.numeric(mat)))
        if (!is.finite(
            value
        ) || value < minimum_covariate_sd) {
            minimum_covariate_sd
        } else {
            value
        }
    }, numeric(1))
    slope_sd <- slope_sd_standardized / cov_sd
    out <- lapply(seq_len(3L), function(k) {
        list(
            meany = intercept_mean[k], meanx1 = 0, meanx2 = 0, meanx3 = 0,
            meanx4 = 0,
            sdy = intercept_sd[k], sdx1 = unname(slope_sd[1]), sdx2 = unname(
                slope_sd[2]
            ),
            sdx3 = unname(slope_sd[3]), sdx4 = unname(slope_sd[4])
        )
    })
    names(out) <- paste0("component", seq_len(3L))
    attr(out, "log1p_covariate_sd") <- stats::setNames(cov_sd, required)
    attr(out, "slope_sd_standardized") <- slope_sd_standardized
    out
}

.hicpotts_ess <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    n <- length(x)
    ## A continuous MCMC parameter whose chain never moves has not been sampled
    ## perfectly - it is frozen (stuck proposal, degenerate scale, numerical
    ## failure). Previously this branch returned the full chain length, so a
    ## constant chain was awarded the maximum possible ESS. Report NA
    ## (non-diagnostic); the reliability gate treats NA as a failure.
    if (n < 4L) {
        return(NA_real_)
    }
    if (stats::var(x) == 0) {
        return(NA_real_)
    }
    rho <- as.numeric(stats::acf(x,
        lag.max = min(n - 1L, floor(n / 2)),
        plot = FALSE, demean = TRUE
    )$acf)[-1L]
    if (!length(rho)) {
        return(as.numeric(n))
    }
    odd <- seq.int(1L, length(rho), by = 2L)
    even <- seq.int(2L, length(rho), by = 2L)
    pairs <- rho[odd] + c(rho[even], 0)[seq_along(odd)]
    first_nonpositive <- which(pairs <= 0)
    keep <- if (length(first_nonpositive)) {
        first_nonpositive[
            1L
        ] - 1L
    } else {
        length(pairs)
    }
    if (keep <= 0L) {
        return(as.numeric(n))
    }
    tau <- 1 + 2 * sum(rho[seq_len(min(length(rho), 2L * keep))])
    max(1, min(n, n / max(tau, 1)))
}

## ---------------------------------------------------------------------------
## Rank-normalized split-chain convergence diagnostics (Vehtari et al. 2021).
##
## The plain split R-hat and autocorrelation ESS above assume roughly normal,
## finite-variance marginals. HiCPotts chains routinely violate that: gamma and
## theta live on bounded intervals and pile up near their edges, and the
## component intercepts are heavy-tailed while the crossing barrier is active.
## On such chains plain R-hat is optimistic -- it can sit below 1.01 while the
## chains disagree about the tails -- and a single ESS number hides the fact
## that tail quantiles are far less well resolved than the mean.
##
## Rank normalisation makes the diagnostic invariant to monotone
## transformations and removes the normality assumption. Bulk ESS (on
## rank-normalized draws) describes the centre; tail ESS (on the folded,
## rank-normalized draws) describes the 5%/95% region, which is what credible
## intervals actually depend on.
## ---------------------------------------------------------------------------

## Split each chain in half, giving 2C sequences of equal length.
.hicpotts_split_chains <- function(chains) {
    chains <- lapply(chains, function(x) as.numeric(x[is.finite(x)]))
    n <- min(vapply(chains, length, integer(1)))
    half <- floor(n / 2L)
    if (half < 4L) {
        return(NULL)
    }
    unlist(lapply(chains, function(x) {
        x <- utils::tail(x, 2L * half)
        list(x[seq_len(half)], x[half + seq_len(half)])
    }), recursive = FALSE)
}

## Rank-normalize pooled draws to approximate normality: average ranks, then
## the inverse normal CDF of the Blom-transformed ranks.
.hicpotts_rank_normalize <- function(splits) {
    pooled <- unlist(splits, use.names = FALSE)
    r <- rank(pooled, ties.method = "average")
    S <- length(pooled)
    z <- stats::qnorm((r - 3 / 8) / (S - 2 * 3 / 8 + 1))
    lens <- vapply(splits, length, integer(1))
    split(z, rep(seq_along(splits), lens))
}

## Fold about the median: |x - median(x)|. Applying the rank normalisation to
## the folded draws yields a diagnostic sensitive to scale/tail behaviour.
.hicpotts_fold <- function(splits) {
    med <- stats::median(unlist(splits, use.names = FALSE))
    lapply(splits, function(x) abs(x - med))
}

## Classic split R-hat on already-prepared equal-length sequences.
.hicpotts_rhat_core <- function(splits) {
    m <- length(splits)
    if (m < 2L) {
        return(NA_real_)
    }
    half <- length(splits[[1L]])
    mat <- vapply(splits, identity, numeric(half))
    w <- mean(apply(mat, 2L, stats::var))
    ## Zero within-chain variance means every split half is constant: a frozen
    ## sampler, not a converged one. NA (non-diagnostic) rather than 1.
    if (!is.finite(w) || w <= 0) {
        return(NA_real_)
    }
    b <- half * stats::var(colMeans(mat))
    sqrt((((half - 1) / half) * w + b / half) / w)
}

## Multi-chain ESS via the Geyer initial-positive-sequence rule.
.hicpotts_ess_core <- function(splits) {
    m <- length(splits)
    if (m < 1L) {
        return(NA_real_)
    }
    n <- length(splits[[1L]])
    if (n < 4L) {
        return(NA_real_)
    }
    mat <- vapply(splits, identity, numeric(n))
    w <- mean(apply(mat, 2L, stats::var))
    if (!is.finite(w) || w <= 0) {
        return(NA_real_)
    }
    max_lag <- min(n - 1L, floor(n / 2))
    acf_per_chain <- lapply(seq_len(m), function(k) {
        as.numeric(stats::acf(mat[, k],
            lag.max = max_lag, plot = FALSE,
            demean = TRUE
        )$acf)[-1L]
    })
    rho <- Reduce(`+`, acf_per_chain) / m
    if (!length(rho)) {
        return(as.numeric(n * m))
    }
    odd <- seq.int(1L, length(rho), by = 2L)
    even <- seq.int(2L, length(rho), by = 2L)
    pairs <- rho[odd] + c(rho[even], 0)[seq_along(odd)]
    first_nonpositive <- which(pairs <= 0)
    keep <- if (length(first_nonpositive)) {
        first_nonpositive[
            1L
        ] - 1L
    } else {
        length(pairs)
    }
    if (keep <= 0L) {
        return(as.numeric(n * m))
    }
    tau <- 1 + 2 * sum(rho[seq_len(min(length(rho), 2L * keep))])
    max(1, min(n * m, (n * m) / max(tau, 1)))
}

## Rank-normalized split R-hat.
#' @noRd
.hicpotts_rhat <- function(chains) {
    splits <- .hicpotts_split_chains(chains)
    if (is.null(splits) || length(splits) < 2L) {
        return(NA_real_)
    }
    if (stats::var(unlist(splits, use.names = FALSE)) <= 0) {
        return(NA_real_)
    }
    bulk <- .hicpotts_rhat_core(.hicpotts_rank_normalize(splits))
    tail <- .hicpotts_rhat_core(.hicpotts_rank_normalize(.hicpotts_fold(
        splits
    )))
    ## Report the more pessimistic of the bulk and tail diagnostics, so a chain
    ## that agrees on the centre but not on the tails is not certified.
    if (!is.finite(bulk) && !is.finite(tail)) {
        return(NA_real_)
    }
    max(c(bulk, tail), na.rm = TRUE)
}

## Bulk effective sample size (rank-normalized).
#' @noRd
.hicpotts_bulk_ess <- function(chains) {
    splits <- .hicpotts_split_chains(chains)
    if (is.null(splits)) {
        return(NA_real_)
    }
    if (stats::var(unlist(splits, use.names = FALSE)) <= 0) {
        return(NA_real_)
    }
    .hicpotts_ess_core(.hicpotts_rank_normalize(splits))
}

## Tail effective sample size: the minimum of the ESS of the 5% and 95%
## indicator series, which is what the credible-interval endpoints rely on.
#' @noRd
.hicpotts_tail_ess <- function(chains) {
    splits <- .hicpotts_split_chains(chains)
    if (is.null(splits)) {
        return(NA_real_)
    }
    pooled <- unlist(splits, use.names = FALSE)
    if (stats::var(pooled) <= 0) {
        return(NA_real_)
    }
    q <- stats::quantile(pooled, c(0.05, 0.95), names = FALSE, na.rm = TRUE)
    ess_indicator <- function(threshold, upper) {
        ind <- lapply(splits, function(x) {
            as.numeric(if (upper) x <= threshold else x >= threshold)
        })
        if (stats::var(unlist(ind, use.names = FALSE)) <= 0) {
            return(NA_real_)
        }
        .hicpotts_ess_core(ind)
    }
    vals <- c(ess_indicator(q[1L], FALSE), ess_indicator(q[2L], TRUE))
    if (all(!is.finite(vals))) {
        return(NA_real_)
    }
    min(vals, na.rm = TRUE)
}

.hicpotts_fit_list <- function(fit) {
    if (is.list(fit) && is.list(fit$chains) && length(fit$chains) == 3L) {
        return(list(fit))
    }
    valid <- is.list(fit) && length(fit) &&
        all(vapply(fit, function(x) {
            is.list(x) && is.list(x$chains) &&
                length(x$chains) == 3L
        }, logical(1)))
    if (!valid) {
        stop(
            "'fit' must be a HiCPotts fit or a list of repeated-chain fits."
        )
    }
    fit
}

## Zero inflation exists only for ZIP and ZINB. Under Poisson or NB
## the sampler carries theta at its starting value and never updates
## it, so the draw sequence is constant and its ESS and R-hat are
## undefined. Prefer the family recorded on the fit; fall back to
## whether theta actually moved for fits produced before the sampler
## settings were stored.
.hicpotts_has_zero_inflation <- function(fit) {
    family <- fit$sampler_settings$distribution
    if (is.character(family) && length(family) == 1L &&
        !is.na(family)) {
        return(family %in% c("ZIP", "ZINB"))
    }
    theta <- fit$theta
    if (!is.numeric(theta) || length(theta) < 2L) {
        return(FALSE)
    }
    length(unique(theta[is.finite(theta)])) > 1L
}

.hicpotts_draw_sets <- function(fits, burnin) {
    beta_names <- c("intercept", "distance", "GC", "TES", "ACC")
    sets <- list()
    for (k in seq_len(3)) {
        for (j in seq_len(5)) {
            sets[[paste0(
                "component", k,
                ":", beta_names[j]
            )]] <- list()
        }
    }
    sets$gamma <- list()
    theta_sets <- list()
    size_sets <- list(
        component1 = list(), component2 = list(),
        component3 = list()
    )
    for (m in seq_along(fits)) {
        f <- fits[[m]]
        n_iter <- nrow(f$chains[[1]])
        b <- if (is.null(burnin)) {
            as.integer(floor(n_iter / 2))
        } else {
            as.integer(
                burnin
            )
        }
        if (length(b) != 1L || is.na(b) || b < 0L || b >= n_iter) {
            stop(
                "'burnin' must be between 0 and the number of iterations ",
                "minus one."
            )
        }
        idx <- seq.int(b + 1L, n_iter)
        for (k in seq_len(3)) {
            for (j in seq_len(5)) {
                sets[[paste0("component", k, ":", beta_names[j])]][[
                    m
                ]] <- f$chains[[k]][idx, j]
            }
        }
        sets$gamma[[m]] <- f$gamma[idx]
        if (.hicpotts_has_zero_inflation(f)) {
            theta_sets[[m]] <- f$theta[idx]
        }
        if (is.matrix(f$size) && nrow(f$size) == 3L && ncol(f$size) >= max(
            idx
        ) &&
            any(f$size[, idx, drop = FALSE] > 0, na.rm = TRUE)) {
            for (k in seq_len(3)) size_sets[[k]][[m]] <- f$size[k, idx]
        }
    }
    ## Only monitor theta when every fit actually sampled it.
    if (length(theta_sets) == length(fits) &&
        all(vapply(theta_sets, length, integer(1)) > 0L)) {
        sets$theta <- theta_sets
    }
    for (k in seq_len(3)) {
        if (length(size_sets[[k]])) {
            sets[[paste0(
                "component",
                k, ":size"
            )]] <- size_sets[[k]]
        }
    }
    sets
}

#' Diagnose HiCPotts parameter estimation
#'
#' Reports posterior estimates, intervals, sign probabilities, ESS and R-hat,
#' plus internal parameter occupancy and reported classification occupancy.
#' Optional probability output adds mean classification entropy.
#'
#' @details
#' The returned \code{reliability_flags} table holds one row per diagnostic
#' criterion, with its configured threshold and observed result.
#' Parameter summaries and cell-level classifications have separate
#' diagnostic outputs. Inspect the criteria relevant to the quantity being
#' reported, together with its posterior interval and Monte Carlo precision.
#' \describe{
#'   \item{\code{mcse_precision}}{The batch-means Monte Carlo standard
#'     error reached its relative threshold before the requested
#'     iterations were exhausted.}
#'   \item{\code{beta_movement}}{Every component accepted at least one
#'     ordinary retained regression-coefficient proposal. Branch-label swaps
#'     are counted separately and cannot satisfy this gate.}
#'   \item{\code{common_posterior_target}}{All production chains use the same
#'     frozen regression-prior hyperparameters, as required for interpretable
#'     multi-chain R-hat and ESS diagnostics.}
#'   \item{\code{gamma_not_boundary}}{In every chain, less than
#'     \code{maximum_gamma_boundary_fraction} of the retained gamma
#'     draws lie within \code{gamma_boundary_tolerance} of zero or one.
#'     This records posterior concentration near the limits of the support.}
#'   \item{\code{gamma_movement}}{Every chain retained at least
#'     \code{minimum_gamma_unique} distinct gamma values. This catches a
#'     frozen gamma chain even when the boundary test passes.}
#'   \item{\code{component_occupancy}}{Every component was allocated at
#'     least \code{minimum_component_cells} lattice cells internally.
#'     Occupancy describes the data contributing to each component.}
#'   \item{\code{effective_sample_size}}{Bulk and tail ESS both reach
#'     \code{minimum_ess} for every parameter. A non-finite ESS marks a
#'     stuck chain and fails the gate rather than being dropped. The tail
#'     is checked separately because a chain can mix well in the bulk
#'     with different precision in its outer quantiles.}
#'   \item{\code{split_rhat}}{Split-chain R-hat is at most
#'     \code{maximum_rhat} for every parameter. This needs at least two
#'     chains; with a single chain the statistic is undefined and the
#'     gate fails.}
#'   \item{\code{independent_chains}}{At least \code{minimum_chains}
#'     independently seeded chains were supplied.}
#'   \item{\code{relabelled}}{Every fit carries a relabelling rule, so
#'     components 2 and 3 mean the same thing in each chain before
#'     chains are pooled.}
#'   \item{\code{coefficient_scale}}{All chains report coefficients on
#'     the original log1p-covariate scale, so estimates from different
#'     chains are comparable.}
#'   \item{\code{covariate_conditioning}}{Carried through from
#'     \code{covariate_diagnostics} when supplied; it flags collinear or
#'     ill-conditioned covariates.}
#' }
#'
#' \code{summarise_hicpotts_parameters()} applies the subset of criteria
#' described on its help page. With \code{require_reliable = FALSE}, it
#' displays estimates together with their diagnostic indicators.
#' \code{classify_hicpotts()} provides the cell-level classification.
#'
#' @param fit A fit, list of independently seeded fits of the same data, or
#'   block-aware result from \code{combine_hicpotts_blocks()}.
#' @param burnin Initial rows to discard; defaults to half.
#' @param ci_level Posterior interval level.
#' @param prob_result Optional output from `compute_HMRFHiC_probabilities()`.
#' @param minimum_component_cells Non-negative occupancy threshold. The default
#'   zero permits empty components and retains occupancy counts for reporting.
#' @param minimum_ess Minimum recommended effective sample size.
#' @param maximum_rhat Maximum recommended split-chain R-hat.
#' @param minimum_chains Minimum recommended number of independent chains.
#' @param gamma_boundary_tolerance Distance from zero or one counted as the
#'   boundary region for the bounded Potts parameter.
#' @param maximum_gamma_boundary_fraction Largest allowed fraction of retained
#'   gamma draws in the boundary region in every chain.
#' @param minimum_gamma_unique Minimum number of distinct retained gamma draws
#'   required in every chain. A frozen gamma chain is rejected explicitly even
#'   when another diagnostic is unavailable.
#' @param covariate_diagnostics Optional pre-computed covariate
#'   conditioning report, used to fill the
#'   \code{covariate_conditioning} reliability gate. It is normally left
#'   as \code{NULL}. \code{\link{summarise_hicpotts_parameters}()} builds
#'   one internally and passes it through when it is given
#'   \code{x_vars} and \code{diagnose_covariates = TRUE} (the default);
#'   for a \code{hicpotts_robust_fit} without \code{x_vars} it reuses the
#'   report already stored on the fit. Supply a value here only to reuse a
#'   report you computed yourself. When \code{NULL} the covariate gate is
#'   recorded as passing, because no conditioning problem has been
#'   demonstrated -- it is not evidence that the covariates are well
#'   conditioned.
#' @param relabel Whether to relabel components 2 and 3 before pooling chains.
#' @return A list of diagnostics: \code{parameters} (estimates,
#'   intervals, ESS and R-hat), \code{reliability_flags} (one row per
#'   gate, with its threshold), \code{component_occupancy},
#'   \code{classification_entropy}, \code{mcse_status},
#'   \code{gamma_diagnostics}, \code{beta_mixing_status},
#'   \code{noise_relationship_probabilities}, the \code{fits} used and
#'   any \code{warnings}.
#' @examples
#' N <- 5
#' y <- matrix(rpois(N * N, 4), N, N)
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' fits <- run_chain_betas(
#'     N = N, iterations = 5, x_vars = x_vars, y = y,
#'     dist = "ZINB", theta_start = 0.5, size_start = c(1, 1, 1),
#'     n_chains = 2, seeds = 1:2, robust = FALSE
#' )
#' d <- diagnose_hicpotts_fit(
#'     fits,
#'     minimum_component_cells = 1, minimum_ess = 1
#' )
#' d$reliability_flags[, c("criterion", "passed")]
#' d$classification_entropy
#' @export
diagnose_hicpotts_fit <- function(
    fit, burnin = NULL, ci_level = 0.95,
    prob_result = NULL, minimum_component_cells = 0L,
    minimum_ess = 200, maximum_rhat = 1.01, minimum_chains = 4L,
    gamma_boundary_tolerance = 0.01,
    maximum_gamma_boundary_fraction = 0.95,
    minimum_gamma_unique = 20L,
    covariate_diagnostics = NULL, relabel = TRUE
) {
    if (inherits(fit, "hicpotts_block_fit")) {
        return(.hicpotts_diagnose_blocks(
            fit, burnin = burnin, ci_level = ci_level,
            prob_result = prob_result,
            minimum_component_cells = minimum_component_cells,
            minimum_ess = minimum_ess, maximum_rhat = maximum_rhat,
            minimum_chains = minimum_chains,
            gamma_boundary_tolerance = gamma_boundary_tolerance,
            maximum_gamma_boundary_fraction =
                maximum_gamma_boundary_fraction,
            minimum_gamma_unique = minimum_gamma_unique,
            covariate_diagnostics = covariate_diagnostics,
            relabel = relabel
        ))
    }
    if (inherits(fit, "hicpotts_robust_fit")) {
        if (is.null(covariate_diagnostics)) {
            covariate_diagnostics <- fit$covariate_diagnostics
        }
        fit <- fit$fits
    }
    fits <- .hicpotts_fit_list(fit)
    if (isTRUE(relabel)) fits <- relabel_hicpotts(fits)
    if (length(ci_level) != 1L || !is.finite(
        ci_level
    ) || ci_level <= 0 || ci_level >= 1) {
        stop("'ci_level' must be strictly between zero and one.")
    }
    minimum_component_cells <- as.integer(minimum_component_cells)
    if (length(minimum_component_cells) != 1L || is.na(
        minimum_component_cells
    ) ||
        minimum_component_cells < 0L) {
        stop("'minimum_component_cells' must be non-negative.")
    }
    minimum_chains <- as.integer(minimum_chains)
    if (!is.finite(minimum_ess) || minimum_ess <= 0) {
        stop("'minimum_ess' must be positive.")
    }
    if (!is.finite(maximum_rhat) || maximum_rhat <= 1) {
        stop("'maximum_rhat' must be greater than one.")
    }
    if (length(minimum_chains) != 1L || is.na(
        minimum_chains
    ) || minimum_chains < 2L) {
        stop("'minimum_chains' must be at least two.")
    }
    if (length(gamma_boundary_tolerance) != 1L ||
        !is.finite(gamma_boundary_tolerance) || gamma_boundary_tolerance <= 0 ||
        gamma_boundary_tolerance >= 0.5) {
        stop(
            "'gamma_boundary_tolerance' must be strictly between zero and 0.5."
        )
    }
    if (length(maximum_gamma_boundary_fraction) != 1L ||
        !is.finite(maximum_gamma_boundary_fraction) ||
        maximum_gamma_boundary_fraction <= 0 ||
        maximum_gamma_boundary_fraction > 1) {
        stop("'maximum_gamma_boundary_fraction' must be in (0, 1].")
    }
    minimum_gamma_unique <- as.integer(minimum_gamma_unique)
    if (length(minimum_gamma_unique) != 1L || is.na(minimum_gamma_unique) ||
        minimum_gamma_unique < 2L) {
        stop("'minimum_gamma_unique' must be an integer of at least two.")
    }
    sets <- .hicpotts_draw_sets(fits, burnin)
    alpha <- (1 - ci_level) / 2
    parameters <- do.call(rbind, lapply(names(sets), function(nm) {
        per_chain <- sets[[nm]]
        draws <- unlist(per_chain, use.names = FALSE)
        q <- stats::quantile(draws, c(alpha, 1 - alpha),
            names = FALSE,
            na.rm = TRUE
        )
        data.frame(
            parameter = nm, estimate = mean(draws, na.rm = TRUE),
            posterior_sd = stats::sd(draws, na.rm = TRUE),
            CI_lower = q[1], CI_upper = q[2], CI_width = q[2] - q[1],
            probability_positive = mean(draws > 0, na.rm = TRUE),
            probability_negative = mean(draws < 0, na.rm = TRUE),
            ## Rank-normalized diagnostics (Vehtari et al. 2021). ESS is
            ## retained as the bulk figure so existing callers keep working,
            ## with tail ESS reported alongside it: credible-interval endpoints
            ## depend on the tail, which the bulk figure can overstate
            ## substantially.
            ESS = .hicpotts_bulk_ess(per_chain),
            bulk_ESS = .hicpotts_bulk_ess(per_chain),
            tail_ESS = .hicpotts_tail_ess(per_chain),
            Rhat = .hicpotts_rhat(per_chain), stringsAsFactors = FALSE
        )
    }))
    rownames(parameters) <- NULL
    occupancy <- do.call(rbind, lapply(seq_along(fits), function(m) {
        z <- fits[[m]]$z_final
        internal <- attr(z, "parameter_component_counts", exact = TRUE)
        reported <- attr(z, "classification_component_counts", exact = TRUE)
        if (is.null(internal)) internal <- rep(NA_integer_, 3L)
        if (is.null(reported)) reported <- tabulate(as.integer(z), nbins = 3L)
        ## Posterior EXPECTED occupancy: the sum of each component's membership
        ## probability over all cells, i.e. the expected number of cells
        ## belonging to it under the posterior. The hard counts above come from
        ## one final allocation and so inherit that sweep's noise; a component
        ## sitting at 40% probability across many cells contributes to the
        ## expectation but may own no cell outright in any single sweep.
        zp <- fits[[m]]$z_probabilities
        expected_cells <- if (is.list(zp) && length(zp) == 3L) {
            vapply(zp, function(mm) sum(as.numeric(mm)), numeric(1))
        } else {
            rep(NA_real_, 3L)
        }
        data.frame(
            chain = m, component = seq_len(3), parameter_cells = as.integer(
                internal
            ),
            reported_classification_cells = as.integer(reported),
            posterior_expected_cells = as.numeric(expected_cells),
            sufficient_parameter_cells = as.integer(
                internal
            ) >= minimum_component_cells,
            relabel_basis = if (is.null(fits[[m]]$relabel_basis)) {
                NA_character_
            } else {
                fits[[m]]$relabel_basis
            }, stringsAsFactors = FALSE
        )
    }))
    ## Mean normalized entropy of the membership probabilities. When the
    ## caller supplies prob_result it is used; otherwise the sampler's own
    ## pooled latent-state frequencies are used, so the field is populated
    ## rather than left NA for every fit that did not pass prob_result.
    entropy <- NA_real_
    if (is.null(prob_result)) {
        stored <- lapply(fits, function(f) f$z_probabilities)
        usable <- vapply(stored, function(z) {
            is.list(z) && length(z) == 3L && all(vapply(z, is.matrix,
                logical(1)))
        }, logical(1))
        if (any(usable)) {
            keep <- stored[usable]
            pooled <- lapply(seq_len(3L), function(k) {
                Reduce(`+`, lapply(keep, `[[`, k)) / length(keep)
            })
            p <- cbind(as.vector(pooled[[1L]]), as.vector(pooled[[2L]]),
                as.vector(pooled[[3L]]))
            totals <- rowSums(p)
            p <- p[is.finite(totals) & totals > 0, , drop = FALSE]
            if (nrow(p)) {
                p <- p / rowSums(p)
                entropy <- mean(
                    -rowSums(ifelse(p > 0, p * log(p), 0)),
                    na.rm = TRUE
                )
            }
        }
    }
    if (!is.null(prob_result)) {
        pcols <- c("prob1", "prob2", "prob3")
        if (!is.data.frame(prob_result) || !all(pcols %in% names(
            prob_result
        ))) {
            stop("'prob_result' must contain prob1, prob2 and prob3.")
        }
        p <- as.matrix(prob_result[, pcols, drop = FALSE])
        p <- p / rowSums(p)
        entropy <- mean(-rowSums(ifelse(p > 0, p * log(p), 0)), na.rm = TRUE)
    }
    notes <- character()
    mcse_status <- do.call(rbind, lapply(seq_along(fits), function(m) {
        gamma <- fits[[m]]$gamma
        get_attribute <- function(name, default = NA) {
            value <- attr(gamma, name, exact = TRUE)
            if (is.null(value) || !length(value)) default else value[[1L]]
        }
        data.frame(
            chain = m,
            enabled = isTRUE(attr(gamma, "mcse_stopping_enabled",
                exact = TRUE
            )),
            converged = isTRUE(attr(gamma, "mcse_converged", exact = TRUE)),
            iterations_completed = as.integer(get_attribute(
                "iterations_completed"
            )),
            iterations_requested = as.integer(get_attribute(
                "iterations_requested"
            )),
            max_relative_mcse = as.numeric(get_attribute("mcse_max_relative")),
            threshold = as.numeric(get_attribute("mcse_relative_threshold"))
        )
    }))
    mcse_ok <- all(!mcse_status$enabled | mcse_status$converged)
    beta_mixing_status <- do.call(rbind, lapply(seq_along(fits), function(m) {
        mixing <- fits[[m]]$beta_mixing
        retained_attempts <- if (is.list(mixing)) {
            mixing$retained_attempts
        } else {
            NULL
        }
        retained_acceptances <- if (is.list(mixing)) {
            mixing$retained_acceptances
        } else {
            NULL
        }
        if (is.null(retained_attempts) || length(retained_attempts) != 3L) {
            retained_attempts <- rep(NA_real_, 3L)
        }
        if (is.null(retained_acceptances) ||
            length(retained_acceptances) != 3L) {
            retained_acceptances <- rep(NA_real_, 3L)
        }
        data.frame(
            chain = m, component = seq_len(3L),
            retained_attempts = as.numeric(retained_attempts),
            retained_acceptances = as.numeric(retained_acceptances),
            retained_acceptance_rate = ifelse(
                retained_attempts > 0,
                retained_acceptances / retained_attempts,
                NA_real_
            ),
            stringsAsFactors = FALSE
        )
    }))
    beta_movement_ok <- all(is.finite(
        beta_mixing_status$retained_acceptances
    )) && all(beta_mixing_status$retained_acceptances > 0)
    prior_objects <- lapply(fits, function(fit) {
        prior <- fit$regression_prior
        if (is.list(prior)) prior$frozen_hyperparameters else NULL
    })
    explicitly_shared <- all(vapply(fits, function(fit) {
        is.list(fit$regression_prior) &&
            isTRUE(fit$regression_prior$shared_across_chains)
    }, logical(1)))
    identical_fixed_priors <- length(prior_objects) > 0L &&
        !is.null(prior_objects[[1L]]) && all(vapply(
            prior_objects,
            function(value) identical(value, prior_objects[[1L]]),
            logical(1)
        ))
    common_target_ok <- explicitly_shared || identical_fixed_priors
    ## Gamma is bounded in [0,1]. A chain that spends nearly all retained draws
    ## within `gamma_boundary_tolerance` of either endpoint is reporting a
    ## truncated/boundary solution, not a well-resolved interior posterior. A
    ## chain with only a handful of unique values is separately rejected. A
    ## gamma chain can approach a boundary and then remain exactly constant
    ## throughout the retained draws. Generic ESS/R-hat already fail a constant
    ## series, but an explicit diagnostic makes this scientifically important
    ## failure impossible to overlook in reporting.
    gamma_diagnostics <- do.call(rbind, lapply(seq_along(sets$gamma), function(
        m
    ) {
        draws <- as.numeric(sets$gamma[[m]])
        draws <- draws[is.finite(draws)]
        at_boundary <- draws <= gamma_boundary_tolerance |
            draws >= 1 - gamma_boundary_tolerance
        data.frame(
            chain = m,
            retained_draws = length(draws),
            unique_draws = length(unique(signif(draws, 12L))),
            boundary_fraction = if (length(draws)) {
                mean(
                    at_boundary
                )
            } else {
                NA_real_
            },
            minimum = if (length(draws)) min(draws) else NA_real_,
            maximum = if (length(draws)) max(draws) else NA_real_,
            movement_fraction = if (length(draws) > 1L) {
                mean(diff(draws) != 0)
            } else {
                NA_real_
            },
            stringsAsFactors = FALSE
        )
    }))
    gamma_boundary_ok <- all(is.finite(gamma_diagnostics$boundary_fraction)) &&
        all(
            gamma_diagnostics$boundary_fraction <
                maximum_gamma_boundary_fraction
        )
    gamma_movement_ok <- all(
        gamma_diagnostics$unique_draws >= minimum_gamma_unique
    )
    occupancy_ok <- !any(occupancy$sufficient_parameter_cells %in% FALSE,
        na.rm = TRUE
    )
    ## NA ESS means a frozen/non-diagnostic chain (see .hicpotts_ess). It must
    ## FAIL the gate: with na.rm = TRUE it was previously dropped, so a stuck
    ## parameter passed silently.
    ess_ok <- all(is.finite(parameters$bulk_ESS)) &&
        all(parameters$bulk_ESS >= minimum_ess) &&
        ## Tail ESS measures precision in the outer quantiles separately
        ## from the bulk used for central posterior summaries.
        all(is.finite(parameters$tail_ESS)) &&
        all(parameters$tail_ESS >= minimum_ess)
    rhat_available <- length(fits) >= 2L && all(is.finite(parameters$Rhat))
    rhat_ok <- rhat_available && all(parameters$Rhat <= maximum_rhat)
    chain_count_ok <- length(fits) >= minimum_chains
    relabel_ok <- all(vapply(
        fits, function(x) !is.null(x$relabel_rule),
        logical(1)
    ))
    relationship_probabilities <- data.frame(
        chain = seq_along(fits),
        enabled = vapply(fits, function(x) {
            is.list(x$noise_relationship) &&
                isTRUE(x$noise_relationship$enabled)
        }, logical(1)),
        probability = vapply(fits, function(x) {
            value <- x$noise_relationship_probability
            if (is.null(value) || length(value) != 1L || !is.finite(value)) {
                NA_real_
            } else {
                as.numeric(value)
            }
        }, numeric(1)),
        threshold_policy = rep(
            "none in package; user chooses downstream probability threshold",
            length(fits)
        ),
        stringsAsFactors = FALSE
    )
    coefficient_scale <- unique(unlist(lapply(fits, function(x) {
        attr(x$chains[[1L]], "coefficient_scale", exact = TRUE)
    })))
    scale_ok <- length(coefficient_scale) == 1L &&
        identical(coefficient_scale, "original log1p-covariate scale")
    covariate_ok <- is.null(covariate_diagnostics) ||
        isTRUE(covariate_diagnostics$reliable)
    flags <- data.frame(
        criterion = c(
            "mcse_precision", "beta_movement", "common_posterior_target",
            "gamma_not_boundary",
            "gamma_movement",
            "component_occupancy", "effective_sample_size", "split_rhat",
            "independent_chains", "relabelled", "coefficient_scale",
            "covariate_conditioning"
        ),
        passed = c(
            mcse_ok, beta_movement_ok, common_target_ok,
            gamma_boundary_ok, gamma_movement_ok,
            occupancy_ok,
            ess_ok, rhat_ok, chain_count_ok, relabel_ok,
            scale_ok, covariate_ok
        ),
        threshold = c(
            "batch-means MCSE/SD threshold before requested iterations",
            ">0 retained ordinary beta acceptances per component and chain",
            "identical frozen regression prior in every production chain",
            paste0(
                "<", maximum_gamma_boundary_fraction, " within ",
                gamma_boundary_tolerance, " of 0 or 1"
            ),
            paste0(
                ">=", minimum_gamma_unique,
                " distinct retained draws per chain"
            ),
            paste0(">=", minimum_component_cells, " cells"),
            paste0(">=", minimum_ess), paste0("<=", maximum_rhat),
            paste0(">=", minimum_chains), "required",
            "original manuscript scale", "correlation/condition thresholds"
        ),
        stringsAsFactors = FALSE
    )
    if (!mcse_ok) {
        notes <- c(
            notes,
            paste0("At least one chain reached its iteration limit before ",
                "satisfying the MCSE precision threshold.")
        )
    }
    if (!beta_movement_ok) {
        notes <- c(
            notes,
            paste0(
                "At least one component accepted no ordinary retained beta ",
                "proposal; branch swaps do not establish coefficient mixing."
            )
        )
    }
    if (!common_target_ok) {
        notes <- c(
            notes,
            paste0(
                "Production chains do not record one common frozen ",
                "regression prior; pooled R-hat is not a same-target ",
                "diagnostic."
            )
        )
    }
    if (!gamma_boundary_ok) {
        notes <- c(notes, paste0(
            "Gamma's boundary fraction exceeds the configured threshold ",
            "in at least one chain. See gamma_diagnostics for the fraction ",
            "and retained-draw movement."
        ))
    }
    if (!gamma_movement_ok) {
        notes <- c(notes, paste0(
            "Gamma has too few distinct retained draws in at least one chain; ",
            "the ",
            "gamma transition is frozen or insufficiently explored."
        ))
    }
    if (!occupancy_ok) {
        notes <- c(
            notes,
            paste0("At least one component has too few internally allocated ",
                "cells for stable component-specific regression.")
        )
    }
    if (!ess_ok) {
        notes <- c(
            notes,
            paste0("At least one parameter has insufficient ESS; run longer ",
                "after checking mode agreement.")
        )
    }
    if (!rhat_ok) {
        notes <- c(
            notes,
            paste0("Split-chain R-hat is unavailable or above threshold; ",
                "investigate mixing and label stability.")
        )
    }
    if (!chain_count_ok) {
        notes <- c(
            notes,
            paste0("Fewer than the recommended number of independent chains ",
                "were supplied.")
        )
    }
    if (!relabel_ok) {
        notes <- c(
            notes,
            "At least one chain was not marked as relabelled before pooling."
        )
    }
    if (!scale_ok) {
        notes <- c(
            notes,
            "Coefficient-scale metadata is missing or inconsistent."
        )
    }
    if (!covariate_ok) {
        notes <- c(notes, covariate_diagnostics$warnings)
    }
    gamma_parameter <- parameters$parameter == "gamma"
    non_gamma_parameter <- !gamma_parameter
    finite_ess <- function(index) {
        length(index) && all(is.finite(parameters$bulk_ESS[index])) &&
            all(parameters$bulk_ESS[index] >= minimum_ess) &&
            all(is.finite(parameters$tail_ESS[index])) &&
            all(parameters$tail_ESS[index] >= minimum_ess)
    }
    finite_rhat <- function(index) {
        length(fits) >= 2L && length(index) &&
            all(is.finite(parameters$Rhat[index])) &&
            all(parameters$Rhat[index] <= maximum_rhat)
    }
    ## The per-gate table is the reported result. Summary verdicts are no
    ## longer returned: a single overall flag collapsed twelve independent
    ## criteria into one word, which read as a blanket failure whenever any
    ## one gate was marginal. Read reliability_flags and decide which gates
    ## matter for the claim being made.
    list(
        parameters = parameters, component_occupancy = occupancy,
        classification_entropy = entropy, reliability_flags = flags,
        coefficient_scale = coefficient_scale,
        mcse_status = mcse_status, gamma_diagnostics = gamma_diagnostics,
        beta_mixing_status = beta_mixing_status,
        noise_relationship_probabilities = relationship_probabilities,
        fits = fits, warnings = unique(notes)
    )
}

#' Report HiCPotts parameter summaries and diagnostic criteria
#'
#' Extracts posterior estimates, intervals and Monte Carlo diagnostics.
#' The reporting criteria check comparable posterior targets, coefficient
#' movement, component alignment, coefficient scales, covariate conditioning
#' and parameter-specific ESS and R-hat. Other recorded criteria provide
#' context for interpreting the fit. Classification is available separately
#' through \code{classify_hicpotts()}.
#'
#' @param fit A HiCPotts fit, list of fits, \code{hicpotts_robust_fit},
#'   block-aware result from \code{combine_hicpotts_blocks()}, or an object
#'   returned by \code{diagnose_hicpotts_fit()}.
#' @param require_reliable Logical. With the default \code{TRUE}, require
#'   a common frozen regression prior, beta movement, aligned components,
#'   compatible coefficient scales, covariate conditioning and the requested
#'   parameter-specific ESS and R-hat thresholds. Gamma boundary fraction,
#'   occupancy and recommended chain count are recorded separately in
#'   \code{reliability_flags}. With \code{FALSE}, return the table together
#'   with a \code{resolved} column indicating each row's ESS/R-hat result.
#' @param x_vars Optional named covariate list. When supplied and
#'   `diagnose_covariates = TRUE`, covariate conditioning is checked inside
#'   this summary call.
#' @param diagnose_covariates Logical. Run the internal covariate diagnostic
#'   when `x_vars` is supplied; otherwise reuse diagnostics stored in a robust
#'   fit when available.
#' @param covariate_diagnostic_args Optional named list of threshold arguments
#'   passed to the internal covariate diagnostic.
#' @param pool_blocks For a block-aware fit, `"none"` retains one row per
#'   block and parameter. `"posterior_weighted"` returns an overall descriptive
#'   posterior summary. Component-specific coefficients and dispersions are
#'   weighted by posterior expected component occupancy; `gamma` and `theta`
#'   are weighted by analysed cells. This does not turn independently fitted
#'   blocks into a joint shared-parameter model.
#' @param pooling_draws Number of independent posterior-draw combinations used
#'   to form empirical credible intervals for a pooled block summary.
#' @param pooling_seed Integer seed for reproducible posterior-draw pooling.
#' @param pooling_min_expected_cells Non-negative minimum posterior expected
#'   component occupancy required for a block to contribute to a
#'   component-specific pooled parameter. A zero-occupancy component therefore
#'   contributes zero weight without invalidating the fitted block.
#' @param ... Arguments passed to \code{diagnose_hicpotts_fit()} when \code{fit}
#'   is a raw fit or fit list.
#'
#' @return The parameter data frame with diagnostic metadata attached. For a
#'   block fit with `pool_blocks = "posterior_weighted"`, the returned
#'   `hicpotts_pooled_parameter_summary` contains overall rows and retains the
#'   per-block table in the `"block_summaries"` attribute.
#' @examples
#' diagnostic <- list(
#'     parameters = data.frame(parameter = "gamma", estimate = 0.30),
#'     reliability_flags = data.frame(
#'         criterion = "example", passed = TRUE, threshold = "illustrative"
#'     ),
#'     classification_entropy = NA_real_,
#'     warnings = character()
#' )
#' summarise_hicpotts_parameters(diagnostic)
#' @seealso \code{\link{diagnose_hicpotts_fit}},
#'   \code{\link{classify_hicpotts}}
#' @export
summarise_hicpotts_parameters <- function(
    fit, require_reliable = TRUE, x_vars = NULL,
    diagnose_covariates = TRUE, covariate_diagnostic_args = list(),
    pool_blocks = c("none", "posterior_weighted"),
    pooling_draws = 100000L, pooling_seed = 1L,
    pooling_min_expected_cells = 0, ...
) {
    pool_blocks <- match.arg(pool_blocks)
    if (!is.logical(require_reliable) || length(require_reliable) != 1L ||
        is.na(require_reliable)) {
        stop("'require_reliable' must be TRUE or FALSE.")
    }
    if (!is.logical(diagnose_covariates) || length(diagnose_covariates) != 1L ||
        is.na(diagnose_covariates)) {
        stop("'diagnose_covariates' must be TRUE or FALSE.")
    }
    if (!is.list(covariate_diagnostic_args) ||
        (length(covariate_diagnostic_args) &&
            is.null(names(covariate_diagnostic_args)))) {
        stop("'covariate_diagnostic_args' must be a named list.")
    }
    pooling_draws <- as.integer(pooling_draws)
    pooling_seed <- as.integer(pooling_seed)
    if (length(pooling_draws) != 1L || is.na(pooling_draws) ||
        pooling_draws < 1000L) {
        stop("'pooling_draws' must be an integer of at least 1000.")
    }
    if (length(pooling_seed) != 1L || is.na(pooling_seed)) {
        stop("'pooling_seed' must be one non-missing integer.")
    }
    if (length(pooling_min_expected_cells) != 1L ||
        !is.finite(pooling_min_expected_cells) ||
        pooling_min_expected_cells < 0) {
        stop("'pooling_min_expected_cells' must be non-negative.")
    }
    if (inherits(fit, "hicpotts_block_fit")) {
        return(.hicpotts_summarise_parameter_blocks(
            fit, require_reliable = require_reliable, x_vars = x_vars,
            diagnose_covariates = diagnose_covariates,
            covariate_diagnostic_args = covariate_diagnostic_args,
            pool_blocks = pool_blocks, pooling_draws = pooling_draws,
            pooling_seed = pooling_seed,
            pooling_min_expected_cells = pooling_min_expected_cells, ...
        ))
    }

    covariate_diagnostics <- NULL
    if (isTRUE(diagnose_covariates) && !is.null(x_vars)) {
        covariate_diagnostics <- do.call(
            diagnose_hicpotts_covariates,
            c(list(x_vars = x_vars), covariate_diagnostic_args)
        )
    } else if (isTRUE(diagnose_covariates) &&
        inherits(fit, "hicpotts_robust_fit")) {
        covariate_diagnostics <- fit$covariate_diagnostics
    }

    is_diagnostics <- is.list(fit) && is.data.frame(fit$parameters) &&
        is.data.frame(fit$reliability_flags)
    diagnostics <- if (inherits(fit, "hicpotts_robust_fit")) {
        fit$diagnostics
    } else if (is_diagnostics) {
        fit
    } else {
        diagnose_hicpotts_fit(
            fit,
            covariate_diagnostics = covariate_diagnostics, ...
        )
    }
    if (!is.list(diagnostics) || !is.data.frame(diagnostics$parameters) ||
        !is.data.frame(diagnostics$reliability_flags)) {
        stop("The supplied object does not contain valid HiCPotts diagnostics.")
    }

    ## Robust and precomputed diagnostic objects may already contain a
    ## covariate gate. When the user supplies x_vars here, make this call's
    ## diagnostic authoritative without re-running all MCMC diagnostics.
    if (!is.null(covariate_diagnostics) &&
        (inherits(fit, "hicpotts_robust_fit") || is_diagnostics)) {
        row <- match(
            "covariate_conditioning",
            diagnostics$reliability_flags$criterion
        )
        passed <- isTRUE(covariate_diagnostics$reliable)
        if (is.na(row)) {
            diagnostics$reliability_flags <- rbind(
                diagnostics$reliability_flags,
                data.frame(
                    criterion = "covariate_conditioning", passed = passed,
                    threshold = "correlation/condition thresholds",
                    stringsAsFactors = FALSE
                )
            )
        } else {
            diagnostics$reliability_flags$passed[row] <- passed
        }
        if (!passed) {
            diagnostics$warnings <- unique(c(
                diagnostics$warnings, covariate_diagnostics$warnings
            ))
        }
    }

    failed <- diagnostics$reliability_flags$criterion[
        !diagnostics$reliability_flags$passed
    ]

    ## Per-parameter resolution: a row is reportable when its own bulk and
    ## tail ESS reach the requested precision and, where more than one
    ## chain is available, its split R-hat is within tolerance. A
    ## non-finite value is a stuck or non-diagnostic chain and fails.
    out <- diagnostics$parameters
    ess_floor <- .hicpotts_flag_threshold(
        diagnostics$reliability_flags, "effective_sample_size", 200
    )
    rhat_ceiling <- .hicpotts_flag_threshold(
        diagnostics$reliability_flags, "split_rhat", 1.01
    )
    has_precision <- all(c("bulk_ESS", "tail_ESS", "Rhat") %in% names(out))
    if (has_precision) {
        ess_ok <- is.finite(out$bulk_ESS) & out$bulk_ESS >= ess_floor &
            is.finite(out$tail_ESS) & out$tail_ESS >= ess_floor
        rhat_ok <- if (any(is.finite(out$Rhat))) {
            is.finite(out$Rhat) & out$Rhat <= rhat_ceiling
        } else {
            rep(TRUE, nrow(out))
        }
        out$resolved <- ess_ok & rhat_ok
        unresolved <- out$parameter[!out$resolved]
    } else {
        ## No precision columns to judge, so make no claim either way.
        out$resolved <- NA
        unresolved <- character(0)
    }

    ## Check the common-target, movement, alignment, scale and conditioning
    ## criteria separately from the row-level Monte Carlo precision criteria.
    preconditions <- c(
        common_posterior_target = paste0(
            "chains do not share a frozen regression prior, so their ",
            "R-hat and ESS are not comparable; refit with ",
            "run_chain_betas(robust = TRUE)"),
        beta_movement = paste0(
            "a component recorded no accepted regression proposals; ",
            "review its trace and proposal settings"),
        relabelled = paste0(
            "components were not aligned across chains, so pooling ",
            "mixes different states"),
        coefficient_scale = paste0(
            "chains report coefficients on different scales, so the ",
            "pooled estimates are not commensurable"),
        covariate_conditioning = paste0(
            "the covariate diagnostics exceed the configured ",
            "conditioning thresholds")
    )
    unmet <- names(preconditions)[!vapply(
        names(preconditions),
        function(g) .hicpotts_flag_passed(diagnostics$reliability_flags, g),
        logical(1)
    )]
    if (isTRUE(require_reliable) && length(unmet)) {
        stop(
            "Parameter summary rejected (", unmet[1L], "): ",
            preconditions[[unmet[1L]]], ".",
            call. = FALSE
        )
    }
    if (isTRUE(require_reliable) && length(unresolved)) {
        template <- paste0(
            "Parameter summary criteria: %d of %d parameters are below ",
            "the configured ESS/R-hat thresholds (%s). Review the ",
            "diagnostics and chain length, or pass ",
            "require_reliable = FALSE for a table with a ",
            "'resolved' column."
        )
        stop(
            sprintf(
                template, length(unresolved), nrow(out),
                paste(utils::head(unresolved, 5), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    attr(out, "resolved") <- all(out$resolved)
    attr(out, "unresolved_parameters") <- unresolved
    attr(out, "failed_reliability_gates") <- failed
    attr(out, "warnings") <- diagnostics$warnings
    attr(out, "covariate_diagnostics") <- covariate_diagnostics
    class(out) <- c("hicpotts_parameter_summary", "data.frame")
    out
}

#' Run independently seeded HiCPotts chains
#'
#' @param N Number of rows and columns in the square Hi-C matrix.
#' @param gamma_prior Initial value for the Potts spatial parameter.
#' @param iterations Positive user-selected maximum MCMC updates per chain. A
#' chain may finish earlier when the MCSE precision rule is satisfied.
#' @param x_vars Named covariate list containing distance, GC, TES and ACC
#' matrices.
#' @param y Numeric N-by-N interaction-count matrix.
#' @param theta_start Optional initial zero-inflation probability for ZIP/ZINB.
#' @param size_start Length-three positive initial NB2 size vector for NB/ZINB.
#' @param use_data_priors Whether to estimate component-specific empirical-
#' Bayes Normal hyperparameters from soft allocation probabilities during
#' warm-up and freeze them before retained posterior draws.
#' @param user_fixed_priors Optional component prior list, required when
#' `use_data_priors = FALSE`.
#' @param dist One of `Poisson`, `NB`, `ZIP` or `ZINB`.
#' @param epsilon Optional positive ABC kernel bandwidth. When NULL, one
#' deterministic prior-predictive calibration is shared by every chain;
#' supplying a value skips calibration simulations entirely.
#' @param distance_metric Retained for source compatibility and ignored.
#' @param seeds Integer seeds, one per chain.
#' @param initialization One or more latent-state initialization strategies:
#' `likelihood_informed`, `count_quantile`, `distance_adjusted`,
#' `noise_anchored_random` or unrestricted `random`. Values are recycled across
#' seeds. Likelihood-informed initialization fits a non-spatial hard mixture
#' under the selected count family; it changes only starting values.
#' @param mcse_stop Whether to stop a chain when all monitored parameters meet
#' the batch-means Monte Carlo standard-error criterion.
#' @param mcse_min_iterations Minimum updates before MCSE stopping is checked.
#' @param mcse_check_interval Updates between MCSE checks.
#' @param mcse_relative_threshold Required maximum MCSE divided by posterior
#' standard deviation across monitored beta, gamma, theta and size chains.
#' @param tempering_warmup Number of optional cyclical latent-state heating
#' updates. Zero disables heating and is the tested general-use default.
#' @param tempering_beta_min Lowest inverse temperature during optional heating.
#' @param tempering_cycle Positive length of each optional heating cycle.
#' @param relabel Whether to relabel components 2 and 3 before returning fits.
#' @param z_probability_burnin Integer iteration after which latent-state
#' membership frequencies are accumulated for classification. The default
#' \code{-1} uses the sampler's automatic rule.
#' @param abc_epsilon_quantile Quantile used to calibrate the Gaussian ABC
#' tolerance; defaults to 0.10.
#' @param gamma_update_interval Positive integer interval between ABC gamma
#' updates. The default 5 retains four cheap allocation/parameter sweeps between
#' expensive auxiliary-Potts updates.
#' @param abc_potts_sweeps Non-negative auxiliary-Potts sweep count; zero uses
#' the lattice-size-dependent default.
#' @param abc_sim_reps Positive number of auxiliary fields averaged per ABC
#' proposal.
#' @param use_noise_relationship_prior,noise_link_sd Switch and link
#' scale for the component-1/3 biological relationship prior.
#' @param noise_order_strength,noise_order_width Ordering strength and
#' width for the component-1/3 biological relationship prior.
#' @param comp23_barrier_kappa,comp23_barrier_w Component-2/3 intercept-barrier
#' settings; enabled by default.
#' @param branch_swap_interval Non-negative interval for the reversible whole
#' component-2/component-3 branch transition; zero disables it.
#' @param signal_block_move_interval Non-negative interval for the connected
#' component-2/component-3 cluster move; zero disables it.
#' @param gamma_large_jump_probability,gamma_large_jump_multiplier Settings for
#' the wide logit-normal component of the gamma proposal mixture.
#' @param gamma_independence_probability Probability of an independent gamma
#' proposal from its configured Beta prior.
#' @param gamma_method Retained for compatibility; only \code{"abc"} is
#' supported.
#' @param mc_cores Positive number of independent chains to execute
#' concurrently. Windows uses a PSOCK cluster; Unix-like systems use forked
#' workers. Explicit per-chain seeds make both paths reproducible.
#' @param verbose Whether each sampler should print periodic progress.
#' @param progress_interval Positive number of iterations between progress
#' messages when `verbose = TRUE`.
#' @return A list of fits suitable for `diagnose_hicpotts_fit()`.
#' @examples
#' N <- 5
#' y <- matrix(rpois(N * N, 4), N, N)
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' fits <- run_hicpotts_chains(
#'     N = N, iterations = 5, x_vars = x_vars, y = y,
#'     dist = "ZINB", theta_start = 0.5, size_start = c(1, 1, 1), seeds = 1:2
#' )
#' length(fits)
#' @noRd
run_hicpotts_chains <- function(
    N, gamma_prior = 0.3, iterations, x_vars, y,
    theta_start = NULL, size_start = NULL, use_data_priors = TRUE,
    user_fixed_priors = NULL, dist = "ZIP", epsilon = NULL,
    distance_metric = "manhattan", seeds = seq_len(4),
    initialization = c(
        "likelihood_informed", "likelihood_informed",
        "distance_adjusted", "noise_anchored_random"
    ),
    mcse_stop = TRUE, mcse_min_iterations = 10000L,
    mcse_check_interval = 500L, mcse_relative_threshold = 0.05,
    tempering_warmup = 0L, tempering_beta_min = 0.30,
    tempering_cycle = 500L,
    z_probability_burnin = -1L,
    gamma_prior_shape1 = 1, gamma_prior_shape2 = 1,
    abc_epsilon_quantile = 0.10, gamma_update_interval = 5L,
    abc_potts_sweeps = 0L, abc_sim_reps = 4L,
    comp23_barrier_kappa = 10, comp23_barrier_w = 0.3,
    use_noise_relationship_prior = TRUE,
    noise_link_sd = 0.5, noise_order_strength = 10,
    noise_order_width = 0.5,
    branch_swap_interval = 1L, signal_block_move_interval = 5L,
    gamma_large_jump_probability = 0.10,
    gamma_large_jump_multiplier = 4,
    gamma_independence_probability = 0.05,
    gamma_method = "abc",
    relabel = TRUE,
    mc_cores = 1L,
    verbose = FALSE,
    progress_interval = 50L
) {
    iterations <- as.integer(iterations)
    if (length(iterations) != 1L || is.na(iterations) || iterations < 1L) {
        stop("'iterations' must be one positive integer.")
    }
    seeds <- as.integer(seeds)
    if (!length(seeds) || anyNA(seeds)) stop("'seeds' must contain integers.")
    initialization <- match.arg(initialization,
        c(
            "likelihood_informed", "count_quantile", "distance_adjusted",
            "noise_anchored_random", "random"
        ),
        several.ok = TRUE
    )
    initialization <- rep(initialization, length.out = length(seeds))
    gamma_method <- match.arg(gamma_method)
    if (length(relabel) != 1L || is.na(relabel)) {
        stop(
            "'relabel' must be TRUE or FALSE."
        )
    }
    .hicpotts_validate_fit_inputs(
        N = N, y = y, x_vars = x_vars, iterations = iterations,
        gamma_prior = gamma_prior, theta_start = theta_start,
        size_start = size_start, dist = dist,
        use_data_priors = use_data_priors,
        user_fixed_priors = user_fixed_priors
    )
    abc_calibration <- .hicpotts_prepare_shared_abc_epsilon(
        epsilon = epsilon, N = N,
        gamma_prior_shape1 = gamma_prior_shape1,
        gamma_prior_shape2 = gamma_prior_shape2,
        abc_epsilon_quantile = abc_epsilon_quantile,
        abc_potts_sweeps = abc_potts_sweeps,
        abc_sim_reps = abc_sim_reps
    )
    initial_z <- Map(function(seed, init_method) {
        make_hicpotts_initial_z(
            y, x_vars = x_vars, method = init_method, seed = seed,
            dist = dist, size_start = size_start
        )
    }, seeds, initialization)
    tasks <- Map(function(seed, init_method, z_start) {
        arguments <- list(
            N = N, gamma_prior = gamma_prior,
            iterations = iterations, x_vars = x_vars, y = y,
            use_data_priors = use_data_priors,
            user_fixed_priors = user_fixed_priors,
            dist = dist, epsilon = abc_calibration$epsilon,
            distance_metric = distance_metric,
            size_start = size_start, theta_start = theta_start,
            mcse_stop = mcse_stop, mcse_min_iterations = mcse_min_iterations,
            mcse_check_interval = mcse_check_interval,
            mcse_relative_threshold = mcse_relative_threshold,
            tempering_warmup = tempering_warmup,
            tempering_beta_min = tempering_beta_min,
            tempering_cycle = tempering_cycle,
            gamma_prior_shape1 = gamma_prior_shape1,
            gamma_prior_shape2 = gamma_prior_shape2,
            abc_epsilon_quantile = abc_epsilon_quantile,
            gamma_update_interval = gamma_update_interval,
            abc_potts_sweeps_arg = abc_potts_sweeps,
            abc_sim_reps = abc_sim_reps,
            z_probability_burnin_arg = z_probability_burnin,
            comp23_barrier_kappa = comp23_barrier_kappa,
            comp23_barrier_w = comp23_barrier_w,
            use_noise_relationship_prior = use_noise_relationship_prior,
            noise_link_sd = noise_link_sd,
            noise_order_strength = noise_order_strength,
            noise_order_width = noise_order_width,
            branch_swap_interval = branch_swap_interval,
            signal_block_move_interval = signal_block_move_interval,
            gamma_large_jump_probability = gamma_large_jump_probability,
            gamma_large_jump_multiplier = gamma_large_jump_multiplier,
            gamma_independence_probability = gamma_independence_probability,
            gamma_method = gamma_method,
            z_start = z_start,
            verbose = verbose,
            progress_interval = progress_interval
        )
        list(
            arguments = arguments, y = y, x_vars = x_vars,
            seed = seed, initialization = init_method
        )
    }, seeds, initialization, initial_z)
    fits <- .hicpotts_parallel_tasks(tasks, mc_cores = mc_cores)
    fits <- .hicpotts_record_shared_abc_calibration(fits, abc_calibration)
    names(fits) <- paste0("chain", seq_along(fits))
    if (isTRUE(relabel)) fits <- relabel_hicpotts(fits)
    fits
}

#' Compare NB and ZINB fits as a sensitivity analysis
#'
#' Fits identical data under NB and ZINB. Use `diagnose_hicpotts_fit()` to
#' identify conclusions sensitive to the zero-inflation assumption.
#'
#' @param seeds Integer seeds used for each family.
#' @inheritParams run_hicpotts_chains
#' @return A named list containing NB and ZINB repeated-chain fits.
#' @examples
#' N <- 6
#' y <- matrix(rpois(N * N, 4), N, N)
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' fam <- run_hicpotts_family_sensitivity(
#'     N = N, iterations = 5, x_vars = x_vars,
#'     y = y, size_start = c(1, 1, 1), seeds = 1:2
#' )
#' names(fam)
#' @noRd
run_hicpotts_family_sensitivity <- function(
    N, gamma_prior = 0.3, iterations,
    x_vars, y, theta_start = 0.5, size_start, use_data_priors = TRUE,
    user_fixed_priors = NULL, epsilon = NULL,
    distance_metric = "manhattan", seeds = seq_len(4),
    initialization = c(
        "likelihood_informed", "likelihood_informed",
        "distance_adjusted", "noise_anchored_random"
    ),
    relabel = TRUE
) {
    common <- list(
        N = N, gamma_prior = gamma_prior, iterations = iterations,
        x_vars = x_vars, y = y, size_start = size_start,
        use_data_priors = use_data_priors,
        user_fixed_priors = user_fixed_priors,
        epsilon = epsilon, distance_metric = distance_metric, seeds = seeds,
        initialization = initialization, relabel = relabel
    )
    list(
        NB = do.call(run_hicpotts_chains, c(common, list(
            dist = "NB",
            theta_start = NULL
        ))),
        ZINB = do.call(run_hicpotts_chains, c(common, list(
            dist = "ZINB",
            theta_start = theta_start
        )))
    )
}

## Read a numeric threshold back out of the reliability_flags table so the
## summary honours whatever limits diagnose_hicpotts_fit() was given,
## rather than re-hardcoding them.
.hicpotts_flag_threshold <- function(flags, criterion, default) {
    row <- match(criterion, flags$criterion)
    if (is.na(row)) {
        return(default)
    }
    ## Check the text is a bare number before coercing, rather than coercing
    ## and muffling the NA warning.
    cleaned <- gsub("[^0-9.]", "", flags$threshold[row])
    if (length(cleaned) != 1L || !grepl("^[0-9]+([.][0-9]+)?$", cleaned)) {
        return(default)
    }
    value <- as.numeric(cleaned)
    if (!is.finite(value)) default else value
}

## TRUE when a named gate passed; TRUE when the gate is absent, since an
## absent gate is not evidence of failure.
.hicpotts_flag_passed <- function(flags, criterion) {
    row <- match(criterion, flags$criterion)
    if (is.na(row)) {
        return(TRUE)
    }
    isTRUE(flags$passed[row])
}

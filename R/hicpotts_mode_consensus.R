.hicpotts_parameter_z <- function(fit) {
    z <- attr(fit$z_final, "z_parameter_final", exact = TRUE)
    if (is.null(z) && is.matrix(
        fit$z_parameter_final
    )) {
        z <- fit$z_parameter_final
    }
    if (is.null(z)) z <- fit$z_final
    z
}

## Post-burn-in membership probabilities as a cells x 3 matrix, or NULL when a
## fit does not carry them.
.hicpotts_soft_z <- function(fit) {
    zp <- fit$z_probabilities
    ok <- is.list(zp) && length(zp) == 3L &&
        all(vapply(zp, function(m) is.matrix(m) && is.numeric(m), logical(1)))
    if (!ok) {
        return(NULL)
    }
    d <- dim(zp[[1L]])
    if (!all(vapply(zp, function(m) identical(dim(m), d), logical(1)))) {
        return(NULL)
    }
    cbind(as.vector(zp[[1L]]), as.vector(zp[[2L]]), as.vector(zp[[3L]]))
}

## Mean per-cell total-variation agreement between two biologically labelled
## soft allocations.
##
## Component 3 is coupled to component 1 by the posterior and component 2 is
## unrestricted, so a 2/3 exchange is a scientific disagreement, not a harmless
## label permutation. Agreement is therefore evaluated directly.
##
## Agreement is 1 - mean total-variation distance, where per cell
##   TV = 0.5 * sum_k |p_ik - p_jk|  in [0, 1].
## This is a soft generalisation of the old "fraction of cells with identical
## hard labels": it degrades smoothly as two chains' probabilities drift apart
## instead of jumping only when an argmax flips, so a cell at 0.51/0.49 no
## longer counts as a full disagreement with a cell at 0.49/0.51.
.hicpotts_soft_agreement <- function(p, q) {
    list(
        agreement = 1 - mean(0.5 * rowSums(abs(p - q))),
        swap23 = FALSE
    )
}

.hicpotts_mode_consensus <- function(fits, agreement_threshold = 0.8,
                                    minimum_mode_chains = 2L) {
    n <- length(fits)
    if (n < 1L) stop("At least one fit is required for mode screening.")
    if (!is.numeric(agreement_threshold) || length(agreement_threshold) != 1L ||
        !is.finite(agreement_threshold) || agreement_threshold <= 0 ||
        agreement_threshold > 1) {
        stop("'agreement_threshold' must be in (0, 1].")
    }
    minimum_mode_chains <- as.integer(minimum_mode_chains)
    if (length(minimum_mode_chains) != 1L || is.na(minimum_mode_chains) ||
        minimum_mode_chains < 1L) {
        stop("'minimum_mode_chains' must be one positive integer.")
    }

    ## Prefer the soft post-burn-in membership probabilities. They use every
    ## retained draw rather than one final allocation, so the mode structure is
    ## estimated from the posterior rather than from a single sweep that may
    ## itself be an excursion. Fits predating z_probabilities fall back to the
    ## hard allocations, and the basis is recorded either way.
    soft <- lapply(fits, .hicpotts_soft_z)
    use_soft <- !any(vapply(soft, is.null, logical(1)))

    swap23 <- matrix(FALSE, n, n)
    if (use_soft) {
        cells <- vapply(soft, nrow, integer(1))
        if (length(unique(cells)) != 1L) {
            stop(
                "All fits must contain equally sized latent-state probability ",
                "fields."
            )
        }
        distance_basis <- "soft z_probabilities (biologically labelled)"
        agreement <- diag(1, n)
        if (n > 1L) {
            for (i in seq_len(n - 1L)) {
                for (j in (i + 1L):n) {
                    cmp <- .hicpotts_soft_agreement(soft[[i]], soft[[j]])
                    agreement[i, j] <- agreement[j, i] <- cmp$agreement
                    swap23[i, j] <- swap23[j, i] <- cmp$swap23
                }
            }
        }
    } else {
        z <- lapply(fits, .hicpotts_parameter_z)
        dims <- lapply(z, dim)
        if (any(vapply(z, function(x) !is.matrix(x), logical(1))) ||
            !all(vapply(dims, identical, logical(1), dims[[1L]]))) {
            stop(
                "All fits must contain equally sized internal latent-state ",
                "matrices."
            )
        }
        distance_basis <- "hard final allocations (no stored z_probabilities)"
        agreement <- diag(1, n)
        if (n > 1L) {
            for (i in seq_len(n - 1L)) {
                for (j in (i + 1L):n) {
                    direct <- mean(z[[i]] == z[[j]])
                    agreement[i, j] <- agreement[j, i] <- direct
                }
            }
        }
    }
    component <- if (n == 1L) {
        1L
    } else {
        stats::cutree(
            stats::hclust(stats::as.dist(1 - agreement), method = "complete"),
            h = 1 - agreement_threshold
        )
    }
    component <- match(component, unique(component))
    sizes <- tabulate(component)
    candidates <- which(sizes == max(sizes))
    if (length(candidates) > 1L) {
        within_score <- vapply(candidates, function(g) {
            idx <- which(component == g)
            mean(agreement[idx, idx, drop = FALSE])
        }, numeric(1))
        chosen <- candidates[which.max(within_score)]
    } else {
        chosen <- candidates
    }
    selected <- which(component == chosen)
    consensus_available <- length(selected) >= minimum_mode_chains
    if (!consensus_available) selected <- seq_len(n)

    ## M4: a materially populated SECOND mode is multimodality, not a set of
    ## outliers to discard. Previously the largest group was selected and the
    ## rest merely warned about, so a 2-vs-2 or 4-vs-2 split could report
    ## consensus_available = TRUE and let within-selected-mode diagnostics pass
    ## while a replicated alternative mode was silently omitted from pooled
    ## uncertainty. Any other mode holding at least `minimum_mode_chains` chains
    ## now marks the fit multimodal, and the reliability gate treats that as a
    ## failure of WHOLE-POSTERIOR reliability. Mode-conditional summaries remain
    ## available via `selected`; they are simply no longer presented as
    ## describing the full posterior.
    other_sizes <- sizes[-chosen]
    competing_modes <- sum(other_sizes >= minimum_mode_chains)
    multimodal <- competing_modes > 0L

    list(
        selected = selected,
        excluded = setdiff(seq_len(n), selected),
        agreement = agreement,
        component = component,
        component_sizes = sizes,
        agreement_threshold = agreement_threshold,
        minimum_mode_chains = minimum_mode_chains,
        consensus_available = consensus_available,
        multimodal = multimodal,
        competing_modes = competing_modes,
        distance_basis = distance_basis,
        ## Retained for structure compatibility. It is always FALSE because a
        ## 2/3 exchange changes the biological interpretation in current fits.
        swap23 = swap23
    )
}

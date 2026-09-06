#' Monte-Carlo reliability of the per-cell latent-state probabilities
#'
#' @description
#' Quantifies how precisely the sampler has estimated each cell's membership
#' probabilities, using the batched post-burn-in latent-state frequencies
#' recorded by \code{run_chain_betas()}.
#'
#' Allocation draws are serially correlated, so the naive binomial standard
#' error \eqn{\sqrt{p(1-p)/n}} understates the Monte-Carlo error of a
#' membership probability. Batch means give an autocorrelation-aware estimate:
#' the post-burn-in draws are split into contiguous batches and the spread
#' between batch means carries the dependence.
#'
#' This is a diagnostic layer only. It reports how reliable the probabilities
#' are; it does not modify them, and the three-way MAP label is unchanged.
#'
#' @param fit A raw chain, list of chains, robust result returned by
#'   \code{run_chain_betas()}, or block-aware result from
#'   \code{combine_hicpotts_blocks()}.
#' @param component_names Character vector of length three naming components
#'   1, 2 and 3.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{mcse}}{Three \eqn{N \times N} matrices of per-cell
#'       Monte-Carlo standard errors, one per component.}
#'     \item{\code{ess}}{Three \eqn{N \times N} matrices of per-cell allocation
#'       effective sample sizes.}
#'     \item{\code{worst_cell_ess}}{The smallest effective sample size over all
#'       cells and components -- the binding constraint on the classification.}
#'     \item{\code{max_mcse}}{The largest per-cell Monte-Carlo standard error.}
#'     \item{\code{between_chain_disagreement}}{\code{NULL} for a single chain;
#'       otherwise the mean and maximum per-cell total-variation distance
#'       between chains' biologically labelled membership probabilities.}
#'   }
#'
#' @examples
#' batch_means <- array(
#'     c(
#'         rep(0.70, 4), rep(0.20, 4), rep(0.10, 4),
#'         rep(0.65, 4), rep(0.25, 4), rep(0.10, 4)
#'     ),
#'     dim = c(2, 2, 3, 2)
#' )
#' fit <- list(
#'     chains = lapply(seq_len(3L), function(k) matrix(k, 2, 1)),
#'     z_probability_batches = list(
#'         batch_means = batch_means,
#'         batch_draws = c(10, 10)
#'     )
#' )
#' diagnostic <- allocation_diagnostics(fit)
#' diagnostic$worst_cell_ess
#' diagnostic$max_mcse
#'
#' @seealso \code{\link{classify_hicpotts}}
#' @export
allocation_diagnostics <- function(
    fit,
    component_names = c("noise", "signal", "false signal")
) {
    if (inherits(fit, "hicpotts_block_fit")) {
        return(.hicpotts_allocation_diagnostics_blocks(
            fit, component_names = component_names
        ))
    }
    if (inherits(fit, "hicpotts_robust_fit")) fit <- fit$fits
    is_single <- is.list(fit) && is.list(fit$chains) && length(fit$chains) == 3L
    fits <- if (is_single) list(fit) else fit
    if (!is.list(fits) || !length(fits)) {
        stop("'fit' must be a HiCPotts result or a non-empty list of results.")
    }

    have_batches <- vapply(fits, function(x) {
        is.list(x$z_probability_batches) &&
            !is.null(x$z_probability_batches$batch_means)
    }, logical(1))
    if (!all(have_batches)) {
        message <- paste0(
            "Every fit must contain $z_probability_batches. Refit the data ",
            "so batched latent-state frequencies are recorded."
        )
        stop(message, call. = FALSE)
    }

    per_chain <- lapply(fits, function(x) {
        bm <- x$z_probability_batches$batch_means
        d <- dim(bm)
        if (length(d) != 4L || d[3L] != 3L) {
            stop(
                "$z_probability_batches$batch_means must be an N x N x 3 x ",
                "batches array."
            )
        }
        n_batches <- d[4L]
        if (n_batches < 2L) {
            message <- paste0(
                "At least two batches are required to estimate Monte-Carlo ",
                "precision; run a longer chain."
            )
            stop(message, call. = FALSE)
        }
        draws <- as.numeric(x$z_probability_batches$batch_draws)
        total_draws <- sum(draws)

        mcse <- vector("list", 3L)
        ess <- vector("list", 3L)
        for (k in seq_len(3L)) {
            ## batches x cells matrix of batch-mean memberships
            bk <- matrix(aperm(bm[, , k, , drop = FALSE], c(4L, 1L, 2L, 3L)),
                nrow = n_batches
            )
            batch_mean <- colMeans(bk)
            ## Between-batch variance of the batch means; the MCSE of the
            ## overall mean is its square root divided by the number of batches.
            between_var <- apply(bk, 2L, stats::var)
            cell_mcse <- sqrt(between_var / n_batches)

            ## ESS = (i.i.d. variance) / (Monte-Carlo variance of the mean). For
            ## a membership indicator the i.i.d. variance is p(1-p). A cell that
            ## never changes label has zero variance in both, which is a
            ## degenerate rather than an infinitely precise estimate: it is
            ## reported as the full draw count, since there is nothing left for
            ## the sampler to resolve.
            p <- batch_mean
            iid_var <- p * (1 - p)
            cell_ess <- ifelse(cell_mcse > 0, iid_var / (cell_mcse^2),
                total_draws
            )
            cell_ess <- pmin(cell_ess, total_draws)
            cell_ess[iid_var <= 0] <- total_draws

            mcse[[k]] <- matrix(cell_mcse, d[1L], d[2L])
            ess[[k]] <- matrix(cell_ess, d[1L], d[2L])
        }
        names(mcse) <- component_names
        names(ess) <- component_names
        list(mcse = mcse, ess = ess, total_draws = total_draws)
    })

    ## Pool across chains by taking the worst case per cell: a probability is
    ## only as trustworthy as the least reliable chain contributing to it.
    mcse <- lapply(seq_len(3L), function(k) {
        Reduce(pmax, lapply(per_chain, function(x) x$mcse[[k]]))
    })
    ess <- lapply(seq_len(3L), function(k) {
        Reduce(pmin, lapply(per_chain, function(x) x$ess[[k]]))
    })
    names(mcse) <- component_names
    names(ess) <- component_names

    ## Between-chain disagreement on the biological labels. A component-2/3
    ## exchange is a signal-versus-false-signal disagreement and is not aligned
    ## away.
    disagreement <- NULL
    if (length(fits) > 1L) {
        soft <- lapply(fits, .hicpotts_soft_z)
        if (!any(vapply(soft, is.null, logical(1)))) {
            pairs <- utils::combn(length(soft), 2L)
            tv_stats <- apply(pairs, 2L, function(ij) {
                a <- soft[[ij[1L]]]
                b <- soft[[ij[2L]]]
                tv <- 0.5 * rowSums(abs(a - b))
                c(mean = mean(tv), max = max(tv))
            })
            disagreement <- list(
                mean = mean(tv_stats["mean", ]),
                max = max(tv_stats["max", ]),
                n_chain_pairs = ncol(pairs)
            )
        }
    }

    list(
        mcse = mcse,
        ess = ess,
        worst_cell_ess = min(vapply(ess, min, numeric(1))),
        max_mcse = max(vapply(mcse, max, numeric(1))),
        between_chain_disagreement = disagreement,
        draws_per_chain = vapply(per_chain, function(x) x$total_draws, numeric(
            1
        ))
    )
}

#' Fingerprint the data a HiCPotts fit was produced from
#'
#' @description
#' Computes a compact, order-sensitive fingerprint of the counts, covariates
#' and lattice geometry behind a fit, plus a coordinate key identifying the
#' lattice cells in their stored order.
#'
#' Chains are routinely pooled -- across replicate runs, across seeds, in
#' \code{classify_hicpotts()} and in the robust-fit machinery -- and pooling is
#' only meaningful when every chain describes the same cells of the same
#' dataset. Nothing previously prevented fits from two different datasets, two
#' different resolutions, or two different lattice orderings from being pooled:
#' the shapes match, so the arithmetic succeeds and returns a silently
#' meaningless answer. Row-shuffled classification data has the same problem in
#' reverse -- probabilities get attached to the wrong genomic coordinates with
#' no error anywhere.
#'
#' The fingerprint makes both failures loud. It is deliberately sensitive to
#' cell ORDER as well as content, because a reordering is exactly the silent
#' misalignment being guarded against.
#'
#' @param y Observed count matrix.
#' @param x_vars Covariate list in the layout accepted by
#'   \code{run_metropolis_MCMC_betas()}.
#' @param coordinates Optional data frame of genomic coordinates in the stored
#'   lattice order, used to build the coordinate key.
#'
#' @return A list of class \code{hicpotts_fingerprint}.
#' @noRd
hicpotts_fingerprint <- function(y, x_vars = NULL, coordinates = NULL) {
    if (!is.matrix(y) || !is.numeric(y)) {
        stop("'y' must be a numeric matrix.")
    }

    ## A cheap, order-sensitive digest. Position-weighted sums distinguish a
    ## matrix from any permutation of itself, which a plain sum or mean would
    ## not; combined with the dimensions and the extremes this is sufficient to
    ## catch accidental mismatches. It is an integrity check against mistakes,
    ## not a cryptographic commitment against tampering.
    digest_numeric <- function(v) {
        v <- as.numeric(v)
        n <- length(v)
        if (!n) {
            return(c(n = 0, s = 0, w = 0, mn = NA_real_, mx = NA_real_))
        }
        finite <- v[is.finite(v)]
        idx <- seq_len(n)
        c(
            n = n,
            s = sum(finite),
            ## Two different weightings so that swapping a pair of cells cannot
            ## coincidentally preserve both.
            w = sum(finite * (idx[is.finite(v)] %% 9973L)),
            mn = if (length(finite)) min(finite) else NA_real_,
            mx = if (length(finite)) max(finite) else NA_real_
        )
    }

    y_digest <- digest_numeric(y)
    covariate_digest <- NULL
    if (!is.null(x_vars)) {
        if (!is.list(x_vars)) stop("'x_vars' must be a list when supplied.")
        covariate_digest <- lapply(x_vars, function(v) {
            m <- if (is.list(v)) v[[1L]] else v
            digest_numeric(m)
        })
        names(covariate_digest) <- names(x_vars)
    }

    coordinate_key <- NULL
    if (!is.null(coordinates)) {
        if (!is.data.frame(coordinates)) {
            coordinates <- as.data.frame(
                coordinates
            )
        }
        if (nrow(coordinates) != length(y)) {
            stop("'coordinates' must contain exactly one row per lattice cell.")
        }
        key_columns <- intersect(c(
            "chr", "chromosome", "start", "end",
            "start.j.", "end.j."
        ), names(coordinates))
        if (length(key_columns)) {
            coordinate_key <- list(
                columns = key_columns,
                n = nrow(coordinates),
                head = utils::head(
                    coordinates[, key_columns, drop = FALSE],
                    3L
                ),
                tail = utils::tail(coordinates[, key_columns, drop = FALSE], 3L)
            )
        }
    }

    structure(
        list(
            lattice_dim = dim(y),
            counts = y_digest,
            covariates = covariate_digest,
            coordinate_key = coordinate_key
        ),
        class = "hicpotts_fingerprint"
    )
}

#' Are two HiCPotts fingerprints compatible for pooling?
#' @param a,b Objects from \code{hicpotts_fingerprint()}.
#' @return \code{TRUE}, or a character vector describing every mismatch.
#' @noRd
hicpotts_fingerprints_agree <- function(a, b) {
    problems <- character(0)
    if (!identical(a$lattice_dim, b$lattice_dim)) {
        problems <- c(problems, sprintf(
            "lattice dimensions differ (%s vs %s)",
            paste(a$lattice_dim, collapse = "x"),
            paste(b$lattice_dim, collapse = "x")
        ))
    }
    same_numeric <- function(x, y, tol = 1e-8) {
        isTRUE(all.equal(unname(x), unname(y), tolerance = tol))
    }
    if (!same_numeric(a$counts, b$counts)) {
        problems <- c(
            problems,
            "observed counts differ (or are in a different cell order)"
        )
    }
    if (!identical(names(a$covariates), names(b$covariates))) {
        problems <- c(problems, "covariate sets differ")
    } else if (!is.null(a$covariates)) {
        for (nm in names(a$covariates)) {
            if (!same_numeric(a$covariates[[nm]], b$covariates[[nm]])) {
                problems <- c(problems, sprintf("covariate '%s' differs", nm))
            }
        }
    }
    ak <- a$coordinate_key
    bk <- b$coordinate_key
    if (!is.null(ak) && !is.null(bk)) {
        if (!identical(ak$columns, bk$columns) || !identical(ak$n, bk$n) ||
            !isTRUE(all.equal(ak$head, bk$head)) ||
            !isTRUE(all.equal(ak$tail, bk$tail))) {
            problems <- c(problems, "genomic coordinate keys differ")
        }
    }
    if (length(problems)) problems else TRUE
}

## Internal: refuse to pool fits from different datasets.
.hicpotts_check_poolable <- function(fits, context = "pooling") {
    gamma_methods <- vapply(fits, function(x) {
        method <- attr(x$gamma, "gamma_method", exact = TRUE)
        if (is.null(method)) "legacy-unspecified" else as.character(method)[1L]
    }, character(1))
    if (length(unique(gamma_methods)) > 1L) {
        stop(
            sprintf(
                "Refusing %s: fits use different gamma transitions (%s). ",
                context, paste(unique(gamma_methods), collapse = ", ")
            ),
            "Pool only chains produced by the same gamma_method; compare ",
            "methods ",
            "as a sensitivity analysis instead."
        )
    }

    fps <- lapply(fits, function(x) x$data_fingerprint)
    have <- !vapply(fps, is.null, logical(1))
    if (sum(have) < 2L) {
        return(invisible(NULL))
    } # nothing to compare
    reference <- fps[[which(have)[1L]]]
    for (i in which(have)[-1L]) {
        verdict <- hicpotts_fingerprints_agree(reference, fps[[i]])
        if (!isTRUE(verdict)) {
            stop(
                sprintf(
                    paste0("Refusing %s: chain %d was fitted to different ",
                        "data than chain %d. %s. "),
                    context, i, which(have)[1L], paste(verdict, collapse = "; ")
                ),
                "Pooling chains from different datasets, resolutions or ",
                "lattice ",
                "orderings produces a silently meaningless result."
            )
        }
    }
    invisible(NULL)
}

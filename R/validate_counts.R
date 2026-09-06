# no roxygen needed as it is a helper function
#'
#' Shared count-matrix validator for every public fitting boundary (M3).
#'
#' Integer validation previously existed only inside \code{process_data()}.
#' Users calling the exported fitters directly - \code{run_chain_betas()},
#' \code{run_metropolis_MCMC_betas()}, \code{hicpotts_robust_estimation()} -
#' could pass normalised/balanced Hi-C values (KR, ICE, VC), which are
#' continuous. Poisson/NB densities then return \code{-Inf} for every
#' non-integer cell, invalidating the likelihood while emitting one warning per
#' element. This validator is applied at each of those boundaries so the same
#' contract holds however the package is entered.
#'
#' @param y A numeric count matrix.
#' @param what Name used in error messages.
#' @param require_square Logical; require an N x N matrix.
#' @return Invisibly \code{TRUE}; errors otherwise.
#' @noRd
.hicpotts_validate_counts <- function(y, what = "y", require_square = TRUE) {
    if (!is.matrix(y) || !is.numeric(y)) {
        stop(sprintf("'%s' must be a numeric matrix of counts.", what),
            call. = FALSE
        )
    }
    if (require_square && nrow(y) != ncol(y)) {
        stop(sprintf(
            "'%s' must be square (N x N); got %d x %d.",
            what, nrow(y), ncol(y)
        ), call. = FALSE)
    }
    if (anyNA(y)) {
        stop(sprintf(
            "'%s' contains NA values; handle missing data before fitting.",
            what
        ), call. = FALSE)
    }
    if (any(!is.finite(y))) {
        stop(sprintf("'%s' contains non-finite values.", what), call. = FALSE)
    }
    if (any(y < 0)) {
        stop(sprintf(
            "'%s' contains negative values; counts must be non-negative.",
            what
        ), call. = FALSE)
    }

    n_noninteger <- sum(abs(y - round(y)) > 1e-8)
    if (n_noninteger > 0L) {
        template <- paste0(
            "'%s' has %d of %d non-integer values. HiCPotts is a count model ",
            "(Poisson/NB/ZIP/ZINB) and requires RAW integer contact counts, ",
            "not normalised/balanced Hi-C values (KR, ICE, VC). Supply the ",
            "raw counts, or round them yourself if discretising is intended."
        )
        stop(sprintf(template, what, n_noninteger, length(y)), call. = FALSE)
    }

    invisible(TRUE)
}

#'
#' Validate covariate matrices against the response (M3).
#' @noRd
.hicpotts_validate_covariates <- function(x_vars, N, index = 1L) {
    required <- c("distance", "GC", "TES", "ACC")
    if (!is.list(x_vars) || !all(required %in% names(x_vars))) {
        stop("'x_vars' must be a named list with entries: ",
            paste(required, collapse = ", "),
            call. = FALSE
        )
    }
    for (nm in required) {
        entry <- x_vars[[nm]]
        m <- if (is.list(entry)) entry[[index]] else entry
        if (!is.matrix(m) || !is.numeric(m)) {
            stop(sprintf("x_vars$%s must be a numeric N x N matrix.", nm),
                call. = FALSE
            )
        }
        if (nrow(m) != N || ncol(m) != N) {
            stop(sprintf(
                "x_vars$%s is %d x %d but the response is %d x %d.",
                nm, nrow(m), ncol(m), N, N
            ), call. = FALSE)
        }
        if (anyNA(m) || any(!is.finite(m))) {
            stop(
                sprintf(
                    "x_vars$%s contains missing or non-finite values.",
                    nm
                ),
                call. = FALSE
            )
        }
    }
    invisible(TRUE)
}

#'
#' Validate sampler starting values (M3).
#' @noRd
.hicpotts_validate_starts <- function(gamma_prior = NULL, theta_start = NULL,
                                        size_start = NULL, dist = NULL) {
    if (!is.null(dist) && !dist %in% c("Poisson", "NB", "ZIP", "ZINB")) {
        stop("'dist' must be one of Poisson, NB, ZIP, ZINB.", call. = FALSE)
    }
    if (!is.null(gamma_prior) &&
        (!is.numeric(gamma_prior) || length(gamma_prior) != 1L ||
            !is.finite(gamma_prior) || gamma_prior <= 0 || gamma_prior >= 1)) {
        stop("'gamma_prior' must be a single value strictly between 0 and 1.",
            call. = FALSE
        )
    }
    if (!is.null(theta_start) &&
        (!is.numeric(theta_start) || length(theta_start) != 1L ||
            !is.finite(theta_start) || theta_start <= 0 || theta_start >= 1)) {
        stop("'theta_start' must be a single value strictly between 0 and 1.",
            call. = FALSE
        )
    }
    if (!is.null(size_start) &&
        (!is.numeric(size_start) || any(!is.finite(size_start)) ||
            any(size_start <= 0))) {
        stop("'size_start' must be positive and finite.", call. = FALSE)
    }
    invisible(TRUE)
}

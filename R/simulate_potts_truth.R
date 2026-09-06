#' Simulate known-gamma HiCPotts truth from the fitted Potts model
#'
#' Generates a three-state latent field from the same four-neighbour Potts
#' auxiliary-field simulator used by the gamma updates, then generates counts
#' from the documented HiCPotts emission model. Component 1 is low-mean
#' zero-inflated noise, component 3 is elevated noise with approximately
#' component-1 slopes, and component 2 is unrestricted true signal.
#'
#' Unlike \code{simulate_hicpotts_truth()}, which deliberately places spatial
#' blocks for allocation stress tests, this function has a known generating
#' gamma and is therefore suitable for gamma bias and coverage experiments.
#'
#' @param N Lattice size.
#' @param gamma Potts coupling in the sampler's supported open interval (0, 1).
#' @param potts_sweeps Positive equilibration sweeps. The default scales with
#'   lattice area and is deliberately more conservative than one ABC proposal.
#' @param intercepts Length-three vector of log baseline intensities.
#' @param noise_slopes Four component-1 slopes.
#' @param signal_slopes Four unrestricted component-2 slopes.
#' @param false_signal_deviation Four deviations added to the component-1
#'   slopes to obtain component 3.
#' @param theta Structural zero-inflation probability for component 1.
#' @param size Length-three NB2 dispersion vector, or \code{NULL} for Poisson.
#' @param seed RNG seed.
#' @return A known-truth list compatible with
#'   \code{validate_hicpotts_simulation()}, with generating gamma and spatial
#'   model recorded in \code{settings}.
#' @examples
#' truth <- simulate_hicpotts_potts_truth(
#'     N = 8, gamma = 0.4,
#'     potts_sweeps = 50, seed = 11
#' )
#' table(truth$z_true)
#' @noRd
simulate_hicpotts_potts_truth <- function(
    N = 24L, gamma = 0.3, potts_sweeps = NULL,
    intercepts = c(log(1), log(14), log(20)),
    noise_slopes = c(-0.15, 0.25, -0.20, 0.15),
    signal_slopes = c(-1.00, 1.00, 0.80, -0.80),
    false_signal_deviation = rep(0, 4L),
    theta = 0.3, size = NULL, seed = 1L
) {
    N <- as.integer(N)
    if (length(N) != 1L || is.na(N) || N < 3L) {
        stop("'N' must be one integer of at least three.")
    }
    if (length(gamma) != 1L || !is.finite(gamma) || gamma <= 0 || gamma >= 1) {
        stop("'gamma' must lie strictly between zero and one.")
    }
    if (is.null(potts_sweeps)) {
        potts_sweeps <- max(500L, as.integer(round(500 * N * N / 400)))
    }
    potts_sweeps <- as.integer(potts_sweeps)
    if (length(potts_sweeps) != 1L || is.na(
        potts_sweeps
    ) || potts_sweeps < 1L) {
        stop("'potts_sweeps' must be one positive integer.")
    }
    if (length(intercepts) != 3L || any(!is.finite(intercepts))) {
        stop("'intercepts' must contain three finite log means.")
    }
    for (argument in c(
        "noise_slopes", "signal_slopes",
        "false_signal_deviation"
    )) {
        value <- get(argument, inherits = FALSE)
        if (length(value) != 4L || any(!is.finite(value))) {
            stop(sprintf(
                "'%s' must contain four finite coefficients.",
                argument
            ))
        }
    }
    if (length(theta) != 1L || !is.finite(theta) || theta < 0 || theta >= 1) {
        stop("'theta' must lie in [0, 1).")
    }
    if (!is.null(size) &&
        (length(size) != 3L || any(!is.finite(size)) || any(size <= 0))) {
        stop("'size' must be NULL or three positive finite values.")
    }

    withr::local_seed(as.integer(seed))
    z <- .hicpotts_simulate_potts_labels_cpp(N, gamma, potts_sweeps)
    storage.mode(z) <- "integer"

    distance <- abs(row(matrix(0, N, N)) - col(matrix(0, N, N)))
    gc <- matrix(stats::runif(N * N), N, N)
    tes <- matrix(stats::runif(N * N), N, N)
    acc <- matrix(stats::runif(N * N), N, N)
    design <- cbind(
        log1p(as.numeric(distance)), log1p(as.numeric(gc)),
        log1p(as.numeric(tes)), log1p(as.numeric(acc))
    )
    beta <- rbind(
        c(intercepts[1L], noise_slopes),
        c(intercepts[2L], signal_slopes),
        c(intercepts[3L], noise_slopes + false_signal_deviation)
    )
    eta_by_component <- vapply(seq_len(3L), function(k) {
        beta[k, 1L] + design %*% beta[k, -1L]
    }, numeric(N * N))
    eta <- eta_by_component[cbind(seq_len(N * N), as.integer(z))]
    mu <- exp(pmin(pmax(eta, -30), 30))

    y <- integer(N * N)
    for (idx in seq_along(y)) {
        component <- z[idx]
        y[idx] <- if (is.null(size)) {
            stats::rpois(1L, mu[idx])
        } else {
            stats::rnbinom(1L, size = size[component], mu = mu[idx])
        }
        if (component == 1L && stats::runif(1L) < theta) y[idx] <- 0L
    }
    dim(y) <- c(N, N)

    x_vars <- list(
        distance = list(distance), GC = list(gc),
        TES = list(tes), ACC = list(acc)
    )
    covariate_sds <- vapply(x_vars, function(x) {
        stats::sd(log1p(as.numeric(x[[1L]])))
    }, numeric(1))
    empirical_proportions <- tabulate(as.integer(z), nbins = 3L) / (N * N)
    list(
        y = y, z_true = z, x_vars = x_vars, beta_true = beta,
        settings = list(
            N = N, spatial_model = "three-state four-neighbour Potts",
            gamma = gamma, potts_sweeps = potts_sweeps,
            proportions = empirical_proportions, intercepts = intercepts,
            noise_slopes = noise_slopes, signal_slopes = signal_slopes,
            false_signal_deviation = false_signal_deviation,
            covariate_sds = covariate_sds, theta = theta, size = size,
            seed = seed
        )
    )
}

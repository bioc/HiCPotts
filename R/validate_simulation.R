#' Simulate a Hi-C lattice with biologically matched three-state truth
#'
#' Generates the component definition used by HiCPotts 1.3.11:
#' \describe{
#'   \item{Component 1, noise}{Low mean with structural zero inflation.}
#'   \item{Component 2, true signal}{Elevated mean and an unrestricted,
#'     different covariate-response pattern.}
#'   \item{Component 3, false signal}{Elevated noise with approximately the
#'     same covariate-response slopes as component 1.}
#' }
#'
#' Thus components 1 and 3 share a noise process in their slopes but differ in
#' baseline mean. Component 2 is separated by its covariate-response pattern,
#' not by a fixed component-2/3 brightness order.
#'
#' @param N Lattice size.
#' @param proportions Length-three vector of state proportions.
#' @param intercepts Length-three vector of log baseline intensities.
#' @param noise_slopes Four coefficients shared by components 1 and 3.
#' @param signal_slopes Four unrestricted component-2 coefficients.
#' @param false_signal_deviation Four deviations added to the shared noise
#'   slopes for component 3. The default zero gives exact equality; non-zero
#'   values test the intended approximate relationship.
#' @param theta Structural zero-inflation probability for component 1.
#' @param size Length-three dispersion vector, or \code{NULL} for Poisson.
#' @param seed RNG seed.
#' @return A list with \code{y}, \code{x_vars}, \code{z_true}, coefficient
#'   truth and generating settings.
#' @noRd
simulate_hicpotts_truth <- function(
    N = 24L,
    proportions = c(0.80, 0.12, 0.08),
    intercepts = c(log(1), log(14), log(20)),
    noise_slopes = c(-0.15, 0.25, -0.20, 0.15),
    signal_slopes = c(-1.00, 1.00, 0.80, -0.80),
    false_signal_deviation = rep(0, 4L),
    theta = 0.3,
    size = NULL,
    seed = 1L
) {
    N <- as.integer(N)
    if (length(N) != 1L || is.na(N) || N < 3L) {
        stop("'N' must be one integer of at least three.")
    }
    if (length(proportions) != 3L || any(!is.finite(proportions)) ||
        any(proportions < 0) || abs(sum(proportions) - 1) > 1e-8) {
        stop("'proportions' must be three non-negative values summing to one.")
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

    z <- matrix(1L, N, N)
    n_blocks <- max(1L, round(N / 4))
    place <- function(field, target, label) {
        placed <- 0L
        guard <- 0L
        while (placed < target && guard < 10000L) {
            guard <- guard + 1L
            width <- sample.int(max(2L, n_blocks), 1L) + 1L
            i0 <- sample.int(N, 1L)
            j0 <- sample.int(N, 1L)
            ii <- i0:min(N, i0 + width - 1L)
            jj <- j0:min(N, j0 + width - 1L)
            block <- field[ii, jj, drop = FALSE]
            free <- block == 1L
            if (!any(free)) next
            block[free] <- label
            field[ii, jj] <- block
            placed <- sum(field == label)
        }
        field
    }
    z <- place(z, round(proportions[2L] * N * N), 2L)
    z <- place(z, round(proportions[3L] * N * N), 3L)

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
    linear_index <- cbind(seq_len(N * N), as.integer(z))
    eta <- eta_by_component[linear_index]
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
    list(
        y = y, z_true = z, x_vars = x_vars, beta_true = beta,
        settings = list(
            N = N, proportions = proportions, intercepts = intercepts,
            noise_slopes = noise_slopes, signal_slopes = signal_slopes,
            false_signal_deviation = false_signal_deviation,
            covariate_sds = covariate_sds, theta = theta, size = size,
            seed = seed
        )
    )
}

#' Validate HiCPotts against a biologically matched simulation
#'
#' Fits the model to known truth and reports parameter recovery, three-class
#' accuracy, probability calibration, component-1/3 relationship recovery and
#' stability across chains.
#'
#' @param truth Output of \code{simulate_hicpotts_truth()}.
#' @param iterations MCMC iterations per chain.
#' @param n_chains Number of independent chains.
#' @param dist Emission family.
#' @param seed Base RNG seed.
#' @param min_draws Passed to \code{classify_hicpotts()}.
#' @param initialization Latent-state initialization methods cycled across
#'   chains. The defaults deliberately use distinct data-informed starts.
#' @param mc_cores Positive number of independent validation chains to execute
#'   concurrently.
#' @param ... Additional sampler arguments.
#' @return An object of class \code{hicpotts_validation}.
#' @noRd
validate_hicpotts_simulation <- function(
    truth, iterations = 2000L, n_chains = 2L, dist = "ZIP",
    seed = 1L, min_draws = 50L,
    initialization = c("likelihood_informed", "distance_adjusted"),
    mc_cores = 1L, ...
) {
    if (!is.list(truth) || is.null(truth$y) || is.null(truth$z_true) ||
        is.null(truth$beta_true)) {
        stop("'truth' must come from simulate_hicpotts_truth().")
    }
    n_chains <- as.integer(n_chains)
    if (length(n_chains) != 1L || is.na(n_chains) || n_chains < 1L) {
        stop("'n_chains' must be at least one.")
    }
    extra <- list(...)
    initialization <- match.arg(
        initialization,
        c(
            "likelihood_informed", "count_quantile", "distance_adjusted",
            "noise_anchored_random", "random"
        ),
        several.ok = TRUE
    )
    initialization <- rep(initialization, length.out = n_chains)
    negative_binomial <- dist %in% c("NB", "ZINB")
    default_size <- if (negative_binomial) {
        if (is.null(truth$settings$size)) rep(5, 3L) else truth$settings$size
    } else {
        NULL
    }

    tasks <- lapply(seq_len(n_chains), function(k) {
        chain_seed <- as.integer(seed + k)
        arguments <- c(list(
            N = nrow(truth$y), gamma_prior = 0.4, iterations = iterations,
            x_vars = truth$x_vars, y = truth$y, use_data_priors = TRUE,
            dist = dist,
            theta_start = if (dist %in% c("ZIP", "ZINB")) {
                truth$settings$theta
            } else {
                NULL
            },
            size_start = default_size
        ), extra)
        arguments <- arguments[!duplicated(names(arguments), fromLast = TRUE)]
        list(
            arguments = arguments, y = truth$y, x_vars = truth$x_vars,
            seed = chain_seed, initialization = initialization[k]
        )
    })
    fits <- .hicpotts_parallel_tasks(tasks, mc_cores = mc_cores)
    canonical_fits <- relabel_hicpotts(fits)

    cls <- classify_hicpotts(canonical_fits,
        relabel = FALSE,
        min_draws = min_draws
    )
    z_true <- as.integer(truth$z_true)
    predicted <- as.integer(cls$map_component)
    confusion <- table(
        factor(z_true, levels = seq_len(3L)),
        factor(predicted, levels = seq_len(3L)),
        dnn = c("true", "predicted")
    )
    per_class <- vapply(seq_len(3L), function(k) {
        tp <- sum(predicted == k & z_true == k)
        tn <- sum(predicted != k & z_true != k)
        sensitivity <- if (sum(z_true == k)) {
            tp / sum(z_true == k)
        } else {
            NA_real_
        }
        specificity <- if (sum(z_true != k)) {
            tn / sum(z_true != k)
        } else {
            NA_real_
        }
        precision <- if (sum(predicted == k)) {
            tp / sum(predicted == k)
        } else {
            NA_real_
        }
        c(
            sensitivity = sensitivity, specificity = specificity,
            precision = precision,
            f1 = if (is.finite(sensitivity) && is.finite(precision) &&
                sensitivity + precision > 0) {
                2 * sensitivity * precision / (sensitivity + precision)
            } else {
                NA_real_
            },
            balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE)
        )
    }, numeric(5))
    colnames(per_class) <- c("noise", "signal", "false signal")

    parameter_recovery <- do.call(rbind, lapply(
        seq_along(canonical_fits),
        function(k) {
            fit <- canonical_fits[[k]]
            do.call(rbind, lapply(seq_len(3L), function(component) {
                draws <- fit$chains[[component]]
                keep <- (floor(nrow(draws) / 2) + 1L):nrow(draws)
                retained <- draws[keep, , drop = FALSE]
                interval <- apply(retained, 2L, stats::quantile,
                    probs = c(0.025, 0.975), names = FALSE
                )
                truth_value <- as.numeric(truth$beta_true[component, ])
                estimate <- colMeans(retained)
                data.frame(
                    chain = k, component = component,
                    component_label = hicpotts_component_definition()$label[
                        component
                    ],
                    parameter = c("intercept", "distance", "GC", "TES", "ACC"),
                    truth = truth_value, estimate = estimate,
                    lower = interval[1L, ], upper = interval[2L, ],
                    bias = estimate - truth_value,
                    covered = interval[1L,
                    ] <= truth_value & truth_value <= interval[2L, ],
                    stringsAsFactors = FALSE
                )
            }))
        }
    ))

    dispersion_recovery <- NULL
    if (!is.null(truth$settings$size) && negative_binomial) {
        dispersion_recovery <- do.call(rbind, lapply(
            seq_along(canonical_fits),
            function(k) {
                draws <- canonical_fits[[k]]$size
                keep <- (floor(ncol(draws) / 2) + 1L):ncol(draws)
                retained <- draws[, keep, drop = FALSE]
                interval <- t(apply(retained, 1L, stats::quantile,
                    probs = c(0.025, 0.975), names = FALSE
                ))
                estimate <- rowMeans(retained)
                data.frame(
                    chain = k, component = seq_len(3L),
                    component_label = hicpotts_component_definition()$label,
                    truth = truth$settings$size, estimate = estimate,
                    lower = interval[, 1L], upper = interval[, 2L],
                    bias = estimate - truth$settings$size,
                    covered = interval[, 1L] <= truth$settings$size &
                        truth$settings$size <= interval[, 2L],
                    stringsAsFactors = FALSE
                )
            }
        ))
    }

    theta_recovery <- NULL
    if (dist %in% c("ZIP", "ZINB")) {
        theta_recovery <- do.call(rbind, lapply(
            seq_along(canonical_fits),
            function(k) {
                draws <- canonical_fits[[k]]$theta
                keep <- (floor(length(draws) / 2) + 1L):length(draws)
                retained <- draws[keep]
                interval <- stats::quantile(retained, c(0.025, 0.975),
                    names = FALSE
                )
                data.frame(
                    chain = k, truth = truth$settings$theta,
                    estimate = mean(retained), lower = interval[1L],
                    upper = interval[2L],
                    bias = mean(retained) - truth$settings$theta,
                    covered = interval[1L] <= truth$settings$theta &&
                        truth$settings$theta <= interval[2L],
                    stringsAsFactors = FALSE
                )
            }
        ))
    }

    p_assigned <- cls$class_probability
    correct <- predicted == z_true
    bins <- cut(p_assigned, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
    calibration <- do.call(rbind, lapply(levels(bins), function(bin) {
        keep <- bins == bin
        if (!any(keep)) {
            return(NULL)
        }
        data.frame(
            bin = bin, n = sum(keep),
            mean_predicted = mean(p_assigned[keep]),
            observed_accuracy = mean(correct[keep]),
            stringsAsFactors = FALSE
        )
    }))
    ece <- if (is.null(calibration)) {
        NA_real_
    } else {
        sum(calibration$n / sum(calibration$n) *
            abs(calibration$mean_predicted - calibration$observed_accuracy))
    }

    relationship_recovery <- do.call(rbind, lapply(
        seq_along(canonical_fits),
        function(k) {
            relabelled <- canonical_fits[[k]]
            n_draw <- nrow(relabelled$chains[[1L]])
            keep <- (floor(n_draw / 2) + 1L):n_draw
            sds <- attr(relabelled$chains[[1L]], "proposal_covariate_sds",
                exact = TRUE
            )
            if (is.null(sds) || length(sds) < 4L) sds <- rep(1, 4L)
            slope_columns <- seq.int(2L, 5L)
            slope_indices <- seq_len(4L)
            b1 <- relabelled$chains[[1L]][keep, slope_columns, drop = FALSE]
            b2 <- relabelled$chains[[2L]][keep, slope_columns, drop = FALSE]
            b3 <- relabelled$chains[[3L]][keep, slope_columns, drop = FALSE]
            d12 <- sqrt(rowSums(sweep(b2 - b1, 2L, sds[slope_indices], "*")^2))
            d13 <- sqrt(rowSums(sweep(b3 - b1, 2L, sds[slope_indices], "*")^2))
            intercept_gap13 <- mean(relabelled$chains[[3L]][keep, 1L] -
                relabelled$chains[[1L]][keep, 1L])
            data.frame(
                chain = k,
                distance_signal_to_noise = mean(d12),
                distance_false_signal_to_noise = mean(d13),
                false_signal_intercept_gap = intercept_gap13,
                correct_relationship = mean(d13) < mean(
                    d12
                ) && intercept_gap13 > 0,
                relationship_probability =
                    relabelled$noise_relationship_probability,
                threshold_policy =
                    paste0("none in package; user chooses downstream ",
                        "probability threshold"),
                stringsAsFactors = FALSE
            )
        }
    ))

    stability <- if (n_chains > 1L) {
        allocation_diagnostics(canonical_fits)$between_chain_disagreement
    } else {
        NULL
    }
    gamma_recovery <- NULL
    if (!is.null(truth$settings$gamma) &&
        length(truth$settings$gamma) == 1L &&
        is.finite(truth$settings$gamma)) {
        gamma_recovery <- do.call(rbind, lapply(
            seq_along(canonical_fits),
            function(k) {
                draws <- canonical_fits[[k]]$gamma
                keep <- (floor(length(draws) / 2) + 1L):length(draws)
                retained <- draws[keep]
                interval <- stats::quantile(retained, c(0.025, 0.975),
                    names = FALSE
                )
                data.frame(
                    chain = k, truth = truth$settings$gamma,
                    estimate = mean(retained), posterior_sd = stats::sd(
                        retained
                    ),
                    lower = interval[1L], upper = interval[2L],
                    bias = mean(retained) - truth$settings$gamma,
                    covered = interval[1L] <= truth$settings$gamma &&
                        truth$settings$gamma <= interval[2L],
                    acceptance_rate = attr(draws, "gamma_acceptance_rate",
                        exact = TRUE
                    ),
                    stringsAsFactors = FALSE
                )
            }
        ))
    }
    structure(
        list(
            accuracy = mean(correct), confusion = confusion,
            per_class = per_class,
            parameter_recovery = parameter_recovery,
            dispersion_recovery = dispersion_recovery,
            theta_recovery = theta_recovery,
            calibration = calibration, expected_calibration_error = ece,
            relationship_recovery = relationship_recovery,
            stability = stability,
            gamma_recovery = gamma_recovery,
            fits = canonical_fits,
            n_chains = n_chains, iterations = iterations
        ),
        class = "hicpotts_validation"
    )
}

#' @export
print.hicpotts_validation <- function(x, ...) {
    cat("HiCPotts known-truth validation\n")
    cat(sprintf("  chains: %d   iterations: %d\n", x$n_chains, x$iterations))
    cat(sprintf("  three-way accuracy: %.3f\n", x$accuracy))
    cat(sprintf(
        "  expected calibration error: %.4f\n",
        x$expected_calibration_error
    ))
    cat("  per-class sensitivity / specificity / precision:\n")
    print(round(x$per_class, 3))
    cat(sprintf(
        "  beta interval inclusion: %.3f\n",
        mean(x$parameter_recovery$covered)
    ))
    if (!is.null(x$dispersion_recovery)) {
        cat(sprintf(
            "  dispersion interval inclusion: %.3f\n",
            mean(x$dispersion_recovery$covered)
        ))
    }
    cat(sprintf(
        "  component-1/3 relationship recovered in %d/%d chains\n",
        sum(x$relationship_recovery$correct_relationship),
        nrow(x$relationship_recovery)
    ))
    if (!is.null(x$stability)) {
        cat(sprintf(
            "  between-chain disagreement: mean %.4f, max %.4f\n",
            x$stability$mean, x$stability$max
        ))
    }
    if (!is.null(x$gamma_recovery)) {
        cat(sprintf(
            "  gamma truth %.3f; chain means %s\n",
            x$gamma_recovery$truth[1L],
            paste(format(round(x$gamma_recovery$estimate, 3), nsmall = 3),
                collapse = ", "
            )
        ))
    }
    invisible(x)
}

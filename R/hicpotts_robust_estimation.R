#' Diagnose HiCPotts covariate identifiability
#'
#' Examines the four `log1p` covariates used by HiCPotts before fitting. It
#' reports their correlation matrix, pairwise high-correlation flags and the
#' condition number of the standardized design matrix. These checks do not
#' modify the model or the supplied covariates.
#'
#' @param x_vars Named covariate list containing distance, GC, TES and ACC.
#' @param dataset_index Dataset index when covariates are stored as lists.
#' @param correlation_threshold Absolute-correlation warning threshold.
#' @param condition_threshold Design-matrix condition-number warning threshold.
#' @param minimum_sd Threshold used to identify effectively constant covariates.
#' @return A list containing the correlation matrix, pairwise diagnostics,
#'   condition number, standard deviations and warnings.
#' @examples
#' N <- 6
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' diagnostics <- diagnose_hicpotts_covariates(x_vars)
#' diagnostics$condition_number
#' @noRd
diagnose_hicpotts_covariates <- function(
    x_vars, dataset_index = 1L, correlation_threshold = 0.8,
    condition_threshold = 30, minimum_sd = 1e-8
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
    if (!is.finite(correlation_threshold) || correlation_threshold <= 0 ||
        correlation_threshold >= 1) {
        stop("'correlation_threshold' must be strictly between zero and one.")
    }
    if (!is.finite(condition_threshold) || condition_threshold <= 1) {
        stop("'condition_threshold' must be greater than one.")
    }

    select_matrix <- function(name) {
        item <- x_vars[[name]]
        if (is.list(item)) {
            if (length(item) < dataset_index) {
                stop("Covariate '", name, "' does not contain dataset_index.")
            }
            item <- item[[dataset_index]]
        }
        if (!is.matrix(item) || !is.numeric(item) || any(!is.finite(item)) ||
            any(item <= -1)) {
            stop(
                "Each selected covariate must be a finite numeric matrix ",
                "greater than -1."
            )
        }
        item
    }
    matrices <- lapply(required, select_matrix)
    dimensions <- vapply(
        matrices, function(x) paste(dim(x), collapse = "x"),
        character(1)
    )
    if (length(unique(dimensions)) != 1L) {
        stop("All selected covariate matrices must have identical dimensions.")
    }
    design <- vapply(matrices, function(x) log1p(as.numeric(x)), numeric(length(
        matrices[[1]]
    )))
    colnames(design) <- required
    covariate_sd <- apply(design, 2L, stats::sd)
    constant <- !is.finite(covariate_sd) | covariate_sd < minimum_sd
    correlation <- matrix(NA_real_, ncol(design), ncol(design),
        dimnames = list(colnames(design), colnames(design))
    )
    if (any(!constant)) {
        correlation[!constant, !constant] <- stats::cor(design[, !constant,
            drop = FALSE
        ])
    }
    pair_index <- which(upper.tri(correlation), arr.ind = TRUE)
    pairs <- data.frame(
        covariate1 = rownames(correlation)[pair_index[, 1L]],
        covariate2 = colnames(correlation)[pair_index[, 2L]],
        correlation = correlation[pair_index], stringsAsFactors = FALSE
    )
    pairs$high_correlation <- abs(pairs$correlation) >= correlation_threshold
    condition_number <- Inf
    if (!any(constant)) {
        singular_values <- svd(scale(design), nu = 0L, nv = 0L)$d
        if (length(singular_values) && min(singular_values) > 0) {
            condition_number <- max(singular_values) / min(singular_values)
        }
    }
    notes <- character()
    if (any(constant)) {
        notes <- c(notes, paste(
            "Effectively constant covariates:",
            paste(required[constant], collapse = ", ")
        ))
    }
    if (any(pairs$high_correlation %in% TRUE, na.rm = TRUE)) {
        notes <- c(
            notes,
            "At least one covariate pair exceeds the correlation threshold."
        )
    }
    if (!is.finite(
        condition_number
    ) || condition_number > condition_threshold) {
        notes <- c(
            notes,
            "The standardized covariate design is poorly conditioned."
        )
    }
    list(
        correlation = correlation, pairs = pairs,
        condition_number = condition_number,
        log1p_covariate_sd = stats::setNames(covariate_sd, required),
        thresholds = c(
            correlation = correlation_threshold,
            condition_number = condition_threshold
        ),
        warnings = unique(notes), reliable = !length(notes)
    )
}

.hicpotts_likelihood_initial_z <- function(
    y, x_vars, dist, size_start, nstart = 6L, max_iter = 15L
) {
    required <- c("distance", "GC", "TES", "ACC")
    if (!is.list(x_vars) || !all(required %in% names(x_vars))) {
        stop("Likelihood-informed initialization requires all four covariates.")
    }
    select_covariate <- function(name) {
        value <- x_vars[[name]]
        if (is.list(value)) value <- value[[1L]]
        if (!is.matrix(value) || !identical(dim(value), dim(y)) ||
            any(!is.finite(value)) || any(value <= -1)) {
            stop(
                "Each initialization covariate must be finite, match y and ",
                "exceed -1."
            )
        }
        log1p(as.numeric(value))
    }
    response <- as.numeric(y)
    design <- cbind(1, vapply(required, select_covariate, numeric(length(y))))
    n <- length(response)
    nstart <- max(1L, as.integer(nstart))
    max_iter <- max(1L, as.integer(max_iter))
    use_nb <- dist %in% c("NB", "ZINB")
    use_zi <- dist %in% c("ZIP", "ZINB")
    if (use_nb) {
        if (is.null(size_start) || length(size_start) != 3L ||
            any(!is.finite(size_start)) || any(size_start <= 0)) {
            stop(
                "Likelihood-informed NB/ZINB initialization requires three ",
                "positive sizes."
            )
        }
        size_start <- as.numeric(size_start)
    }
    count_start <- as.integer(pmin(3L, pmax(
        1L,
        ceiling(3 * rank(log1p(response), ties.method = "random") / n)
    )))
    starts <- vector("list", nstart)
    starts[[1L]] <- count_start
    if (nstart > 1L) {
        for (s in 2:nstart) {
            low <- rank(log1p(response), ties.method = "random") <= ceiling(
                n / 3
            )
            candidate <- rep.int(1L, n)
            candidate[!low] <- sample(rep(c(2L, 3L), length.out = sum(!low)))
            starts[[s]] <- candidate
        }
    }
    best_score <- -Inf
    best_labels <- count_start
    for (s in seq_along(starts)) {
        labels <- starts[[s]]
        log_probability <- matrix(-Inf, n, 3L)
        failed_fit <- FALSE
        for (iteration in seq_len(max_iter)) {
            for (component in seq_len(3)) {
                index <- which(labels == component)
                if (length(index) < ncol(design) + 2L) {
                    labels <- count_start
                    index <- which(labels == component)
                }
                fit <- tryCatch(
                    withCallingHandlers(
                        stats::glm.fit(
                            design[index, , drop = FALSE], response[index],
                            family = stats::poisson()
                        ),
                        warning = function(w) invokeRestart("muffleWarning")
                    ),
                    error = function(e) NULL
                )
                if (is.null(fit) ||
                    length(fit$coefficients) != ncol(design)) {
                    failed_fit <- TRUE
                    break
                }
                beta <- fit$coefficients
                beta[!is.finite(beta)] <- 0
                eta <- pmax(-20, pmin(20, as.numeric(design %*% beta)))
                mu <- exp(eta)
                if (use_nb) {
                    log_density <- stats::dnbinom(response,
                        mu = mu,
                        size = size_start[component], log = TRUE
                    )
                    p_zero <- stats::dnbinom(0, mu = mu, size = size_start[
                        component
                    ])
                } else {
                    log_density <- stats::dpois(response,
                        lambda = mu,
                        log = TRUE
                    )
                    p_zero <- stats::dpois(0, lambda = mu)
                }
                if (component == 1L && use_zi) {
                    observed_zero <- mean(response[index] == 0)
                    model_zero <- mean(p_zero[index])
                    theta <- pmax(0.01, pmin(
                        0.95,
                        (observed_zero - model_zero) / pmax(
                            1e-8,
                            1 - model_zero
                        )
                    ))
                    zero <- response == 0
                    log_density[zero] <- log(theta + (1 - theta) * p_zero[zero])
                    log_density[!zero] <- log1p(-theta) + log_density[!zero]
                }
                weight <- max(length(index) / n, 1 / n)
                log_probability[, component] <- log(weight) + log_density
            }
            if (failed_fit) break
            jitter <- matrix(stats::runif(3L * n, -1e-10, 1e-10), n, 3L)
            updated <- max.col(log_probability + jitter, ties.method = "first")
            if (identical(updated, labels)) break
            if (any(tabulate(updated, nbins = 3L) < ncol(design) + 2L)) {
                labels <- count_start
                break
            }
            labels <- updated
        }
        if (failed_fit) next
        score <- sum(log_probability[cbind(seq_len(n), labels)])
        if (is.finite(score) && score > best_score) {
            best_score <- score
            best_labels <- labels
        }
    }
    zero_fraction <- vapply(seq_len(3), function(k) {
        mean(response[best_labels == k] == 0)
    }, numeric(1))
    mean_count <- vapply(seq_len(3), function(k) {
        mean(log1p(response[best_labels == k]))
    }, numeric(1))
    component1 <- order(-zero_fraction, mean_count)[1L]
    remaining <- setdiff(seq_len(3), component1)
    remaining <- remaining[order(mean_count[remaining])]
    map <- integer(3L)
    map[component1] <- 1L
    map[remaining] <- 2:3
    matrix(map[best_labels], nrow(y), ncol(y))
}

#' Construct an initial latent-state allocation
#'
#' Creates an N-by-N three-state starting allocation for independently started
#' chains. This affects computation and mode exploration only; the Gibbs target
#' and the fitted generative model are unchanged. Observed zeros are not forced
#' to component 1 in these parameter-inference initializations.
#'
#' @param y Numeric square interaction-count matrix.
#' @param x_vars Covariates required for `distance_adjusted` initialization.
#' @param method One of `likelihood_informed`, `random`, `count_quantile`,
#'   `distance_adjusted` or `noise_anchored_random`. The likelihood-informed
#'   method fits a non-spatial hard mixture under the requested count family.
#' @param seed Optional integer seed.
#' @param dist Count family used by likelihood-informed initialization.
#' @param size_start Optional length-three NB2 size starting values.
#' @return An integer matrix with labels 1, 2 and 3.
#' @examples
#' N <- 6
#' y <- matrix(rpois(N * N, 3), N, N)
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' z0 <- make_hicpotts_initial_z(y, x_vars = x_vars, method = "count_quantile")
#' table(z0)
#' @noRd
make_hicpotts_initial_z <- function(
    y, x_vars = NULL,
    method = c(
        "likelihood_informed", "count_quantile", "distance_adjusted",
        "noise_anchored_random", "random"
    ), seed = NULL,
    dist = "ZIP", size_start = NULL
) {
    method <- match.arg(method)
    dist <- match.arg(dist, c("Poisson", "NB", "ZIP", "ZINB"))
    ## M3: the shared count contract. This wrapper previously checked only
    ## non-negativity and finiteness, so it accepted continuous normalised
    ## values that make the Poisson/NB likelihood -Inf.
    .hicpotts_validate_counts(y, what = "y")
    .hicpotts_validate_starts(size_start = size_start, dist = dist)
    if (!is.null(seed)) {
        seed <- as.integer(seed)
        if (length(seed) != 1L || is.na(seed)) {
            stop(
                "'seed' must be one integer."
            )
        }
        withr::local_seed(seed)
    }
    n <- length(y)
    if (method == "likelihood_informed") {
        return(.hicpotts_likelihood_initial_z(y, x_vars, dist, size_start))
    }
    if (method == "random") {
        return(matrix(sample.int(3L, n, replace = TRUE), nrow(y), ncol(y)))
    }
    raw_score <- log1p(as.numeric(y))
    if (method == "count_quantile") {
        ranks <- rank(raw_score, ties.method = "random")
        labels <- pmin(3L, pmax(1L, ceiling(3 * ranks / n)))
        return(matrix(as.integer(labels), nrow(y), ncol(y)))
    }
    raw_ranks <- rank(raw_score, ties.method = "random")
    noise <- raw_ranks <= ceiling(n / 3)
    labels <- rep.int(1L, n)
    if (method == "distance_adjusted") {
        if (!is.list(x_vars) || is.null(x_vars$distance)) {
            stop(
                "'x_vars$distance' is required for distance_adjusted ",
                "initialization."
            )
        }
        distance <- if (is.list(x_vars$distance)) {
            x_vars$distance[[
                1L
            ]]
        } else {
            x_vars$distance
        }
        if (!is.matrix(distance) || !identical(dim(distance), dim(y)) ||
            any(!is.finite(distance)) || any(distance <= -1)) {
            stop(
                "The distance covariate must be a finite matrix matching y ",
                "and greater than -1."
            )
        }
        score <- stats::lm.fit(
            cbind(1, log1p(as.numeric(distance))),
            raw_score
        )$residuals
    }
    remaining <- which(!noise)
    if (method == "noise_anchored_random") {
        remaining <- sample(remaining)
        labels[remaining] <- rep(c(2L, 3L), length.out = length(remaining))
    } else {
        residual_ranks <- rank(score[remaining], ties.method = "random")
        labels[remaining] <- ifelse(
            residual_ranks <= ceiling(length(remaining) / 2), 2L, 3L
        )
    }
    matrix(as.integer(labels), nrow(y), ncol(y))
}

#' Run standardized-prior sensitivity fits
#'
#' Repeats independently initialized, relabelled chains under standardized
#' slope prior SDs of 0.25, 0.5 and 1 by default.
#'
#' @param slope_sds Positive standardized-scale slope prior SDs.
#' @param burnin Initial draws discarded by diagnostics; defaults to half.
#' @param minimum_component_cells Non-negative occupancy threshold. Zero
#'   permits empty components while retaining occupancy reporting.
#' @param minimum_ess Minimum recommended effective sample size.
#' @param maximum_rhat Maximum recommended R-hat.
#' @inheritParams run_hicpotts_chains
#' @return A named list containing priors, relabelled fits and diagnostics for
#'   each prior scale.
#' @examples
#' N <- 5
#' y <- matrix(rpois(N * N, 4), N, N)
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' sens <- run_hicpotts_prior_sensitivity(
#'     N = N, iterations = 5, x_vars = x_vars,
#'     y = y, dist = "ZINB", theta_start = 0.5, size_start = c(1, 1, 1),
#'     seeds = 1:2, slope_sds = c(0.5, 1)
#' )
#' names(sens)
#' @noRd
run_hicpotts_prior_sensitivity <- function(
    N, gamma_prior = 0.3, iterations, x_vars, y,
    theta_start = NULL, size_start = NULL, dist = "ZIP", epsilon = NULL,
    distance_metric = "manhattan", seeds = seq_len(4),
    initialization = c(
        "likelihood_informed", "likelihood_informed",
        "distance_adjusted", "noise_anchored_random"
    ),
    slope_sds = c(0.25, 0.5, 1), burnin = NULL,
    minimum_component_cells = 0L, minimum_ess = 200,
    maximum_rhat = 1.01
) {
    if (!length(slope_sds) || any(!is.finite(slope_sds)) || any(
        slope_sds <= 0
    )) {
        stop("'slope_sds' must contain finite positive values.")
    }
    covariates <- diagnose_hicpotts_covariates(x_vars)
    out <- lapply(slope_sds, function(prior_sd) {
        priors <- make_hicpotts_scaled_priors(
            x_vars,
            slope_sd_standardized = prior_sd
        )
        fits <- run_hicpotts_chains(
            N = N, gamma_prior = gamma_prior, iterations = iterations,
            x_vars = x_vars, y = y, theta_start = theta_start,
            size_start = size_start, use_data_priors = FALSE,
            user_fixed_priors = priors, dist = dist, epsilon = epsilon,
            distance_metric = distance_metric, seeds = seeds,
            initialization = initialization, relabel = TRUE
        )
        diagnostics <- diagnose_hicpotts_fit(
            fits,
            burnin = burnin,
            minimum_component_cells = minimum_component_cells,
            minimum_ess = minimum_ess, maximum_rhat = maximum_rhat,
            covariate_diagnostics = covariates, relabel = FALSE
        )
        list(
            prior_sd = prior_sd, priors = priors, fits = fits,
            diagnostics = diagnostics
        )
    })
    names(out) <- paste0("slope_sd_", format(slope_sds, trim = TRUE))
    out
}

.hicpotts_shared_pilot_priors <- function(pilot_fits, y, x_vars) {
    if (!is.list(pilot_fits) || !length(pilot_fits)) {
        stop("The empirical-Bayes pilot did not return any chains.")
    }
    ## Pool posterior membership probabilities only after applying the same
    ## biological relabelling rule to every pilot.  The resulting three weight
    ## matrices define one empirical-Bayes prior shared by every production
    ## chain, which restores the common-target assumption required by R-hat.
    pilot_fits <- relabel_hicpotts(pilot_fits)
    probabilities <- lapply(pilot_fits, function(fit) {
        value <- fit$z_probabilities
        if (!is.list(value) || length(value) != 3L ||
            !all(vapply(value, is.matrix, logical(1))) ||
            !all(vapply(value, function(z) identical(dim(z), dim(y)),
                logical(1)))) {
            stop("Every pilot chain must contain three N-by-N probabilities.")
        }
        value
    })
    weights <- lapply(seq_len(3L), function(component) {
        Reduce(`+`, lapply(probabilities, `[[`, component)) /
            length(probabilities)
    })
    total <- Reduce(`+`, weights)
    if (any(!is.finite(total)) || any(total <= 0)) {
        stop("Pooled pilot membership probabilities are invalid.")
    }
    weights <- lapply(weights, function(value) value / total)

    selected_covariate <- function(name) {
        value <- x_vars[[name]]
        if (is.list(value) && !is.matrix(value)) value <- value[[1L]]
        if (!is.matrix(value) || !identical(dim(value), dim(y))) {
            stop("Every pilot covariate must be an N-by-N matrix.")
        }
        value
    }
    covariates <- lapply(
        c("distance", "GC", "TES", "ACC"), selected_covariate
    )
    priors <- .hicpotts_soft_eb_priors_cpp(
        y, covariates[[1L]], covariates[[2L]], covariates[[3L]],
        covariates[[4L]], weights[[1L]], weights[[2L]], weights[[3L]]
    )
    attr(priors, "method") <-
        "shared soft-allocation empirical Bayes pilot"
    attr(priors, "shared_across_production_chains") <- TRUE
    attr(priors, "pilot_chains") <- length(pilot_fits)
    attr(priors, "posterior_expected_cells") <- vapply(
        weights, function(value) sum(as.numeric(value)), numeric(1)
    )
    list(priors = priors, weights = weights, pilot_fits = pilot_fits)
}

.hicpotts_pilot_progress_message <- function(
    n_chains, pilot_iterations, production_iterations
) {
    sprintf(
        paste0(
            "Running %d shared-prior pilot chains for %d pilot iterations; ",
            "production will then run up to %d requested iterations per chain."
        ),
        n_chains, pilot_iterations, production_iterations
    )
}

#' Robust general-use HiCPotts estimation workflow
#'
#' Runs four independently seeded chains from varied latent-state starts, uses
#' soft-allocation empirical-Bayes priors by default, relabels before pooling
#' and returns convergence, occupancy and collinearity reliability flags. The
#' default 20,000 iterations and 5,000 burn-in are production recommendations,
#' not changes to the underlying statistical model.
#'
#' @param burnin Initial draws removed from summaries.
#' @param use_data_priors Whether to estimate one soft-allocation
#' empirical-Bayes regression prior in a multi-chain pilot and freeze that same
#' prior for every production chain. The default is \code{TRUE}.
#' @param slope_sd_standardized Standardized-scale slope prior SD used only when
#' \code{use_data_priors = FALSE}.
#' @param minimum_component_cells Non-negative occupancy threshold per
#' component. Zero permits empty components while retaining occupancy reports.
#' @param minimum_ess Minimum recommended effective sample size.
#' @param maximum_rhat Maximum recommended R-hat.
#' @param gamma_boundary_tolerance Distance from zero or one counted as the
#' gamma boundary region.
#' @param maximum_gamma_boundary_fraction Largest allowed retained gamma
#' fraction in that region in every chain.
#' @param minimum_gamma_unique Minimum distinct retained gamma draws required
#' per chain.
#' @param mode_screen Whether parameter summaries use the largest group of
#' chains with mutually coherent internal allocations.
#' @param mode_agreement_threshold Minimum allocation agreement connecting
#' chains.
#' @param minimum_mode_chains Minimum coherent chains needed before screening
#' excludes any chain. Every chain remains available in `all_fits`.
#' @inheritParams run_hicpotts_chains
#' @return A list containing relabelled fits, diagnostics, covariate checks,
#'   priors and settings.
#' @examples
#' N <- 6
#' y <- matrix(rpois(N * N, 4), N, N)
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' fit <- fit_hicpotts_robust(
#'     N = N, x_vars = x_vars, y = y, dist = "ZINB",
#'     iterations = 5, burnin = 2, seeds = 1:2, theta_start = 0.5,
#'     size_start = c(1, 1, 1)
#' )
#' names(fit)
#' @noRd
fit_hicpotts_robust <- function(
    N, x_vars, y, dist = "ZIP", gamma_prior = 0.3,
    iterations = 20000L, burnin = 5000L, theta_start = NULL,
    size_start = NULL, epsilon = NULL, distance_metric = "manhattan",
    seeds = seq_len(4),
    initialization = c(
        "likelihood_informed", "likelihood_informed",
        "distance_adjusted", "noise_anchored_random"
    ),
    use_data_priors = TRUE, user_fixed_priors = NULL,
    slope_sd_standardized = 0.5,
    minimum_component_cells = 0L,
    minimum_ess = 200, maximum_rhat = 1.01,
    gamma_boundary_tolerance = 0.01,
    maximum_gamma_boundary_fraction = 0.95,
    minimum_gamma_unique = 20L,
    mcse_stop = TRUE, mcse_min_iterations = 10000L,
    mcse_check_interval = 500L, mcse_relative_threshold = 0.05,
    tempering_warmup = 0L, tempering_beta_min = 0.30,
    tempering_cycle = 500L, mode_screen = TRUE,
    gamma_prior_shape1 = 1, gamma_prior_shape2 = 1,
    abc_epsilon_quantile = 0.10, gamma_update_interval = 5L,
    abc_potts_sweeps = 0L, abc_sim_reps = 4L,
    z_probability_burnin = NULL,
    mode_agreement_threshold = 0.8, minimum_mode_chains = 2L,
    use_noise_relationship_prior = TRUE,
    noise_link_sd = 0.5, noise_order_strength = 10,
    noise_order_width = 0.5,
    comp23_barrier_kappa = 10, comp23_barrier_w = 0.3,
    branch_swap_interval = 1L, signal_block_move_interval = 5L,
    gamma_large_jump_probability = 0.10,
    gamma_large_jump_multiplier = 4,
    gamma_independence_probability = 0.05,
    gamma_method = "abc", relabel = TRUE, mc_cores = 1L,
    verbose = FALSE, progress_interval = 50L
) {
    iterations <- as.integer(iterations)
    gamma_method <- match.arg(gamma_method)
    burnin <- as.integer(burnin)
    if (length(iterations) != 1L || is.na(iterations) || iterations < 1L) {
        stop("'iterations' must be one positive integer.")
    }
    if (length(burnin) != 1L || is.na(
        burnin
    ) || burnin < 0L || burnin >= iterations + 1L) {
        stop("'burnin' must be between zero and iterations.")
    }
    if (is.null(z_probability_burnin)) z_probability_burnin <- burnin
    if (iterations < 15000L) {
        warning(
            "Fewer than 15,000 iterations requested; rely on ESS and R-hat ",
            "before interpreting estimates."
        )
    }
    if (length(seeds) < 4L) {
        warning(
            "Fewer than four independent chains requested; four are ",
            "recommended for general use."
        )
    }
    if (!is.logical(use_data_priors) || length(use_data_priors) != 1L ||
        is.na(use_data_priors)) {
        stop("'use_data_priors' must be TRUE or FALSE.")
    }
    covariates <- diagnose_hicpotts_covariates(x_vars)
    run_phase <- function(
        phase_iterations, phase_z_burnin, phase_use_data_priors,
        phase_priors, phase_epsilon, phase_mcse_stop, phase_relabel
    ) {
        run_hicpotts_chains(
        N = N, gamma_prior = gamma_prior, iterations = phase_iterations,
        x_vars = x_vars, y = y, theta_start = theta_start,
        size_start = size_start,
        use_data_priors = phase_use_data_priors,
        user_fixed_priors = phase_priors,
        dist = dist, epsilon = phase_epsilon,
        distance_metric = distance_metric, seeds = seeds,
        initialization = initialization,
        mcse_stop = phase_mcse_stop,
        mcse_min_iterations = mcse_min_iterations,
        mcse_check_interval = mcse_check_interval,
        mcse_relative_threshold = mcse_relative_threshold,
        tempering_warmup = tempering_warmup,
        tempering_beta_min = tempering_beta_min,
        tempering_cycle = tempering_cycle,
        gamma_prior_shape1 = gamma_prior_shape1,
        gamma_prior_shape2 = gamma_prior_shape2,
        abc_epsilon_quantile = abc_epsilon_quantile,
        gamma_update_interval = gamma_update_interval,
        abc_potts_sweeps = abc_potts_sweeps,
        abc_sim_reps = abc_sim_reps,
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
        z_probability_burnin = phase_z_burnin,
        relabel = phase_relabel,
        mc_cores = mc_cores, verbose = verbose,
        progress_interval = progress_interval
        )
    }

    pilot_summary <- NULL
    production_epsilon <- epsilon
    if (isTRUE(use_data_priors)) {
        ## The pilot is capped at 5,000 iterations so a 15,000-iteration
        ## production analysis has the same maximum total work as the previous
        ## 20,000-iteration single-stage workflow.  Its only purpose is to
        ## estimate one shared prior; no pilot draws enter final summaries.
        pilot_iterations <- max(2L, min(5000L, max(2L, burnin)))
        pilot_burnin <- min(
            pilot_iterations - 1L, floor(pilot_iterations / 2L)
        )
        if (isTRUE(verbose)) {
            message(.hicpotts_pilot_progress_message(
                length(seeds), pilot_iterations, iterations
            ))
        }
        pilot_fits <- run_phase(
            phase_iterations = pilot_iterations,
            phase_z_burnin = pilot_burnin,
            phase_use_data_priors = TRUE,
            phase_priors = NULL,
            phase_epsilon = epsilon,
            phase_mcse_stop = FALSE,
            phase_relabel = TRUE
        )
        shared <- .hicpotts_shared_pilot_priors(pilot_fits, y, x_vars)
        priors <- shared$priors
        if (is.null(production_epsilon)) {
            production_epsilon <-
                shared$pilot_fits[[1L]]$sampler_settings$abc_epsilon
        }
        pilot_summary <- list(
            method = "multi-chain soft-allocation empirical-Bayes pilot",
            iterations = pilot_iterations,
            burnin = pilot_burnin,
            chains = length(shared$pilot_fits),
            seeds = as.integer(seeds),
            posterior_expected_cells = attr(
                priors, "posterior_expected_cells", exact = TRUE
            ),
            abc_epsilon = production_epsilon,
            elapsed_seconds = sum(vapply(
                shared$pilot_fits,
                function(fit) fit$performance$elapsed_seconds,
                numeric(1)
            )),
            excluded_from_posterior_summaries = TRUE
        )
        rm(pilot_fits, shared)
        if (isTRUE(verbose)) {
            message(
                "Shared prior frozen; starting production chains with a ",
                "common posterior target."
            )
        }
    } else if (!is.null(user_fixed_priors)) {
        priors <- user_fixed_priors
    } else {
        priors <- make_hicpotts_scaled_priors(
            x_vars,
            slope_sd_standardized = slope_sd_standardized
        )
    }

    ## Production chains always use one fixed prior object.  With empirical
    ## Bayes this is the pooled pilot result; with user/scaled priors it is the
    ## supplied deterministic object.  Consequently every chain entering
    ## split-R-hat targets the same posterior.
    fits <- run_phase(
        phase_iterations = iterations,
        phase_z_burnin = z_probability_burnin,
        phase_use_data_priors = FALSE,
        phase_priors = priors,
        phase_epsilon = production_epsilon,
        phase_mcse_stop = mcse_stop,
        phase_relabel = relabel
    )
    if (isTRUE(use_data_priors)) {
        fits <- lapply(fits, function(fit) {
            fit$regression_prior$method <-
                "shared_soft_allocation_empirical_Bayes_pilot_frozen"
            fit$regression_prior$shared_across_chains <- TRUE
            fit$regression_prior$pilot_iterations <- pilot_summary$iterations
            fit$regression_prior$pilot_chains <- pilot_summary$chains
            fit
        })
    }
    all_fits <- fits
    mode_selection <- .hicpotts_mode_consensus(
        all_fits,
        agreement_threshold = mode_agreement_threshold,
        minimum_mode_chains = minimum_mode_chains
    )
    if (isTRUE(mode_screen) && isTRUE(mode_selection$consensus_available)) {
        fits <- all_fits[mode_selection$selected]
    } else {
        mode_selection$selected <- seq_along(all_fits)
    }
    mode_selection$excluded <- setdiff(
        seq_along(all_fits),
        mode_selection$selected
    )
    mode_selection$screening_enabled <- isTRUE(mode_screen)
    mode_selection$chain_table <- data.frame(
        chain = seq_along(all_fits), mode_group = mode_selection$component,
        selected_for_parameters = seq_along(
            all_fits
        ) %in% mode_selection$selected
    )
    diagnostics <- diagnose_hicpotts_fit(
        fits,
        burnin = burnin,
        minimum_component_cells = minimum_component_cells,
        minimum_ess = minimum_ess, maximum_rhat = maximum_rhat,
        gamma_boundary_tolerance = gamma_boundary_tolerance,
        maximum_gamma_boundary_fraction = maximum_gamma_boundary_fraction,
        minimum_gamma_unique = minimum_gamma_unique,
        covariate_diagnostics = covariates, relabel = FALSE
    )
    mode_pass <- !isTRUE(mode_screen) || isTRUE(
        mode_selection$consensus_available
    )
    mode_flag <- data.frame(
        criterion = "replicated_allocation_mode",
        passed = mode_pass,
        threshold = paste0(">=", minimum_mode_chains, " coherent chains"),
        stringsAsFactors = FALSE
    )
    ## M4: a second replicated mode fails WHOLE-POSTERIOR reliability even when
    ## the selected mode is internally healthy. `reliable` describes the full
    ## posterior; conditioning it on one mode while a rival mode is equally well
    ## replicated would overstate what the fit supports.
    unimodal_flag <- data.frame(
        criterion = "unimodal_posterior",
        passed = !isTRUE(mode_selection$multimodal),
        threshold = paste0(
            "no rival mode with >=", minimum_mode_chains,
            " chains"
        ),
        stringsAsFactors = FALSE
    )
    diagnostics$reliability_flags <- rbind(
        mode_flag, unimodal_flag,
        diagnostics$reliability_flags
    )
    if (isTRUE(mode_selection$multimodal)) {
        diagnostics$warnings <- unique(c(
            sprintf(
                paste0(
                    "Posterior is multimodal: %d rival allocation mode(s) ",
                    "each hold at ",
                    "least %d chains. Summaries below are CONDITIONAL on the ",
                    "selected ",
                    "mode and do not describe the full posterior."
                ),
                mode_selection$competing_modes, minimum_mode_chains
            ),
            diagnostics$warnings
        ))
    }
    ## Mode screening is already recorded as its own gates in
    ## reliability_flags, so there is nothing further to fold in here.
    if (length(mode_selection$excluded)) {
        diagnostics$warnings <- unique(c(
            sprintf(
                paste0("Chains %s were retained in all_fits but excluded from ",
                    "parameter pooling as incoherent modes."),
                paste(mode_selection$excluded, collapse = ", ")
            ),
            diagnostics$warnings
        ))
    }
    if (isTRUE(mode_screen) && !isTRUE(mode_selection$consensus_available)) {
        diagnostics$warnings <- unique(c(
            diagnostics$warnings,
            "No replicated allocation mode was found; no chains were excluded."
        ))
    }
    structure(
        list(
            fits = fits, all_fits = all_fits, mode_selection = mode_selection,
            diagnostics = diagnostics,
            covariate_diagnostics = covariates, priors = priors,
            empirical_bayes_pilot = pilot_summary,
            settings = list(
                N = N, dist = dist, gamma_start = gamma_prior,
                gamma_prior = gamma_prior,
                iterations = iterations, burnin = burnin, seeds = as.integer(
                    seeds
                ),
                z_probability_burnin = z_probability_burnin,
                initialization = initialization,
                mcse_stop = mcse_stop,
                mcse_min_iterations = mcse_min_iterations,
                mcse_check_interval = mcse_check_interval,
                mcse_relative_threshold = mcse_relative_threshold,
                default_iterations = 20000L,
                tempering_warmup = tempering_warmup,
                tempering_beta_min = tempering_beta_min,
                tempering_cycle = tempering_cycle,
                mode_screen = mode_screen,
                mode_agreement_threshold = mode_agreement_threshold,
                minimum_mode_chains = minimum_mode_chains,
                use_data_priors = use_data_priors,
                shared_empirical_bayes_prior = isTRUE(use_data_priors),
                empirical_bayes_pilot_iterations = if (is.null(
                    pilot_summary
                )) {
                    0L
                } else {
                    pilot_summary$iterations
                },
                production_abc_epsilon = production_epsilon,
                use_noise_relationship_prior = use_noise_relationship_prior,
                noise_link_sd = noise_link_sd,
                noise_order_strength = noise_order_strength,
                noise_order_width = noise_order_width,
                comp23_barrier_kappa = comp23_barrier_kappa,
                comp23_barrier_w = comp23_barrier_w,
                branch_swap_interval = branch_swap_interval,
                signal_block_move_interval = signal_block_move_interval,
                gamma_large_jump_probability = gamma_large_jump_probability,
                gamma_large_jump_multiplier = gamma_large_jump_multiplier,
                gamma_independence_probability = gamma_independence_probability,
                gamma_method = gamma_method,
                gamma_prior_shape1 = gamma_prior_shape1,
                gamma_prior_shape2 = gamma_prior_shape2,
                gamma_boundary_tolerance = gamma_boundary_tolerance,
                maximum_gamma_boundary_fraction =
                    maximum_gamma_boundary_fraction,
                minimum_gamma_unique = minimum_gamma_unique,
                minimum_component_cells = minimum_component_cells,
                minimum_ess = minimum_ess,
                maximum_rhat = maximum_rhat,
                slope_sd_standardized = if (use_data_priors) {
                    NA_real_
                } else {
                    slope_sd_standardized
                },
                relabel = relabel, mc_cores = as.integer(mc_cores),
                coefficient_scale = "original log1p-covariate scale"
            )
        ),
        class = "hicpotts_robust_fit"
    )
}

.hicpotts_matrix_summary <- function(y) {
    ly <- log1p(y)
    horizontal <- if (ncol(y) > 1L) {
        cbind(
            as.numeric(ly[, -ncol(y)]), as.numeric(ly[, -1L])
        )
    } else {
        NULL
    }
    vertical <- if (nrow(y) > 1L) {
        cbind(
            as.numeric(ly[-nrow(y), ]), as.numeric(ly[-1L, ])
        )
    } else {
        NULL
    }
    neighbour_pairs <- rbind(horizontal, vertical)
    neighbour_correlation <- if (is.null(neighbour_pairs) ||
        any(apply(neighbour_pairs, 2L, stats::sd) == 0)) {
        0
    } else {
        stats::cor(neighbour_pairs[, 1L], neighbour_pairs[, 2L])
    }
    diagonal_distance <- abs(row(y) - col(y))
    decay <- tapply(as.numeric(ly), as.numeric(diagonal_distance), mean)
    near_diagonal <- diagonal_distance <= 1L
    mean_y <- mean(y)
    list(
        zero_fraction = mean(y == 0),
        near_diagonal_zero_fraction = mean(y[near_diagonal] == 0),
        neighbour_correlation = unname(neighbour_correlation),
        diagonal_decay = decay,
        variance_to_mean = if (mean_y > 0) {
            stats::var(as.numeric(
                y
            )) / mean_y
        } else {
            NA_real_
        }
    )
}

#' Conditional posterior-predictive diagnostics for HiCPotts
#'
#' Generates replicated count matrices from posterior parameter draws,
#' conditional on each chain's final internal (unforced) latent allocation.
#' It reports cellwise log-count error, neighbour-correlation error,
#' diagonal-decay error, overall and near-diagonal zero-fraction error, and
#' variance-to-mean error. The reported forced-zero `z_final` is deliberately
#' not used for these parameter diagnostics.
#'
#' @param fit One fit, a list of independently seeded fits, or block-aware
#'   result from \code{combine_hicpotts_blocks()}.
#' @param x_vars Named HiCPotts covariates. May be omitted when a block-aware
#'   fit stored the structured output of \code{process_data()}.
#' @param y Observed count matrix. May be omitted when a block-aware fit stored
#'   the structured output of \code{process_data()}.
#' @param dist One of `Poisson`, `NB`, `ZIP` or `ZINB`.
#' @param burnin Initial draws to discard; defaults to half of each chain.
#' @param n_rep Number of replicated matrices.
#' @param seed Integer simulation seed.
#' @return A list with replicate discrepancies, interval summaries and observed
#'   matrix summaries.
#' @examples
#' N <- 6
#' y <- matrix(rpois(N * N, 4), N, N)
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' fits <- run_chain_betas(
#'     N = N, iterations = 5, x_vars = x_vars, y = y,
#'     dist = "ZINB", theta_start = 0.5, size_start = c(1, 1, 1),
#'     n_chains = 2, seeds = 1:2, robust = FALSE
#' )
#' pp <- posterior_predictive_hicpotts(
#'     fits, x_vars, y,
#'     dist = "ZINB", n_rep = 5
#' )
#' pp$summary
#' @export
posterior_predictive_hicpotts <- function(
    fit, x_vars = NULL, y = NULL, dist = "ZINB", burnin = NULL,
    n_rep = 100L, seed = 1L
) {
    dist <- match.arg(dist, c("Poisson", "NB", "ZIP", "ZINB"))
    if (inherits(fit, "hicpotts_block_fit")) {
        return(.hicpotts_posterior_predictive_blocks(
            fit, x_vars = x_vars, y = y, dist = dist, burnin = burnin,
            n_rep = n_rep, seed = seed
        ))
    }
    if (inherits(fit, "hicpotts_robust_fit")) fit <- fit$fits
    fits <- .hicpotts_fit_list(fit)
    fits <- lapply(fits, function(x) {
        if (is.null(x$relabel_rule)) {
            relabel_hicpotts(x)
        } else {
            x
        }
    })
    if (!is.matrix(y) || !is.numeric(y) || nrow(y) != ncol(y) ||
        any(!is.finite(y)) || any(y < 0)) {
        stop("'y' must be a finite non-negative square numeric matrix.")
    }
    n_rep <- as.integer(n_rep)
    if (length(n_rep) != 1L || is.na(n_rep) || n_rep < 1L) {
        stop("'n_rep' must be a positive integer.")
    }
    required <- c("distance", "GC", "TES", "ACC")
    if (!is.list(x_vars) || !all(required %in% names(x_vars))) {
        stop("'x_vars' must contain distance, GC, TES and ACC.")
    }
    select_covariate <- function(name) {
        value <- x_vars[[name]]
        if (is.list(value)) value <- value[[1L]]
        if (!is.matrix(value) || !identical(dim(value), dim(y)) ||
            any(!is.finite(value)) || any(value <= -1)) {
            stop(
                "Every covariate must be a finite matrix matching y and ",
                "greater than -1."
            )
        }
        log1p(as.numeric(value))
    }
    design <- cbind(1, vapply(required, select_covariate, numeric(length(y))))
    observed <- .hicpotts_matrix_summary(y)
    burnins <- vapply(fits, function(x) {
        n_iter <- nrow(x$chains[[1L]])
        value <- if (is.null(burnin)) {
            as.integer(floor(
                n_iter / 2L
            ))
        } else {
            as.integer(burnin)
        }
        if (length(value) != 1L || is.na(
            value
        ) || value < 0L || value >= n_iter) {
            stop(
                "'burnin' must be between zero and the number of draws minus ",
                "one."
            )
        }
        value
    }, integer(1))
    get_parameter_z <- function(x) {
        if (!is.null(x$z_parameter_final)) {
            return(x$z_parameter_final)
        }
        attr(x$z_final, "z_parameter_final", exact = TRUE)
    }
    if (any(vapply(fits, function(x) is.null(get_parameter_z(x)), logical(
        1
    )))) {
        stop(
            "The fit lacks z_parameter_final; refit with the updated sampler ",
            "for unforced diagnostics."
        )
    }
    withr::local_seed(as.integer(seed))
    discrepancies <- vector("list", n_rep)
    replicated <- vector("list", n_rep)
    for (r in seq_len(n_rep)) {
        chain_id <- sample.int(length(fits), 1L)
        current <- fits[[chain_id]]
        draw <- sample.int(nrow(current$chains[[1L]]) - burnins[chain_id], 1L) +
            burnins[chain_id]
        z <- get_parameter_z(current)
        if (!is.matrix(z) || !identical(dim(z), dim(y)) || !all(z %in% seq_len(
            3
        ))) {
            stop(
                "Each z_parameter_final must match y and contain labels 1, 2 ",
                "and 3."
            )
        }
        mu <- numeric(length(y))
        z_vector <- as.integer(z)
        for (component in seq_len(3)) {
            index <- which(z_vector == component)
            if (length(index)) {
                eta <- as.numeric(design[index, , drop = FALSE] %*%
                    current$chains[[component]][draw, seq_len(5)])
                mu[index] <- exp(pmin(eta, 20))
            }
        }
        if (dist %in% c("NB", "ZINB")) {
            size <- current$size[, draw]
            if (length(size) != 3L || any(!is.finite(size)) || any(size <= 0)) {
                stop(
                    "NB/ZINB posterior-predictive simulation requires ",
                    "positive size draws."
                )
            }
            y_rep <- stats::rnbinom(length(y), mu = mu, size = size[z_vector])
        } else {
            y_rep <- stats::rpois(length(y), lambda = mu)
        }
        if (dist %in% c("ZIP", "ZINB")) {
            theta <- current$theta[draw]
            if (!is.finite(theta) || theta < 0 || theta > 1) {
                stop(
                    "ZIP/ZINB posterior-predictive simulation requires theta ",
                    "in [0,1]."
                )
            }
            structural <- z_vector == 1L & stats::runif(length(y)) < theta
            y_rep[structural] <- 0
        }
        y_rep <- matrix(y_rep, nrow(y), ncol(y))
        rep_summary <- .hicpotts_matrix_summary(y_rep)
        common_decay <- intersect(
            names(observed$diagonal_decay),
            names(rep_summary$diagonal_decay)
        )
        discrepancies[[r]] <- data.frame(
            replicate = r, chain = chain_id, iteration = draw,
            S1_log_count = mean(abs(log1p(y) - log1p(y_rep))),
            S2_neighbour = abs(observed$neighbour_correlation -
                rep_summary$neighbour_correlation),
            diagonal_decay = sqrt(mean((observed$diagonal_decay[common_decay] -
                rep_summary$diagonal_decay[common_decay])^2)),
            zero_fraction = abs(
                observed$zero_fraction - rep_summary$zero_fraction
            ),
            near_diagonal_zero_fraction = abs(
                observed$near_diagonal_zero_fraction -
                    rep_summary$near_diagonal_zero_fraction
            ),
            variance_to_mean = abs(observed$variance_to_mean -
                rep_summary$variance_to_mean)
        )
        replicated[[r]] <- unlist(rep_summary[c(
            "zero_fraction",
            "near_diagonal_zero_fraction", "neighbour_correlation",
            "variance_to_mean"
        )])
    }
    discrepancy_table <- do.call(rbind, discrepancies)
    metric_names <- setdiff(
        names(discrepancy_table),
        c("replicate", "chain", "iteration")
    )
    summary <- do.call(rbind, lapply(metric_names, function(metric) {
        values <- discrepancy_table[[metric]]
        q <- stats::quantile(values, c(0.025, 0.5, 0.975),
            na.rm = TRUE,
            names = FALSE
        )
        data.frame(metric = metric, median = q[2L], lower = q[1L], upper = q[
            3L
        ])
    }))
    rownames(summary) <- NULL
    list(
        discrepancies = discrepancy_table, summary = summary,
        observed = observed, replicated = do.call(rbind, replicated),
        conditioning = "final internal unforced Z"
    )
}

#' Compare posterior-predictive NB and ZINB performance
#'
#' @param family_fits Named output from `run_hicpotts_family_sensitivity()`.
#' @inheritParams posterior_predictive_hicpotts
#' @return A list containing family-specific diagnostics and a combined table;
#'   smaller discrepancy values indicate better reproduction of the data.
#' @examples
#' N <- 6
#' y <- matrix(rpois(N * N, 4), N, N)
#' mk <- function() list(matrix(runif(N * N), N, N))
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' fam <- run_hicpotts_family_sensitivity(
#'     N = N, iterations = 5, x_vars = x_vars,
#'     y = y, size_start = c(1, 1, 1), seeds = 1:2
#' )
#' comparison <- compare_hicpotts_families(fam, x_vars, y, n_rep = 5)
#' comparison$comparison
#' @noRd
compare_hicpotts_families <- function(
    family_fits, x_vars, y, burnin = NULL, n_rep = 100L, seed = 1L
) {
    if (!is.list(family_fits) ||
        !all(c("NB", "ZINB") %in% names(family_fits))) {
        stop("'family_fits' must contain named NB and ZINB fits.")
    }
    results <- lapply(c("NB", "ZINB"), function(family) {
        posterior_predictive_hicpotts(family_fits[[family]], x_vars, y,
            dist = family, burnin = burnin, n_rep = n_rep,
            seed = as.integer(seed) + match(family, c("NB", "ZINB")) - 1L
        )
    })
    names(results) <- c("NB", "ZINB")
    comparison <- do.call(rbind, lapply(names(results), function(family) {
        data.frame(
            family = family, results[[family]]$summary,
            stringsAsFactors = FALSE
        )
    }))
    rownames(comparison) <- NULL
    list(
        comparison = comparison, diagnostics = results,
        interpretation = paste0("Smaller posterior-predictive discrepancy is ",
            "better; prefer ZINB only when zero diagnostics improve ",
            "materially and theta is identified.")
    )
}

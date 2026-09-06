#' Fit one or more Hi-C datasets with HiCPotts
#'
#' `run_chain_betas()` is the single public fitting entry point. It can run a
#' user-selected number of independently seeded chains for one matrix or for
#' each matrix in a list. With `robust = TRUE`, it also constructs priors,
#' checks covariates, screens replicated allocation modes and applies the
#' parameter-reliability gates. With `robust = FALSE`, it returns the raw
#' relabelled chains. These choices change orchestration and reporting, not the
#' fitted likelihood or sampler target.
#'
#' @param N Number of rows and columns in every square interaction matrix.
#' @param x_vars Named list containing `distance`, `GC`, `TES` and `ACC`.
#' Entries may be matrices for one dataset or lists of matrices matching `y`.
#' @param y One numeric count matrix, or a list of independently fitted count
#' matrices.
#' @param dist One of `"Poisson"`, `"NB"`, `"ZIP"` or `"ZINB"`.
#' @param gamma_start Initial value of the Potts spatial parameter. This is a
#' starting value, not the parameters of its Beta prior.
#' @param iterations Positive maximum MCMC updates per chain.
#' @param burnin Initial iterations excluded by robust diagnostics. `NULL` uses
#' half the requested iterations, capped at 5,000.
#' @param n_chains Positive number of independently seeded chains per dataset.
#' The default of one preserves the historical `run_chain_betas()` result; four
#' or more chains are recommended when convergence will be assessed.
#' @param seeds Optional integer vector with exactly `n_chains` values. When
#' omitted, `seq_len(n_chains)` is used; when `n_chains` is omitted, the length
#' of `seeds` determines it.
#' @param robust Logical. If `TRUE`, return the diagnostic and mode-screened
#' robust workflow; if `FALSE`, return raw fitted chains.
#' @param initialization Starting allocation strategy. `"auto"` uses two
#' likelihood-informed starts followed by distance-adjusted and
#' noise-anchored-random starts, recycled over the requested chains. A single
#' method is used for every chain; alternatively supply exactly one method per
#' chain. Available methods are `"likelihood_informed"`, `"count_quantile"`,
#' `"distance_adjusted"`, `"noise_anchored_random"` and `"random"`.
#' @param theta_start Optional starting zero-inflation probability.
#' @param size_start Optional three-element vector of starting NB2 dispersions.
#' @param use_data_priors Whether to estimate empirical-Bayes regression-prior
#' parameters. With `robust = TRUE`, a multi-chain pilot constructs one prior
#' that is shared and frozen across all production chains; without the robust
#' workflow each chain performs its historical warm-up update.
#' @param user_fixed_priors Optional fixed regression priors used when
#' `use_data_priors = FALSE`.
#' @param gamma_prior_shape Positive two-element numeric vector containing the
#' shape parameters of the Beta prior on gamma. The default `c(1, 1)` is
#' uniform on `(0, 1)`.
#' @param slope_sd_standardized Standardized slope prior SD used when fixed
#' priors are constructed automatically by the robust workflow.
#' @param epsilon Optional fixed positive ABC bandwidth. `NULL` calibrates one
#' shared bandwidth for all chains of a dataset.
#' @param distance_metric Retained for source compatibility and ignored.
#' @param relabel Whether to apply the model-based component relabelling rule.
#' @param mc_cores Positive number of chains to execute concurrently.
#' @param abc_potts_sweeps Non-negative number of Gibbs sweeps used to
#' equilibrate each simulated auxiliary Potts field. Zero selects the
#' lattice-dependent default, `max(100, 100 * N^2 / 400)`, which keeps the
#' auxiliary fields equally well equilibrated as `N` grows. That default
#' makes the ABC step scale as `N^4` and dominates run time on large
#' lattices: at `N = 100` it is 2500 sweeps. A smaller value is much
#' faster but under-equilibrates the auxiliary fields, which biases the
#' estimated spatial coupling upward.
#' @param abc_sim_reps Positive number of auxiliary fields averaged per ABC
#' proposal. Cost is linear in this value; fewer replicates give a noisier
#' acceptance decision for the spatial coupling.
#' @param verbose Logical. Print sampler progress when `TRUE`.
#' @param progress_interval Positive number of iterations between progress
#' messages when `verbose = TRUE`.
#' @param abc_epsilon_quantile Quantile in `(0, 1)` used to calibrate the
#' Gaussian ABC bandwidth when `epsilon = NULL`.
#' @param gamma_update_interval Positive number of iterations between ABC gamma
#' updates. Larger values reduce computation but produce fewer gamma moves.
#' @param z_probability_burnin Iteration after which latent-state membership
#' frequencies are accumulated. `NULL` uses `burnin`.
#' @param mcse_stop Whether the sampler may stop early after satisfying the
#' relative Monte Carlo standard-error rule for every monitored parameter.
#' @param mcse_min_iterations Minimum completed iterations before MCSE stopping
#' is considered.
#' @param mcse_check_interval Positive number of iterations between MCSE checks.
#' @param mcse_relative_threshold Maximum batch-means MCSE divided by posterior
#' standard deviation required for early stopping.
#' @param diagnostic_control Named list overriding robust diagnostic defaults.
#' Recognised entries are `minimum_component_cells`, `minimum_ess`,
#' `maximum_rhat`, `gamma_boundary_tolerance`,
#' `maximum_gamma_boundary_fraction`, `minimum_gamma_unique`,
#' `mode_agreement_threshold` and `minimum_mode_chains`. These settings affect
#' reliability reporting and mode selection, not the sampler target.
#' @param gamma_prior Deprecated compatibility alias for `gamma_start`. New
#' code should use `gamma_start`.
#'
#' @details
#' Initialisation, progress reporting, MCSE stopping, gamma-update frequency,
#' ABC accuracy/workload and diagnostic thresholds are public controls. The
#' warm-up tempering schedule, component-defining relationship prior,
#' component-2/3 barrier and reversible allocation-move settings remain fixed
#' to preserve the package's component definitions and default sampler.
#'
#' `abc_potts_sweeps`, `abc_sim_reps` and `gamma_update_interval` are the main
#' ABC cost controls. Their defaults reproduce the package workflow; reducing
#' sweeps or replicates trades auxiliary-simulation fidelity for speed, while
#' increasing the update interval reduces the number of gamma transitions.
#'
#' @return For one matrix, a `hicpotts_robust_fit` when `robust = TRUE` or a
#'   list of fitted chains otherwise. For a list of matrices, a named list with
#'   one such result per dataset. For backward compatibility, a non-robust fit
#'   of a list with `n_chains = 1` is returned as one flat fit per dataset.
#' @examples
#' N <- 5
#' y <- matrix(rpois(N * N, 4), N, N)
#' mk <- function() matrix(runif(N * N), N, N)
#' x_vars <- list(distance = mk(), GC = mk(), TES = mk(), ACC = mk())
#' fits <- run_chain_betas(
#'     N = N, x_vars = x_vars, y = y, dist = "Poisson",
#'     iterations = 5, n_chains = 2, robust = FALSE
#' )
#' length(fits)
#' @export
run_chain_betas <- function(
    N, x_vars, y, dist = "ZIP", gamma_start = 0.3,
    iterations = 20000L, burnin = NULL, n_chains = 1L, seeds = NULL,
    robust = FALSE,
    initialization = "auto",
    theta_start = NULL, size_start = NULL,
    use_data_priors = TRUE, user_fixed_priors = NULL,
    gamma_prior_shape = c(1, 1),
    epsilon = NULL,
    abc_potts_sweeps = 0L, abc_sim_reps = 4L,
    distance_metric = "manhattan", relabel = TRUE, mc_cores = 1L,
    verbose = FALSE, progress_interval = 50L,
    slope_sd_standardized = 0.5,
    abc_epsilon_quantile = 0.10, gamma_update_interval = 5L,
    z_probability_burnin = NULL,
    mcse_stop = TRUE, mcse_min_iterations = 10000L,
    mcse_check_interval = 500L, mcse_relative_threshold = 0.05,
    diagnostic_control = list(), gamma_prior = NULL
) {
    gamma_start_was_missing <- missing(gamma_start)
    if (!is.null(gamma_prior)) {
        if (!is.numeric(gamma_prior) || length(gamma_prior) != 1L ||
            !is.finite(gamma_prior) || gamma_prior <= 0 || gamma_prior >= 1) {
            stop("'gamma_prior' must be strictly between 0 and 1.",
                call. = FALSE
            )
        }
        if (!gamma_start_was_missing &&
            !isTRUE(all.equal(as.numeric(gamma_start),
                as.numeric(gamma_prior)))) {
            stop(
                "Supply only 'gamma_start'; the deprecated 'gamma_prior' ",
                "alias was also supplied with a different value.",
                call. = FALSE
            )
        }
        gamma_start <- as.numeric(gamma_prior)
        warning(
            "'gamma_prior' is deprecated; use 'gamma_start' for the initial ",
            "Potts parameter value.",
            call. = FALSE
        )
    }
    if (!is.numeric(gamma_start) || length(gamma_start) != 1L ||
        !is.finite(gamma_start) || gamma_start <= 0 || gamma_start >= 1) {
        stop("'gamma_start' must be strictly between 0 and 1.",
            call. = FALSE
        )
    }
    gamma_start <- as.numeric(gamma_start)

    if (!is.numeric(gamma_prior_shape) || length(gamma_prior_shape) != 2L ||
        any(!is.finite(gamma_prior_shape)) || any(gamma_prior_shape <= 0)) {
        stop("'gamma_prior_shape' must contain two positive finite values.",
            call. = FALSE
        )
    }
    gamma_prior_shape <- as.numeric(gamma_prior_shape)
    gamma_prior_shape1 <- gamma_prior_shape[[1L]]
    gamma_prior_shape2 <- gamma_prior_shape[[2L]]

    diagnostic_defaults <- list(
        minimum_component_cells = 0L,
        minimum_ess = 200,
        maximum_rhat = 1.01,
        gamma_boundary_tolerance = 0.01,
        maximum_gamma_boundary_fraction = 0.95,
        minimum_gamma_unique = 20L,
        mode_agreement_threshold = 0.8,
        minimum_mode_chains = 2L
    )
    if (!is.list(diagnostic_control) ||
        (length(diagnostic_control) && is.null(names(diagnostic_control)))) {
        stop("'diagnostic_control' must be a named list.", call. = FALSE)
    }
    if (length(diagnostic_control) &&
        (any(!nzchar(names(diagnostic_control))) ||
            anyDuplicated(names(diagnostic_control)))) {
        stop(
            "'diagnostic_control' must have unique, non-empty names.",
            call. = FALSE
        )
    }
    unknown_diagnostics <- setdiff(
        names(diagnostic_control), names(diagnostic_defaults)
    )
    if (length(unknown_diagnostics)) {
        stop(
            "Unknown 'diagnostic_control' setting(s): ",
            paste(unknown_diagnostics, collapse = ", "), ".",
            call. = FALSE
        )
    }
    diagnostics <- utils::modifyList(
        diagnostic_defaults, diagnostic_control, keep.null = FALSE
    )
    positive_integer <- function(x, name, minimum = 1L) {
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
            x < minimum || x != as.integer(x)) {
            stop(sprintf("'%s' must be one integer of at least %d.",
                name, minimum
            ), call. = FALSE)
        }
        as.integer(x)
    }
    positive_number <- function(x, name) {
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
            stop(sprintf("'%s' must be one positive finite number.", name),
                call. = FALSE
            )
        }
        as.numeric(x)
    }
    minimum_component_cells <- positive_integer(
        diagnostics$minimum_component_cells, "minimum_component_cells",
        minimum = 0L
    )
    minimum_ess <- positive_number(diagnostics$minimum_ess, "minimum_ess")
    maximum_rhat <- positive_number(diagnostics$maximum_rhat, "maximum_rhat")
    if (maximum_rhat <= 1) {
        stop("'maximum_rhat' must be greater than one.", call. = FALSE)
    }
    gamma_boundary_tolerance <- positive_number(
        diagnostics$gamma_boundary_tolerance, "gamma_boundary_tolerance"
    )
    if (gamma_boundary_tolerance >= 0.5) {
        stop("'gamma_boundary_tolerance' must be below 0.5.", call. = FALSE)
    }
    maximum_gamma_boundary_fraction <- positive_number(
        diagnostics$maximum_gamma_boundary_fraction,
        "maximum_gamma_boundary_fraction"
    )
    if (maximum_gamma_boundary_fraction > 1) {
        stop("'maximum_gamma_boundary_fraction' must not exceed one.",
            call. = FALSE
        )
    }
    minimum_gamma_unique <- positive_integer(
        diagnostics$minimum_gamma_unique, "minimum_gamma_unique", 2L
    )
    mode_agreement_threshold <- positive_number(
        diagnostics$mode_agreement_threshold, "mode_agreement_threshold"
    )
    if (mode_agreement_threshold > 1) {
        stop("'mode_agreement_threshold' must not exceed one.",
            call. = FALSE
        )
    }
    minimum_mode_chains <- positive_integer(
        diagnostics$minimum_mode_chains, "minimum_mode_chains"
    )

    ## ---- fixed component-definition and low-level proposal settings ------
    mode_screen <- TRUE
    tempering_warmup <- 0L
    tempering_beta_min <- 0.30
    tempering_cycle <- 500L
    comp23_barrier_kappa <- 10
    comp23_barrier_w <- 0.3
    use_noise_relationship_prior <- TRUE
    noise_link_sd <- 0.5
    noise_order_strength <- 10
    noise_order_width <- 0.5
    branch_swap_interval <- 1L
    signal_block_move_interval <- 5L
    gamma_large_jump_probability <- 0.10
    gamma_large_jump_multiplier <- 4
    gamma_independence_probability <- 0.05
    gamma_method <- "abc"

    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
    }
    progress_interval <- positive_integer(
        progress_interval, "progress_interval"
    )
    if (!is.logical(mcse_stop) || length(mcse_stop) != 1L ||
        is.na(mcse_stop)) {
        stop("'mcse_stop' must be TRUE or FALSE.", call. = FALSE)
    }
    mcse_min_iterations <- positive_integer(
        mcse_min_iterations, "mcse_min_iterations"
    )
    mcse_check_interval <- positive_integer(
        mcse_check_interval, "mcse_check_interval"
    )
    mcse_relative_threshold <- positive_number(
        mcse_relative_threshold, "mcse_relative_threshold"
    )
    if (!is.numeric(abc_epsilon_quantile) ||
        length(abc_epsilon_quantile) != 1L ||
        !is.finite(abc_epsilon_quantile) || abc_epsilon_quantile <= 0 ||
        abc_epsilon_quantile >= 1) {
        stop("'abc_epsilon_quantile' must be strictly between 0 and 1.",
            call. = FALSE
        )
    }
    abc_epsilon_quantile <- as.numeric(abc_epsilon_quantile)
    gamma_update_interval <- positive_integer(
        gamma_update_interval, "gamma_update_interval"
    )

    if (!is.numeric(abc_potts_sweeps) || length(abc_potts_sweeps) != 1L ||
        !is.finite(abc_potts_sweeps) || abc_potts_sweeps < 0 ||
        abc_potts_sweeps != as.integer(abc_potts_sweeps)) {
        stop(
            "'abc_potts_sweeps' must be one non-negative integer; zero ",
            "selects the lattice-dependent default.",
            call. = FALSE
        )
    }
    abc_potts_sweeps <- as.integer(abc_potts_sweeps)
    if (!is.numeric(abc_sim_reps) || length(abc_sim_reps) != 1L ||
        !is.finite(abc_sim_reps) || abc_sim_reps < 1 ||
        abc_sim_reps != as.integer(abc_sim_reps)) {
        stop("'abc_sim_reps' must be one positive integer.", call. = FALSE)
    }
    abc_sim_reps <- as.integer(abc_sim_reps)
    if (!is.numeric(iterations) || length(iterations) != 1L ||
        !is.finite(iterations) || iterations < 1 ||
        iterations != as.integer(iterations)) {
        stop("'iterations' must be one positive integer.", call. = FALSE)
    }
    iterations <- as.integer(iterations)
    if (missing(n_chains) && !is.null(seeds)) n_chains <- length(seeds)
    if (!is.numeric(n_chains) || length(n_chains) != 1L ||
        !is.finite(n_chains) || n_chains < 1 ||
        n_chains != as.integer(n_chains)) {
        stop("'n_chains' must be one positive integer.", call. = FALSE)
    }
    n_chains <- as.integer(n_chains)
    if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
        stop("'robust' must be TRUE or FALSE.", call. = FALSE)
    }
    if (is.null(seeds)) {
        seeds <- seq_len(n_chains)
    } else {
        seeds <- as.integer(seeds)
        if (length(seeds) != n_chains || anyNA(seeds)) {
            stop("'seeds' must contain exactly 'n_chains' integers.",
                call. = FALSE
            )
        }
    }
    allowed_initializations <- c(
        "likelihood_informed", "count_quantile", "distance_adjusted",
        "noise_anchored_random", "random"
    )
    automatic_initialization <- c(
        "likelihood_informed", "likelihood_informed",
        "distance_adjusted", "noise_anchored_random"
    )
    if (!is.character(initialization) || !length(initialization) ||
        anyNA(initialization) || any(!nzchar(initialization))) {
        stop("'initialization' must contain one or more method names.",
            call. = FALSE
        )
    }
    if ("auto" %in% initialization) {
        if (length(initialization) != 1L) {
            stop("'auto' cannot be combined with explicit initializations.",
                call. = FALSE
            )
        }
        initialization <- rep(
            automatic_initialization, length.out = n_chains
        )
    } else {
        unknown_initializations <- setdiff(
            initialization, allowed_initializations
        )
        if (length(unknown_initializations)) {
            stop(
                "Unknown initialization method(s): ",
                paste(unknown_initializations, collapse = ", "), ".",
                call. = FALSE
            )
        }
        if (!length(initialization) %in% c(1L, n_chains)) {
            stop(
                "'initialization' must contain one method or exactly ",
                "'n_chains' methods.",
                call. = FALSE
            )
        }
        initialization <- rep(initialization, length.out = n_chains)
    }
    if (is.null(burnin)) burnin <- min(5000L, floor(iterations / 2L))
    if (!is.numeric(burnin) || length(burnin) != 1L || !is.finite(burnin) ||
        burnin != as.integer(burnin)) {
        stop("'burnin' must be one integer.", call. = FALSE)
    }
    burnin <- as.integer(burnin)
    if (length(burnin) != 1L || is.na(burnin) || burnin < 0L ||
        burnin >= iterations + 1L) {
        stop("'burnin' must be between zero and iterations.", call. = FALSE)
    }
    if (is.null(z_probability_burnin)) z_probability_burnin <- burnin
    if (!is.numeric(z_probability_burnin) ||
        length(z_probability_burnin) != 1L ||
        !is.finite(z_probability_burnin) ||
        z_probability_burnin != as.integer(z_probability_burnin) ||
        z_probability_burnin < 0L || z_probability_burnin >= iterations) {
        stop(
            "'z_probability_burnin' must be an integer between zero and ",
            "iterations - 1.",
            call. = FALSE
        )
    }
    z_probability_burnin <- as.integer(z_probability_burnin)

    y_was_list <- is.list(y) && !is.matrix(y)
    datasets <- if (y_was_list) y else list(y)
    if (!length(datasets)) {
        stop("'y' must contain at least one matrix.",
            call. = FALSE
        )
    }
    required <- c("distance", "GC", "TES", "ACC")
    if (!is.list(x_vars) || !all(required %in% names(x_vars))) {
        stop(
            "'x_vars' must be a named list containing distance, GC, TES and ",
            "ACC.",
            call. = FALSE
        )
    }

    dataset_covariates <- function(i) {
        out <- lapply(required, function(nm) {
            entry <- x_vars[[nm]]
            if (is.list(entry) && !is.matrix(entry)) {
                if (length(entry) < i) {
                    stop(sprintf("x_vars$%s has fewer matrices than 'y'.", nm),
                        call. = FALSE
                    )
                }
                entry[[i]]
            } else {
                if (length(datasets) > 1L) {
                    template <- paste0(
                        "x_vars$%s must be a list with one matrix ",
                        "per dataset."
                    )
                    stop(sprintf(template, nm), call. = FALSE)
                }
                entry
            }
        })
        ## The native sampler consumes the historical one-matrix-per-entry list
        ## layout even when the public interface accepts plain matrices.
        stats::setNames(lapply(out, list), required)
    }

    fit_one <- function(i) {
        yi <- datasets[[i]]
        xi <- dataset_covariates(i)
        common <- list(
            N = N, x_vars = xi, y = yi, dist = dist,
            gamma_prior = gamma_start, iterations = iterations,
            theta_start = theta_start, size_start = size_start,
            epsilon = epsilon, distance_metric = distance_metric,
            seeds = seeds, initialization = initialization,
            use_data_priors = use_data_priors,
            user_fixed_priors = user_fixed_priors,
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
            abc_potts_sweeps = abc_potts_sweeps,
            abc_sim_reps = abc_sim_reps,
            z_probability_burnin = z_probability_burnin,
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
            gamma_method = gamma_method, relabel = relabel,
            mc_cores = mc_cores, verbose = verbose,
            progress_interval = progress_interval
        )
        if (isTRUE(robust)) {
            return(do.call(fit_hicpotts_robust, c(common, list(
                burnin = burnin,
                slope_sd_standardized = slope_sd_standardized,
                minimum_component_cells = minimum_component_cells,
                minimum_ess = minimum_ess, maximum_rhat = maximum_rhat,
                gamma_boundary_tolerance = gamma_boundary_tolerance,
                maximum_gamma_boundary_fraction =
                    maximum_gamma_boundary_fraction,
                minimum_gamma_unique = minimum_gamma_unique,
                mode_screen = mode_screen,
                mode_agreement_threshold = mode_agreement_threshold,
                minimum_mode_chains = minimum_mode_chains
            ))))
        }
        do.call(run_hicpotts_chains, common)
    }

    result <- lapply(seq_along(datasets), fit_one)
    names(result) <- names(datasets)
    if (length(result) == 1L && !y_was_list) {
        return(result[[1L]])
    }
    if (!isTRUE(robust) && n_chains == 1L) {
        result <- lapply(result, `[[`, 1L)
    }
    if (is.null(names(result)) || any(!nzchar(names(result)))) {
        names(result) <- paste0("dataset", seq_along(result))
    }
    result
}

#' Structured sensitivity analysis for a HiCPotts fit
#'
#' @description
#' Refits the model under systematic perturbations of the choices that are not
#' determined by the data, and reports how far the parameters, the membership
#' probabilities and the three-way classification move as a result.
#'
#' Several HiCPotts settings are analyst choices rather than data: the priors,
#' the ABC tolerance and simulation effort, the component-1/3 relationship
#' prior, and whether chains are pooled. A result that changes materially when
#' these are varied is conditional on them, and that dependence should be
#' measured and reported rather than assumed away.
#'
#' Each scenario is a full refit under one perturbation, compared against a
#' baseline fit under the current defaults. Classification agreement is reported
#' as the fraction of contacts keeping the same three-way label, and probability
#' movement as the mean per-cell total-variation distance, with the the
#' biological relabelling applied first. The comparison does not subsequently
#' swap components 2/3 to improve agreement, because doing so would erase a real
#' disagreement about signal versus false signal.
#'
#' @param N,x_vars,y,dist,iterations Arguments forwarded to
#'   \code{run_metropolis_MCMC_betas()} to build each fit.
#' @param scenarios Which families of perturbation to run. Any subset of
#'   \code{"priors"}, \code{"abc"}, \code{"relationship"},
#'   \code{"pooling"}, or the backwards-compatibility check \code{"barrier"}.
#' @param seed Base RNG seed; each scenario is run from the same seed so that
#'   differences reflect the setting, not the randomness.
#' @param min_draws Passed to \code{classify_hicpotts()}.
#' @param ... Further arguments forwarded to \code{run_metropolis_MCMC_betas()}.
#'
#' @return An object of class \code{hicpotts_sensitivity}: a list with a
#'   \code{summary} data frame (one row per scenario) and the per-scenario
#'   parameter tables.
#'
#' @seealso \code{\link{classify_hicpotts}},
#' \code{\link{validate_hicpotts_simulation}}
#' @noRd
hicpotts_sensitivity <- function(
    N, x_vars, y, dist = "ZIP", iterations = 2000L,
    scenarios = c("priors", "abc", "relationship", "pooling"),
    seed = 1L, min_draws = 50L, ...
) {
    scenarios <- match.arg(
        scenarios,
        c("priors", "abc", "relationship", "pooling", "barrier"),
        several.ok = TRUE
    )
    extra <- list(...)

    run_one <- function(overrides, chain_seed = seed) {
        args <- c(list(
            N = N, gamma_prior = 0.4, iterations = iterations,
            x_vars = x_vars, y = y, use_data_priors = TRUE,
            dist = dist
        ), extra, overrides)
        args <- args[!duplicated(names(args), fromLast = TRUE)]
        withr::local_seed(as.integer(chain_seed))
        do.call(run_metropolis_MCMC_betas, args)
    }

    ## Baseline under current defaults.
    baseline_fit <- run_one(list())
    baseline_class <- classify_hicpotts(baseline_fit,
        relabel = TRUE,
        min_draws = min_draws
    )
    baseline_probs <- as.matrix(baseline_class[, c("prob1", "prob2", "prob3")])

    ## Both fits have already been biologically relabelled. A second, agreement-
    ## maximising 2/3 swap would conceal scientific instability.
    compare <- function(cls) {
        p <- as.matrix(cls[, c("prob1", "prob2", "prob3")])
        tv <- 0.5 * rowSums(abs(p - baseline_probs))
        lab <- as.integer(cls$map_component)
        list(
            mean_tv = mean(tv), max_tv = max(tv),
            label_agreement = mean(lab == as.integer(
                baseline_class$map_component
            )),
            permutation_aligned = FALSE
        )
    }

    definitions <- list()
    if ("priors" %in% scenarios) {
        definitions <- c(definitions, list(
            list(
                family = "priors", label = "gamma prior Beta(2,2)",
                overrides = list(gamma_prior_shape1 = 2, gamma_prior_shape2 = 2)
            ),
            list(
                family = "priors", label = "data-driven priors off",
                overrides = list(use_data_priors = FALSE)
            )
        ))
    }
    if ("abc" %in% scenarios) {
        definitions <- c(definitions, list(
            list(
                family = "abc", label = "tighter ABC tolerance (q=0.05)",
                overrides = list(
                    gamma_method = "abc",
                    abc_epsilon_quantile = 0.05
                )
            ),
            list(
                family = "abc", label = "looser ABC tolerance (q=0.30)",
                overrides = list(
                    gamma_method = "abc",
                    abc_epsilon_quantile = 0.30
                )
            ),
            list(
                family = "abc", label = "more ABC simulations (reps=8)",
                overrides = list(gamma_method = "abc", abc_sim_reps = 8L)
            )
        ))
    }
    if ("relationship" %in% scenarios) {
        definitions <- c(definitions, list(
            list(
                family = "relationship", label = "relationship prior disabled",
                overrides = list(use_noise_relationship_prior = FALSE)
            ),
            list(
                family = "relationship", label = "tighter slope link (sd=0.25)",
                overrides = list(noise_link_sd = 0.25)
            ),
            list(
                family = "relationship", label = "looser slope link (sd=1.0)",
                overrides = list(noise_link_sd = 1.0)
            ),
            list(
                family = "relationship", label = "weaker elevation preference",
                overrides = list(noise_order_strength = 2)
            )
        ))
    }
    if ("barrier" %in% scenarios) {
        definitions <- c(definitions, list(
            list(
                family = "barrier", label = "legacy barrier enabled (kappa=2)",
                overrides = list(comp23_barrier_kappa = 2)
            ),
            list(
                family = "barrier", label = "legacy strong barrier (kappa=25)",
                overrides = list(comp23_barrier_kappa = 25)
            ),
            list(
                family = "barrier", label = "legacy wide barrier",
                overrides = list(
                    comp23_barrier_kappa = 10,
                    comp23_barrier_w = 0.6
                )
            )
        ))
    }

    rows <- list()
    parameter_tables <- list()
    for (d in definitions) {
        fit <- run_one(d$overrides)
        cls <- classify_hicpotts(fit, relabel = TRUE, min_draws = min_draws)
        cmp <- compare(cls)
        keep <- (floor(iterations / 2) + 1L):(iterations + 1L)
        parameter_tables[[d$label]] <- data.frame(
            parameter = c(
                "intercept1", "intercept2", "intercept3", "gamma",
                "theta"
            ),
            estimate = c(
                mean(fit$chains[[1]][keep, 1]), mean(fit$chains[[2]][keep, 1]),
                mean(fit$chains[[3]][keep, 1]), mean(fit$gamma[keep]),
                mean(fit$theta[keep])
            ), stringsAsFactors = FALSE
        )
        rows[[length(rows) + 1L]] <- data.frame(
            family = d$family, scenario = d$label,
            label_agreement = cmp$label_agreement,
            mean_probability_shift = cmp$mean_tv,
            max_probability_shift = cmp$max_tv,
            permutation_aligned = cmp$permutation_aligned,
            relabel_basis = if (is.null(fit$relabel_basis)) {
                NA_character_
            } else {
                fit$relabel_basis
            },
            stringsAsFactors = FALSE
        )
    }

    ## Chain pooling: does the answer depend on pooling replicate chains?
    if ("pooling" %in% scenarios) {
        second <- run_one(list(), chain_seed = seed + 1000L)
        pooled_class <- classify_hicpotts(list(baseline_fit, second),
            relabel = TRUE,
            min_draws = min_draws
        )
        cmp <- compare(pooled_class)
        rows[[length(rows) + 1L]] <- data.frame(
            family = "pooling", scenario = "two chains pooled vs single chain",
            label_agreement = cmp$label_agreement,
            mean_probability_shift = cmp$mean_tv,
            max_probability_shift = cmp$max_tv,
            permutation_aligned = cmp$permutation_aligned,
            relabel_basis = NA_character_, stringsAsFactors = FALSE
        )
    }

    summary <- do.call(rbind, rows)
    rownames(summary) <- NULL
    structure(
        list(
            summary = summary,
            parameters = parameter_tables,
            baseline_classification = baseline_class,
            iterations = iterations, seed = seed
        ),
        class = "hicpotts_sensitivity"
    )
}

#' @export
print.hicpotts_sensitivity <- function(x, ...) {
    cat("HiCPotts sensitivity analysis\n")
    cat(sprintf(
        "  iterations: %d   base seed: %d   scenarios: %d\n",
        x$iterations, x$seed, nrow(x$summary)
    ))
    s <- x$summary
    cat(sprintf(
        "  worst label agreement: %.3f (%s)\n",
        min(s$label_agreement), s$scenario[which.min(s$label_agreement)]
    ))
    cat(sprintf(
        "  largest mean probability shift: %.4f (%s)\n",
        max(s$mean_probability_shift),
        s$scenario[which.max(s$mean_probability_shift)]
    ))
    print(s[, c(
        "family", "scenario", "label_agreement",
        "mean_probability_shift"
    )], row.names = FALSE)
    invisible(x)
}

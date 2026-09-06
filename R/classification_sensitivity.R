#' Compare official and parameter-based HiCPotts classifications
#'
#' Treats the post-burn-in latent-state frequencies returned by
#' \code{classify_hicpotts()} as the official three-component classification,
#' then recomputes a secondary allocation from retained parameter draws using
#' \code{compute_HMRFHiC_probabilities()} with the Potts mean-field term.  The
#' comparison is a sensitivity analysis only: it never replaces, averages, or
#' otherwise changes the official noise/signal/false-signal labels.
#'
#' @param fit A HiCPotts fit, a list of fits, or a
#'   \code{hicpotts_robust_fit}.
#' @param data Data frame in the same column-major lattice order used for
#'   fitting, with the columns required by
#'   \code{compute_HMRFHiC_probabilities()}.
#' @param N Positive lattice dimension satisfying \code{nrow(data) == N^2}.
#' @param dist Count family used to fit the model.
#' @param use For a robust fit, compare the selected allocation mode or all
#'   retained chains. Passed to \code{classify_hicpotts()}.
#' @param relabel Whether to apply the package's component-2/3 semantic
#'   relabelling to both calculations.
#' @param min_draws Minimum pooled latent-state draws required by the official
#'   classifier.
#' @param reflect Reflection handling passed to \code{classify_hicpotts()}.
#'   When reflection averaging is applied there, the same averaging is applied
#'   to the parameter-based probabilities before comparison.
#' @param method Parameter sensitivity method, \code{"integrated"} by default
#'   or the faster approximate \code{"plugin"} calculation.
#' @param n_draws Number of retained parameter draws evaluated by the
#'   integrated method.
#' @param potts_iterations Number of deterministic mean-field Potts sweeps.
#' @param component_names The three biological labels, in component order.
#'
#' @return An object containing \code{official}, \code{parameter_sensitivity},
#'   per-cell \code{comparison}, an overall \code{summary}, component counts and
#'   the 3-by-3 \code{confusion} table. Every allocation has exactly the three
#'   labels noise, signal and false signal; no additional class is introduced.
#' @examples
#' \dontrun{
#' check <- compare_hicpotts_classifications(
#'     robust_fit,
#'     data = lattice_data, N = 40, dist = "ZINB"
#' )
#' check$summary
#' check$confusion
#' }
#' @seealso \code{\link{classify_hicpotts}},
#'   \code{\link{compute_HMRFHiC_probabilities}}
#' @noRd
compare_hicpotts_classifications <- function(
    fit,
    data,
    N,
    dist = "ZINB",
    use = c("selected", "all"),
    relabel = TRUE,
    min_draws = 100L,
    reflect = c("auto", "always", "never"),
    method = c("integrated", "plugin"),
    n_draws = 200L,
    potts_iterations = 5L,
    component_names = c("noise", "signal", "false signal")
) {
    use <- match.arg(use)
    reflect <- match.arg(reflect)
    method <- match.arg(method)
    dist <- match.arg(dist, c("Poisson", "NB", "ZIP", "ZINB"))
    if (!is.data.frame(data)) data <- as.data.frame(data)
    N <- as.integer(N)
    if (length(N) != 1L || is.na(N) || N < 1L || N * N != nrow(data)) {
        stop("'N' must be positive and satisfy nrow(data) == N^2.")
    }
    if (!is.character(component_names) || length(component_names) != 3L ||
        anyNA(component_names) || any(!nzchar(component_names)) ||
        anyDuplicated(component_names)) {
        stop("'component_names' must contain three distinct non-empty labels.")
    }

    robust <- inherits(fit, "hicpotts_robust_fit")
    fits <- if (robust) {
        if (identical(use, "all")) fit$all_fits else fit$fits
    } else {
        is_single <- is.list(fit) && is.list(fit$chains) &&
            length(fit$chains) == 3L
        if (is_single) list(fit) else fit
    }
    if (!is.list(fits) || !length(fits) ||
        !all(vapply(
            fits, function(x) {
                is.list(x) && is.list(x$chains) && length(x$chains) == 3L
            },
            logical(1)
        ))) {
        stop("'fit' must be a HiCPotts fit, non-empty fit list, or robust fit.")
    }

    official <- classify_hicpotts(
        fit = fit,
        data = data,
        component_names = component_names,
        use = use,
        relabel = relabel,
        min_draws = min_draws,
        reflect = reflect
    )

    ## `iterations` is retained in the lower-level API for compatibility but
    ## indexing is based on each fit's actual length. Supply a valid scalar
    ## derived from the fit rather than asking the user for redundant metadata.
    nominal_iterations <- max(
        2L,
        max(vapply(fits, function(x) nrow(x$chains[[1L]]) - 1L, integer(1)))
    )
    parameter <- compute_HMRFHiC_probabilities(
        data = data,
        chain_betas = fits,
        iterations = nominal_iterations,
        dist = dist,
        consistent_dist = FALSE,
        relabel = relabel,
        N = N,
        potts_iterations = potts_iterations,
        method = method,
        n_draws = n_draws
    )

    pcols <- c("prob1", "prob2", "prob3")
    parameter_probs <- as.matrix(parameter[, pcols, drop = FALSE])
    if (isTRUE(attr(official, "reflection_averaged"))) {
        parameter_probs <- vapply(seq_len(3L), function(k) {
            m <- matrix(parameter_probs[, k], N, N)
            as.vector((m + t(m)) / 2)
        }, numeric(N * N))
        colnames(parameter_probs) <- pcols
    }
    parameter_probs <- parameter_probs / rowSums(parameter_probs)
    parameter_component <- max.col(parameter_probs, ties.method = "first")
    parameter_label <- component_names[parameter_component]

    parameter$prob1 <- parameter_probs[, 1L]
    parameter$prob2 <- parameter_probs[, 2L]
    parameter$prob3 <- parameter_probs[, 3L]
    parameter$map_component <- parameter_component
    parameter$map_label <- parameter_label
    parameter$classification <- factor(parameter_label,
        levels = component_names
    )
    attr(parameter, "classification_source") <- paste0(
        "secondary parameter-", method,
        " allocation with deterministic Potts mean-field smoothing"
    )
    class(parameter) <- c("hicpotts_parameter_sensitivity", "data.frame")

    official_probs <- as.matrix(official[, pcols, drop = FALSE])
    same_label <- official$map_component == parameter_component
    total_variation <- rowSums(abs(official_probs - parameter_probs)) / 2
    comparison <- data.frame(
        official_component = official$map_component,
        official_label = as.character(official$classification),
        parameter_component = parameter_component,
        parameter_label = parameter_label,
        same_label = same_label,
        total_variation = total_variation,
        stringsAsFactors = FALSE
    )
    confusion <- table(
        official = factor(official$map_component,
            levels = seq_len(3L),
            labels = component_names
        ),
        parameter = factor(parameter_component,
            levels = seq_len(3L),
            labels = component_names
        )
    )
    component_summary <- data.frame(
        component = seq_len(3L),
        label = component_names,
        official_n = tabulate(official$map_component, nbins = 3L),
        parameter_n = tabulate(parameter_component, nbins = 3L),
        stringsAsFactors = FALSE
    )
    component_summary$difference_n <-
        component_summary$parameter_n - component_summary$official_n
    summary <- data.frame(
        official_source = attr(official, "classification_source"),
        sensitivity_source = attr(parameter, "classification_source"),
        label_agreement = mean(same_label),
        mean_total_variation = mean(total_variation),
        maximum_total_variation = max(total_variation),
        cells_compared = length(same_label),
        stringsAsFactors = FALSE
    )

    structure(
        list(
            official = official,
            parameter_sensitivity = parameter,
            comparison = comparison,
            summary = summary,
            component_summary = component_summary,
            confusion = confusion
        ),
        class = "hicpotts_classification_sensitivity"
    )
}

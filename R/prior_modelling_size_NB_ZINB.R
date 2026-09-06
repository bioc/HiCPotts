#' @title Size Prior for the Negative Binomial-Type Distributions
#'
#' @description This function computes the prior contribution of the \code{size}
#' parameter (also known as the dispersion parameter) for a specified component
#' in models using a Negative Binomial (NB) or Zero-Inflated Negative Binomial
#' (ZINB) distribution. The prior follows a component-group-specific Gamma
#' distribution.
#'
#' @usage
#' size_prior(size_value, component)
#'
#' @param size_value A numeric value representing the \code{size} (dispersion)
#' parameter of the NB or ZINB distribution. This parameter controls the
#' variance of the distribution, with larger values implying less
#' overdispersion.
#'
#' @param component An integer specifying which component of the mixture model
#' to consider (e.g., \code{1}, \code{2}, or \code{3}). Component 1 has its own
#' prior scale; components 2 and 3 share a common prior.
#'
#' @details In NB and ZINB distributions, the \code{size} parameter controls the
#' variance relative to the mean. Bayesian inference often places a prior on
#' this parameter to regularize its estimation. The prior is placed directly on
#' the positive natural scale of \code{size}.
#'
#' Component 1 uses Gamma(shape=3, rate=1), while components 2 and 3 both use
#' the broader Gamma(shape=2, rate=0.2). Thus the prior means are 3 and 10;
#' the elevated-component prior has standard deviation \eqn{\sqrt{50}} and a
#' long upper tail. This preserves the component-group distinction without
#' imposing a component-2 versus component-3 dispersion ordering.
#'
#' The returned value is the log-density of the Gamma prior at the given
#' \code{size_value}.
#'
#' @return A numeric value representing the log of the Gamma prior density for
#' the \code{size_value} parameter defined for the specified component.
#'
#' @examples
#'
#' # Example: Compute the size prior for size_value = 2.5 in component 1
#' log_prior_comp1 <- size_prior(size_value = 2.5, component = 1)
#' # log_prior_comp1
#'
#' @noRd
#'
# Component 1 is the structurally distinct ZINB/noise component and receives
# Gamma(3, 1). Components 2 and 3 share the broader Gamma(2, 0.2): mean 10,
# SD sqrt(50), and a sufficiently long upper tail for weakly overdispersed
# elevated components. The shared prior cannot itself determine which elevated
# component is signal or false signal.
size_prior <- function(size_value, component) {
    if (!is.numeric(size_value) || length(size_value) != 1L || !is.finite(
        size_value
    ) ||
        size_value <= 0) {
        stop("size_value must be a positive, finite numeric scalar.")
    }
    if (!component %in% seq_len(3)) {
        stop("Invalid component specified. Must be 1, 2, or 3.")
    }
    shape <- c(3, 2, 2)[component]
    rate <- c(1, 0.2, 0.2)[component]
    stats::dgamma(size_value, shape = shape, rate = rate, log = TRUE)
}

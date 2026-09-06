#' @title Proposal Density Function for Model Components
#'
#' @description This function computes the proposal density under a given
#' proposal distribution for a specified model component. In a Bayesian MCMC
#' framework, a proposal distribution is used to generate candidate parameter
#' values in each iteration. By evaluating the log-density of the proposed
#' parameters under this distribution, the MCMC algorithm can determine the
#' acceptance or rejection of the new draw.
#'
#' @usage
#' proposaldensity_combined(params, component)
#'
#' @param params A numeric vector of parameters for which the proposal density
#' is evaluated. The length and meaning of these parameters are determined by
#' the model and the selected component.
#'
#' @param component An integer or identifier specifying which component’s
#' proposal distribution should be used. Each component has its own set of means
#' and standard deviations defining the proposal distribution. Valid components
#' should match those defined within the function (e.g., 1, 2, or 3).
#'
#' @details The function defines three components, each with its own proposal
#' distribution parameters (means and standard deviations) for a set of
#' parameters. It then evaluates the log-density of the \code{params} vector
#' under a multivariate normal distribution (assuming independence between
#' parameters), using the component-specific means and standard deviations.
#'
#' Steps performed by the function:
#' \enumerate{
#'   \item Verifies that the specified \code{component} is valid and corresponds
#'   to one of the defined distributions.
#'   \item Extracts the mean and standard deviation vectors associated with that
#'   component.
#'   \item Checks that the length of \code{params} matches the expected length
#'   of the means and standard deviations.
#'   \item Evaluates the log of the normal probability density function (PDF)
#'   for each parameter using the corresponding mean and standard deviation.
#'   \item Sums these log-densities to produce the overall log-density of the
#'   proposal distribution for the parameter vector.
#' }
#'
#' This function is used internally in our MCMC algorithms where proposals are
#' drawn from component-specific Normal distributions. By returning a
#' log-probability, it simplifies calculations involved in Metropolis-Hastings
#' updates and related MCMC procedures.
#'
#' @return A single numeric value representing the sum of the log proposal
#' densities of all parameters for the specified component.
#'
#' @examples
#'
#' # Example usage:
#' # Suppose we have a parameter vector and want to compute its
#' # proposal density under component 2:
#' params <- c(0, 3, 5, 6, 1) # parameter vector
#' component <- 2
#' log_proposal <- proposaldensity_combined(params, component)
#' # print(log_proposal)
#'
#' @noRd
#'
proposaldensity_combined <- function(params, component) {
    densities <- list(
        component1 = list(means = c(5, 1, 2, 0, 1), sds = c(
            1000, 1000, 1000,
            1000, 1000
        )),
        component2 = list(means = c(300, 2, 4, 5, 1), sds = c(
            5000, 7000, 1000,
            9000, 1000
        )),
        component3 = list(means = c(700, 2, 8, 1, 2), sds = c(
            2000, 5000, 5000,
            2000, 1000
        ))
    )
    key <- paste0("component", component)
    if (!key %in% names(densities)) {
        stop("Invalid component specified: ", component)
    }
    d <- densities[[key]]
    if (length(params) != length(d$means)) {
        stop(
            "The length of params does not match the expected length for ",
            "component ",
            component, "."
        )
    }
    sds <- pmax(d$sds, 1e-6)
    sum(dnorm(params, d$means, sds, log = TRUE))
}

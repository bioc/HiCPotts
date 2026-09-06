#' @title Prior Value for the Potts Model Interaction Parameter
#'
#' @description This function generates a prior value for the interaction
#' parameter in a Potts model from a Beta distribution. The Potts model is a
#' spatial statistical model where the interaction parameter influences the
#' tendency of neighboring sites on a lattice to take on similar states. By
#' drawing the interaction parameter from a Beta distribution, we impose a prior
#' belief on the range and likely values of this parameter.
#'
#' @usage
#' gamma_prior_value()
#'
#' @details
#' The function samples a single random value from a \code{Beta(2, 2)}
#' distribution. This distribution is symmetric around 0.5, weighting middling
#' values of gamma over the extremes, and represents a weakly-informative
#' prior about the direction and strength of spatial clustering.
#'
#' @return A numeric value between \[0,1\] representing the prior draw for the
#' interaction parameter. This is a random draw from the specified Beta
#' distribution.
#'
#' @examples
#' #
#' # Generate a single prior value for the Potts model interaction parameter
#' prior_val <- gamma_prior_value()
#' # prior_val
#'
#' @noRd
#
gamma_prior_value <- function() {
    # draw
    value <- stats::rbeta(1, 2, 2)

    # sanity‐check output
    stopifnot(is.finite(value), value >= 0, value <= 1)
    value
}

#' @title Prediction Function for Combined Model
#'
#' @description This function computes the predicted linear predictor
#' (\eqn{\lambda}) values for a specified model component. It is designed to
#' correct for various sources of bias by incorporating multiple covariates
#' through a log-transformed relationship. The resulting predictions serves as
#' inputs to downstream likelihood or posterior calculations in the Bayesian
#' framework.
#'
#' @usage
#' pred_combined(params, z, x_vars, component, N)
#'
#' @param params A numeric vector of parameters for the model. Typically,
#' \code{params} = \eqn{(a, b, c, d, e)}, where \eqn{a} is the intercept and
#' \eqn{b, c, d, e} are coefficients for the log-transformed covariates.
#'
#' @param z A matrix used to indicate component membership of each observation.
#' The function uses \code{z == component} to select the subset of observations
#' belonging to the specified component.
#'
#' @param x_vars A list of covariates, each element itself a list or vector of
#' length \code{N}, corresponding to different sources of bias or predictors.
#' For example:
#'   \itemize{
#'     \item \code{x_vars[[1]]}: The first covariate (wrapped in a list),
#'     accessible via \code{x_vars[[1]][[1]]}.
#'     \item \code{x_vars[[2]]}, \code{x_vars[[3]]}, \code{x_vars[[4]]}:
#'     Additional covariates in similar nested-list structures.
#'   } Each covariate vector must be of length \code{N}, and the values are
#'   expected to be non-negative since a \eqn{\log(x+1)} transformation is
#'   applied.
#'
#' @param component An integer specifying which component of the mixture or
#' hierarchical model to use when subsetting the data. For the three components,
#' \code{component = 1} will extract the subset of data where \code{z == 1}.
#'
#' @param N An integer specifying the total number of observations. This should
#' match the length of the \code{y} vector and each covariate vector.
#'
#' @details
#' The predicted values (\eqn{\lambda}) are computed as:
#'
#' \deqn{\lambda_i = a + b \cdot \log(x1_i + 1) + c \cdot \log(x2_i + 1) + d
#' \cdot \log(x3_i + 1) + e \cdot \log(x4_i + 1)}
#'
#' for each observation \eqn{i} in the specified component. This log-transform
#' often stabilizes variance and normalizes skewed distributions of covariates.
#'
#' The resulting \eqn{\lambda} values serves as linear predictors (in a Poisson,
#' NB, ZIP, or ZINB model) and as inputs into the likelihood and posterior
#' computation steps.
#'
#' @return A numeric vector containing the predicted values (\eqn{\lambda}) for
#' observations belonging to the specified component.
#'
#' @examples
#' #
#' \donttest{
#' # Assume we have:
#' # N = 3 observations
#' # z indicates component membership for each observation
#' # x_vars is a list of four covariates, each nested in a list:
#' # x_vars[[i]][[1]]
#' N <- 3
#' params <- c(1, 2, 3, 4, 5) # Example parameter values
#' # Example component matrix
#' z <- matrix(c(1, 1, 2, 3, 1, 3, 2, 3, 1), nrow = 3)
#' x_vars <- list(
#'     list(matrix(runif(9, 1, 10), nrow = N)), # first covariate
#'     list(matrix(runif(9, 1, 10), nrow = N)), # second covariate
#'     list(matrix(runif(9, 1, 10), nrow = N)), # third covariate
#'     list(matrix(runif(9, 1, 10), nrow = N)) # fourth covariate
#' )
#'
#'
#' # Get predicted values for component 1:
#' preds_comp1 <- pred_combined(params, z, x_vars, component = 1, N = N)
#' print(preds_comp1)
#' #
#' }
#'
#' @noRd
#'
pred_combined <- function(params, z, x_vars, component, N) {
    # Validate inputs
    if (is.null(x_vars)) {
        stop("x_vars cannot be NULL")
    }
    if (!is.list(x_vars) || length(x_vars) != 4L) {
        stop("x_vars must be a list of four matrices")
    }

    if (length(params) != 5L) {
        stop("params must be a numeric vector of length 5")
    }
    bad <- vapply(x_vars, function(x) {
        if (!is.list(x) || length(x) < 1L) {
            return(TRUE)
        }
        m <- x[[1]]
        !is.matrix(m) || nrow(m) != N || ncol(m) != N
    }, logical(1))
    if (any(bad)) {
        stop("Each x_vars[[k]][[1]] must be an N x N matrix")
    }

    mask <- z == component
    if (!any(mask)) {
        return(numeric(0L))
    } # empty component: caller must handle

    a <- params[1]
    b <- params[2]
    cc <- params[3]
    d <- params[4]
    e <- params[5]

    x1 <- x_vars[[1]][[1]][mask]
    x2 <- x_vars[[2]][[1]][mask]
    x3 <- x_vars[[3]][[1]][mask]
    x4 <- x_vars[[4]][[1]][mask]

    eta <- a + b * log1p(x1) + cc * log1p(x2) + d * log1p(x3) + e * log1p(x4)

    # Numerical guard: clamp the linear predictor before it is exponentiated by
    # the callers (likelihood_combined, pz_123). With raw covariates such as
    # genomic distance, log1p(x) can be ~12, so a moderately large coefficient
    # makes exp(eta) overflow to Inf, yielding -Inf/NaN log-likelihoods that
    # stall the Metropolis step. Capping eta at +/-30 bounds mu to roughly
    # [1e-13, 1e13] without affecting any realistic parameter region.
    pmin(pmax(eta, -30), 30)
}

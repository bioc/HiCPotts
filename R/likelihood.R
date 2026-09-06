#' @title Likelihood Function for a Specified Model Component and Distribution
#'
#' @description This function computes the likelihood contribution of a specific
#' model component using a chosen distribution. It is designed to work within a
#' mixture-modeling or hierarchical modeling framework, where each interaction
#' (such as genomic interactions in a Hi-C experiment) can be modeled by a
#' combination of regression parameters and potentially zero-inflation or
#' overdispersion parameters.
#'
#' @usage
#' likelihood_combined(pred_combined, params, z, y, x_vars, component,
#'                     theta, size, N, dist)
#'
#' @param pred_combined A numeric vector of predictor values obtained from a
#' prediction function. This typically represents the linear predictor
#' \eqn{\lambda} for the given component.
#'
#' @param params A numeric vector of parameters associated with the model’s
#' linear predictors. Typically includes intercepts and regression coefficients
#' related to the covariates in \code{x_vars}.
#'
#' @param z A matrix or array representing the latent state indicators or other
#' structural variables that influence the model. \code{z} helps map each
#' observation to a particular mixture component.
#'
#' @param x_vars A list of covariates used as predictors in the linear model.
#' Each element in the list corresponds to a covariate vector, and all must be
#' of length \code{N}, the number of observations.
#'
#' @param component An integer or factor specifying which component of the
#' mixture model is currently being evaluated. In the 3-component mixture, this
#' is \code{1}, \code{2}, or \code{3}.
#'
#' @param N An integer specifying the number of observations. This should match
#' the length of the response variable \code{y} and the covariates in
#' \code{x_vars}.
#'
#' @param y A numeric vector of observed interaction counts (the response
#' variable). Each element corresponds to one observation (e.g., interaction
#' count between a pair of genomic loci).
#'
#' @param dist A character string specifying the distribution to be used for
#' modeling the interaction counts. Options may include:
#'   \itemize{
#'     \item \code{"Poisson"}: Poisson distribution
#'     \item \code{"NB"}: Negative Binomial distribution
#'     \item \code{"ZIP"}: Zero-Inflated Poisson distribution
#'     \item \code{"ZINB"}: Zero-Inflated Negative Binomial distribution
#'   }
#'
#' @param theta (Optional) A numeric value for the zero-inflation parameter.
#' Required for ZIP and ZINB models. This parameter controls the probability of
#' excess zeros not explained by the Poisson or NB component.
#'
#' @param size (Optional) A numeric value for the size (overdispersion)
#' parameter. Required for NB and ZINB models. This parameter captures variance
#' that exceeds that of a Poisson distribution.
#'
#' @details This function calculates the likelihood for a single component given
#' a set of parameters and covariates. The steps typically involved are:
#' \enumerate{
#'   \item Compute the linear predictor \eqn{\lambda} from the supplied
#'   parameters and covariates. The \code{pred_combined} represents
#'   \eqn{\log(\lambda)}.
#'   \item Depending on \code{dist}, convert the linear predictor into a mean
#'   parameter (e.g., \eqn{\lambda} for Poisson or NB).
#'   \item Compute the likelihood of each observed count \code{y[i]} under the
#'   chosen distribution with the given parameters (\eqn{\lambda}, \eqn{\theta},
#'   \eqn{size}, etc.).
#'   \item For zero-inflated models (ZIP, ZINB), the likelihood incorporates the
#'   probability of an extra zero.
#'   \item Return the computed likelihood values. This is used internally to
#'   evaluate and update model parameters during the MCMC steps.
#' }
#'
#' The function is a building block in a larger modeling framework (MCMC
#' inference for mixture models). While end-users might not call it directly, it
#' enables flexible specification of distributions and model components for
#' complex hierarchical models.
#'
#' @return
#' The function returns the computed likelihood under the specified model setup.
#' It returns a numeric vector of likelihood contributions for each observation.
#'
#' @examples
#'
#' # Example setup
#' N <- 3
#' z <- matrix(c(1, 1, 2, 3, 1, 2, 3, 1, 2), nrow = 3, byrow = TRUE)
#' y <- matrix(c(0, 3, 5, 1, 0, 4, 6, 2, 0), nrow = 3, byrow = TRUE)
#' params <- c(0.1, 0.2, 0.3, 0.4, 0.5)
#' x_vars <- x_vars <- list(
#'     list(matrix(runif(9, 1, 10), nrow = N)), # first covariate
#'     list(matrix(runif(9, 1, 10), nrow = N)), # second covariate
#'     list(matrix(runif(9, 1, 10), nrow = N)), # third covariate
#'     list(matrix(runif(9, 1, 10), nrow = N)) # fourth covariate
#' )
#'
#' theta <- 0.2
#' size <- 10
#' #' # Compute likelihood under a Poisson model, component 1
#' ll_values <- likelihood_combined(
#'     pred_combined = pred_combined,
#'     params = params,
#'     z = z,
#'     y = y,
#'     x_vars = x_vars,
#'     component = 1,
#'     theta,
#'     size,
#'     N = N,
#'     dist = "Poisson"
#' )
#'
#' # print(sum(ll_values)) # sum of likelihood contributions
#'
#' @seealso \code{\link{dpois}}, \code{\link{dnbinom}} for related probability
#' mass functions.
#'
#' @importFrom stats dgamma dnbinom dnorm dpois rbeta rgamma rnorm sd var
#'
#'
#'
#' @noRd
#'
likelihood_combined <- function(
    pred_combined, params, z, y, x_vars, component,
    theta, size, N, dist
) {
    # Subset the data based on the component
    yc <- y[z == component]

    ## early exit: if no observations in component, likelihood contribution is 0
    if (length(yc) == 0L) {
        return(0)
    }

    ## linear predictor vector (length == sum(z == component))
    eta <- pred_combined(params, z, x_vars, component, N)
    lambda <- exp(eta)

    ## guard against mismatched lengths
    if (length(lambda) != length(yc)) {
        template <- paste0(
            "Length mismatch in likelihood_combined: lambda has length ",
            "%d but yc has length %d."
        )
        stop(sprintf(template, length(lambda), length(yc)))
    }

    # Calculate the likelihood based on the specified distribution
    if (component == 1L) {
        ll <- switch(dist,
            "ZIP" = ifelse(yc == 0,
                log(theta + (1 - theta) * exp(-lambda)),
                log1p(-theta) + dpois(yc, lambda = lambda, log = TRUE)
            ),
            "ZINB" = ifelse(yc == 0,
                log(theta + (1 - theta) * dnbinom(0, size = size, mu = lambda)),
                log1p(-theta) + dnbinom(yc,
                    size = size, mu = lambda,
                    log = TRUE
                )
            ),
            "Poisson" = dpois(yc, lambda = lambda, log = TRUE),
            "NB" = dnbinom(yc, size = size, mu = lambda, log = TRUE),
            stop("Invalid distribution specified.")
        )
    } else if (component == 2L || component == 3L) {
        ll <- switch(dist,
            "Poisson" = ,
            "ZIP" = dpois(yc, lambda = lambda, log = TRUE),
            "NB" = ,
            "ZINB" = dnbinom(yc, size = size, mu = lambda, log = TRUE),
            stop("Invalid distribution specified for component > 1.")
        )
    } else {
        stop(sprintf("Invalid component: %s", component))
    }

    ## process_data() supplies the structured complete N x N lattice used by
    ## both the emission model and the Potts field. Each stored cell contributes
    ## once, including both reflected cells when the matrix is symmetric.
    sum(ll)
}

#' Is this contact matrix a mirrored (symmetric) full matrix?
#'
#' Cached on the object's address for the duration of a call stack would be
#' ideal, but a direct check on a square numeric matrix is cheap relative to
#' the density evaluations it guards.
#' @noRd
.hicpotts_symmetric_pairs <- function(y) {
    if (!is.matrix(y)) {
        return(FALSE)
    }
    if (nrow(y) != ncol(y) || nrow(y) < 2L) {
        return(FALSE)
    }
    ## Exact comparison, not all.equal(): mirrored cells are literal copies of
    ## one measured value, so they are bit-identical. Exact equality is both the
    ## correct test and roughly 20x cheaper here. If floating-point drift ever
    ## made a genuinely mirrored matrix compare unequal, the effect is to skip
    ## the halving -- the conservative direction.
    if (anyNA(y)) {
        return(FALSE)
    }
    all(y == t(y))
}

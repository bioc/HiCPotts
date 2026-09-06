#' @title Compute the Posterior for the Model Components
#'
#' @description This function computes the combined posterior value for a
#' specified component in the model utilizing a variety of distributions
#' (Poisson, Negative Binomial, Zero-Inflated Poisson, Zero-Inflated Negative
#' Binomial). It incorporates both the likelihood of the observed data under the
#' given model parameters and the prior distributions for those parameters. When
#' applicable, it also includes a prior on the dispersion parameter (size) used
#' in Negative Binomial or Zero-Inflated Negative Binomial models.
#'
#' @usage
#' posterior_combined(pred_combined, params, z, y, x_vars, component, theta, N,
#'                    use_data_priors, user_fixed_priors, dist, size)
#'
#' @param pred_combined A numeric vector of predictor values generated from a
#' prediction function. Typically represents the linear predictors on a log
#' scale.
#'
#' @param params A numeric vector of parameters for the model. This includes the
#' regression coefficients.
#'
#'
#' @param z A matrix or data structure representing latent variables,
#' components, or states in the model. For example, in the mixture model,
#' \code{z} it indicates component membership of observations.
#'
#' @param y A numeric vector of observed counts or interaction values. Each
#' element corresponds to an observation (e.g., interaction counts between
#' genomic loci).
#'
#' @param x_vars A list of covariates used in the model. Each element
#' corresponds to a predictor variable, and all should be of length \code{N},
#' the number of observations.
#'
#' @param component An integer or factor indicating the component (e.g., mixture
#' component) for which the posterior is being computed. Different components
#' might have different parameter sets or priors.
#'
#' @param theta (Optional) A numeric value representing the zero-inflation
#' parameter for Zero-Inflated models. \eqn{\theta} governs the probability of
#' extra zeros beyond what the Poisson or Negative Binomial distribution would
#' predict.
#'
#' @param N An integer specifying the number of observations. This should match
#' the length of \code{y}.
#'
#' @param use_data_priors A logical value indicating whether to use data-driven
#' priors for each component. If \code{TRUE}, the prior distributions may be
#' estimated from the data.
#'
#' @param user_fixed_priors A logical value indicating whether to use
#' user-specified priors for each component. If \code{TRUE}, a user-supplied
#' prior distribution is applied instead of a data-driven one.
#'
#' @param dist A character string specifying the distribution family. One of:
#'   \itemize{
#'     \item \code{"Poisson"}: Poisson distribution
#'     \item \code{"NB"}: Negative Binomial distribution
#'     \item \code{"ZIP"}: Zero-Inflated Poisson distribution
#'     \item \code{"ZINB"}: Zero-Inflated Negative Binomial distribution
#'   }
#'
#' @param size (Optional) A numeric value for the dispersion parameter in the NB
#' or ZINB distribution. This parameter models overdispersion, allowing the
#' variance to exceed the mean.
#'
#' @details The posterior is computed as: \deqn{\text{Posterior} =
#' \text{Likelihood} + \text{Prior} (+ \text{Size Prior if NB or ZINB})}
#'
#' Steps involved:
#' \enumerate{
#'   \item Compute the likelihood of the data (\code{y}) given the parameters,
#'   latent structure (\code{z}), covariates (\code{x_vars}), and chosen
#'   distribution (\code{dist}) using the provided \code{pred_combined} values.
#'   This is done via the \code{likelihood_combined} function.
#'   \item Compute the prior distribution over the parameters (\code{params}).
#'   Depending on \code{use_data_priors} and \code{user_fixed_priors}, these
#'   priors might be estimated from data or supplied by the user. This step uses
#'   the \code{prior_combined} function.
#'   \item If the distribution is \code{"NB"} or \code{"ZINB"}, incorporate a
#'   separate prior on the dispersion parameter (\code{size}). This step uses
#'   the \code{size_prior} function.
#' }
#'
#' By combining the likelihood and priors, the function returns the posterior
#' value, which is used in the Bayesian inference and MCMC steps to update
#' parameters and perform parameters estimations.
#'
#' @return A numeric value representing the posterior log-probability for the
#' given component under the specified distribution and modeling assumptions.
#'
#' @examples
#'
#' # Example usage of posterior:
#' N <- 5
#' y <- matrix(rpois(N * N, lambda = 5), nrow = N, ncol = N)
#' x_vars <- list(
#'     list(matrix(runif(25), nrow = 5)),
#'     list(matrix(runif(25), nrow = 5)),
#'     list(matrix(runif(25), nrow = 5)),
#'     list(matrix(runif(25), nrow = 5))
#' )
#' params <- c(1, 2, 3, 4, 5)
#' z <- matrix(sample(1:3, N * N, replace = TRUE), nrow = N, ncol = N)
#' preds_comp1 <- pred_combined(params, z, x_vars, component = 1, N = N)
#'
#' # Compute posterior for component 1 using a Poisson distribution
#' post_val <- posterior_combined(
#'     pred_combined = preds_comp1,
#'     params = params,
#'     z = z,
#'     y = y,
#'     x_vars = x_vars,
#'     component = 1,
#'     theta = NULL,
#'     N = N,
#'     use_data_priors = TRUE,
#'     user_fixed_priors = NULL,
#'     dist = "Poisson",
#'     size = NULL
#' )
#' # print(post_val) # Display result
#'
#' @noRd
#'
posterior_combined <- function(
    pred_combined, params, z, y, x_vars, component,
    theta, N,
    use_data_priors, user_fixed_priors, dist, size
) {
    if (!is.character(dist) || length(dist) != 1L) {
        stop("'dist' must be a single character string.")
    }

    lp_size <- 0 # <- move here
    if (dist %in% c("NB", "ZINB")) {
        if (is.null(size) || !is.finite(size) || size <= 0) {
            stop(
                "A positive finite 'size' is required for NB/ZINB ",
                "distributions."
            )
        }
        lp_size <- size_prior(size, component)
    }

    ll <- likelihood_combined(
        pred_combined = pred_combined, params = params, z = z, y = y,
        x_vars = x_vars, component = component, theta = theta, size = size,
        N = N, dist = dist
    )

    lp <- prior_combined(
        params = params, component = component, y = y, x_vars = x_vars, z = z,
        use_data_priors = use_data_priors, user_fixed_priors = user_fixed_priors
    )

    ll + lp + lp_size
}

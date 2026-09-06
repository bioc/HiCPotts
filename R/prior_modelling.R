#' @title Prior Function for Model Components
#'
#' @description This function computes the prior contribution to the
#' log-posterior for the parameters associated with a specified model component.
#' It can use either data-derived empirical-Bayes priors or user-specified fixed
#' priors. By integrating these priors into the Bayesian inference framework,
#' the function helps ensure that parameter estimates remain grounded in
#' reasonable values.
#'
#' @usage
#' prior_combined(params, component, y, x_vars, z, use_data_priors = TRUE,
#'                user_fixed_priors)
#'
#' @param params A numeric vector of parameters for the model. Typically,
#' \code{params} will include a set of parameters associated with the regression
#' coefficients for the model. For example, \code{params} might be \eqn{(a, b,
#' c, d, e)} representing coefficients for the covariates in \code{x_vars}.
#'
#' @param component An integer specifying which component of the model to
#' consider. The model is a mixture with multiple components, this argument
#' indicates which component’s parameters and priors are being evaluated.
#'
#' @param y A numeric vector of observed responses (e.g., interaction counts).
#'   Used with the supplied allocation to construct the empirical-Bayes prior.
#'
#' @param x_vars A list of covariate information. Each element of \code{x_vars}
#' is typically another list or vector representing a single covariate. The same
#' four log1p-transformed covariates as the count model are used in the
#' empirical regression.
#'
#' @param z A matrix or similar structure used to assign observations to
#' components. If \code{z == component}, that observation belongs to the given
#' component. For a direct call to this helper, hyperparameters are calculated
#' from this supplied hard allocation field.
#'
#' @param use_data_priors A logical value indicating whether to use the
#'   data-derived empirical-Bayes priors (\code{TRUE}) or
#'   user-specified fixed priors (\code{FALSE}).
#'
#' @param user_fixed_priors A list of user-defined priors for each component.
#' This argument is only used when \code{use_data_priors = FALSE}. The list
#' should contain entries named \code{component1}, \code{component2}, and
#' \code{component3}, each defining \code{meany}, \code{meanx1}, \code{meanx2},
#' \code{meanx3}, \code{meanx4}, and corresponding \code{sdy}, \code{sdx1},
#' \code{sdx2}, \code{sdx3}, \code{sdx4}.
#'
#' @details
#' This function supports two modes:
#' \enumerate{
#'   \item \strong{Data-derived empirical-Bayes priors:} When
#'         \code{use_data_priors = TRUE}, component-specific Normal centres
#'         and scales are calculated from cells allocated to each component.
#'         The production sampler instead uses soft membership weights during
#'         warm-up and freezes the resulting hyperparameters before retaining
#'         posterior draws. A small-component global fallback and scale floors
#'         avoid singular or degenerate priors.
#'
#'   \item \strong{User-Specified Fixed Priors:} When \code{use_data_priors =
#'   FALSE}, the function uses priors provided by the user through
#'   \code{user_fixed_priors}.
#' }
#'
#' In both modes, the priors are assumed to be normal distributions for each
#' parameter. The function returns the sum of the log of these priors.
#'
#' @return A single numeric value representing the sum of the log-priors for all
#' parameters in \code{params}. This value can be integrated into a Bayesian
#' inference framework to update parameter estimates based on the chosen priors.
#'
#' @examples
#'
#' # Example:
#' N <- 5
#' y <- rpois(N * N, lambda = 5)
#' x_vars <- list(
#'     list(runif(N * N)), # x1
#'     list(rnorm(N * N)), # x2
#'     list(rnorm(N * N)), # x3
#'     list(rnorm(N * N)) # x4
#' )
#' z <- sample(1:3, N * N, replace = TRUE)
#' params <- c(a = 1, b = 0.5, c = -0.2, d = 0.1, e = 0.3)
#'
#' # Using data-driven priors for component 1
#' #' prior_combined(params, 1, y, x_vars, z, TRUE, NULL)
#'
#' # Using user-specified fixed priors for component 2
#' # user_fixed_priors <- list(
#' #  component1 = list(
#' #    meany = 5, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0,
#' #    sdy = 1, sdx1 = 1, sdx2 = 1, sdx3 = 1, sdx4 = 1
#' #  ),
#' #  component2 = list(
#' #    meany = 10, meanx1 = 1, meanx2 = 1, meanx3 = 1, meanx4 = 1,
#' #    sdy = 2, sdx1 = 2, sdx2 = 2, sdx3 = 2, sdx4 = 2
#' #  ),
#' #  component3 = list(
#' #    meany = 2, meanx1 = 2, meanx2 = 2, meanx3 = 2, meanx4 = 2,
#' #    sdy = 1, sdx1 = 1, sdx2 = 1, sdx3 = 1, sdx4 = 1
#' #  )
#' # )
#' #
#' # prior_val_fixed <- prior_combined(
#' #  params = params,
#' #  component = 2,
#' #  y = y,
#' #  x_vars = x_vars,
#' #  z = z,
#' #  use_data_priors = FALSE,
#' #  user_fixed_priors = user_fixed_priors
#' # )
#'
#' @importFrom stats dnorm sd
#'
#' @noRd
#'
# Prior for the sources-of-biases regression parameters. Under iterated
# empirical Bayes, the sampler updates soft-membership hyperparameters only
# during warm-up and then freezes them for all retained transitions. Direct
# calls obtain deterministic hard-allocation hyperparameters through the
# compiled helper below. The separate joint prior
# that couples component 3 to component 1 remains part of the sampler target.
prior_combined <- function(
    params, component, y, x_vars, z,
    use_data_priors = TRUE, user_fixed_priors = NULL
) {
    if (length(params) != 5L) {
        stop("params must be a numeric vector of length 5.")
    }
    if (!component %in% seq_len(3)) {
        stop("Invalid component specified. Must be 1, 2, or 3.")
    }

    a <- params[1]
    b <- params[2]
    c_ <- params[3]
    d <- params[4]
    e <- params[5]
    epsilon <- 1e-6 # floor on sd to avoid dnorm(..., sd = 0)

    if (isTRUE(use_data_priors)) {
        priors <- .hicpotts_iterated_eb_priors(y, x_vars, z)[[
            paste0("component", component)
        ]]
        meany <- priors$meany
        sdy <- max(priors$sdy, epsilon)
        meanx1 <- priors$meanx1
        sdx1 <- max(priors$sdx1, epsilon)
        meanx2 <- priors$meanx2
        sdx2 <- max(priors$sdx2, epsilon)
        meanx3 <- priors$meanx3
        sdx3 <- max(priors$sdx3, epsilon)
        meanx4 <- priors$meanx4
        sdx4 <- max(priors$sdx4, epsilon)
    } else {
        ## ---- user-fixed priors
        ## --------------------------------------------------
        if (is.null(user_fixed_priors)) {
            stop(
                "user_fixed_priors must be provided when use_data_priors = ",
                "FALSE."
            )
        }
        priors <- user_fixed_priors[[paste0("component", component)]]
        if (is.null(priors)) {
            stop(
                "user_fixed_priors must contain entry 'component", component,
                "'."
            )
        }
        req <- c(
            "meany", "meanx1", "meanx2", "meanx3", "meanx4",
            "sdy", "sdx1", "sdx2", "sdx3", "sdx4"
        )
        missing_fields <- setdiff(req, names(priors))
        if (length(missing_fields)) {
            stop(
                "user_fixed_priors$component", component,
                " is missing fields: ", paste(missing_fields, collapse = ", ")
            )
        }
        meany <- priors$meany
        sdy <- max(priors$sdy, epsilon)
        meanx1 <- priors$meanx1
        sdx1 <- max(priors$sdx1, epsilon)
        meanx2 <- priors$meanx2
        sdx2 <- max(priors$sdx2, epsilon)
        meanx3 <- priors$meanx3
        sdx3 <- max(priors$sdx3, epsilon)
        meanx4 <- priors$meanx4
        sdx4 <- max(priors$sdx4, epsilon)
    }

    dnorm(a, meany, sdy, log = TRUE) +
        dnorm(b, meanx1, sdx1, log = TRUE) +
        dnorm(c_, meanx2, sdx2, log = TRUE) +
        dnorm(d, meanx3, sdx3, log = TRUE) +
        dnorm(e, meanx4, sdx4, log = TRUE)
}

#' @noRd
.hicpotts_iterated_eb_priors <- function(y, x_vars, z) {
    if (!is.list(x_vars) || length(x_vars) != 4L) {
        stop("x_vars must contain four covariate matrices.")
    }
    covariates <- lapply(x_vars, function(value) {
        matrix_value <- if (is.list(value)) value[[1L]] else value
        if (!is.matrix(matrix_value) || !is.numeric(matrix_value)) {
            stop("Each covariate must contain a numeric matrix.")
        }
        matrix_value
    })
    reference_dim <- dim(covariates[[1L]])
    if (!all(vapply(covariates, function(value) {
        identical(dim(value), reference_dim)
    }, logical(1)))) {
        stop("All empirical-Bayes covariates must have identical dimensions.")
    }
    as_matching_matrix <- function(value, name) {
        if (is.matrix(value) && identical(dim(value), reference_dim)) {
            return(value)
        }
        if (length(value) != prod(reference_dim)) {
            stop(name, " must match the covariate matrix dimensions.")
        }
        matrix(value, nrow = reference_dim[1L], ncol = reference_dim[2L])
    }
    y <- as_matching_matrix(y, "y")
    z <- as_matching_matrix(z, "z")
    .hicpotts_iterated_eb_priors_cpp(
        y, covariates[[1L]], covariates[[2L]], covariates[[3L]],
        covariates[[4L]], z
    )
}

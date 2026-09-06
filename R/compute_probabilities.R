#' @title Compute HiCPotts Probabilities of Assigning an Interaction to Each
#' Component
#'
#' @description
#' Computes a secondary parameter-based allocation of genomic interactions to
#' each of three HiCPotts components. The official classification is returned
#' by \code{classify_hicpotts()} from sampled latent-state frequencies. This
#' function is retained for parameter-plus-Potts sensitivity analysis and
#' backwards compatibility.
#'
#' Posterior mean regression parameters are extracted from MCMC output and used
#' to evaluate the component-wise likelihood of each observed interaction count.
#' The resulting probabilities are normalized to sum to 1 for each interaction.
#'
#' Optionally, the MCMC output can be relabeled before posterior summaries are
#' computed. This helps reduce label-switching effects across MCMC draws and
#' chains.
#'
#' @usage
#' compute_HMRFHiC_probabilities(
#'   data, chain_betas, iterations,
#'   dist = "ZINB",
#'   max_interactions = NA_integer_,
#'   consistent_dist = FALSE,
#'   relabel = FALSE,
#'   N = NULL,
#'   potts_iterations = 5L,
#'   method = c("integrated", "plugin"),
#'   n_draws = 200L,
#'   component_definition = hicpotts_component_definition()
#' )
#'
#' @param data A \code{data.frame}, or block data accepted by
#'   \code{combine_hicpotts_blocks()}, with required columns:
#'   \code{start}, \code{end}, \code{interactions}, \code{GC}, \code{TES},
#'   and \code{ACC}.
#' @param chain_betas A list of MCMC chain results, robust fit, or block-aware
#'   result from \code{combine_hicpotts_blocks()}.
#' @param iterations Total number of MCMC iterations.
#' @param dist Distribution used for count modeling. One of
#'   \code{"Poisson"}, \code{"NB"}, \code{"ZIP"}, or \code{"ZINB"}.
#' @param max_interactions Optional cap on interactions before density
#'   evaluation. Defaults to \code{NA}, meaning no cap.
#' @param consistent_dist Logical; controls the family used for components 2
#'   and 3. Defaults to \code{FALSE}, which matches the fitted model: the
#'   sampler applies structural zero inflation to component 1 only, so
#'   components 2 and 3 use \code{"Poisson"} or \code{"NB"} for their
#'   emission density, matching the likelihood the sampler evaluates.
#'   \code{TRUE} gives all three components the
#'   zero-inflated family and therefore does NOT correspond to any model this
#'   package fits; with equal means it drives every observed zero to a flat
#'   1/3, 1/3, 1/3 instead of the decisive component-1 call the fitted
#'   likelihood implies. It is retained only for backwards comparison.
#' @param method Either \code{"integrated"} (default) or \code{"plugin"}.
#'   \code{"integrated"} computes normalised component probabilities for each
#'   retained posterior draw and averages them, which is the integrated
#'   posterior membership probability. \code{"plugin"} reproduces the legacy
#'   behaviour: average the parameters first, then evaluate a single density.
#'   Since \eqn{p(z \mid y, E[\theta]) \neq E[p(z \mid y, \theta)]}, the
#'   plug-in values understate uncertainty and should be described as
#'   approximate scores, not posterior probabilities.
#' @param n_draws Integer; number of posterior draws to thin to when
#'   \code{method = "integrated"}. Cost is linear in this value.
#' @param component_definition Three-row data frame used to label the
#'   returned probability columns, with columns \code{component} (1, 2,
#'   3), \code{label} and \code{definition}. The default is the package's
#'   canonical definition, equivalent to
#'   \preformatted{data.frame(
#'   component  = 1:3,
#'   label      = c("noise", "signal", "false signal"),
#'   definition = c("low-mean noise",
#'                  "true interaction with an unrestricted
#'                   covariate-response pattern",
#'                  "elevated noise with a covariate-response pattern
#'                   resembling component 1")
#' )}
#'   It is built internally, which is why the usage above shows a function
#'   call. Supplying your own wording renames the reported classes only; it
#'   does not change the fitted component identities.
#' @param relabel Logical; if \code{TRUE}, relabel each fitted chain result
#'   before posterior means are computed.
#' @param N Optional integer lattice dimension. If \code{data} represents an
#'   \eqn{N \times N} lattice in the same row-order convention
#'   \code{process_data()} uses (i.e. \code{nrow(data) == N^2}), supplying
#'   \code{N} enables a Potts/neighbour spatial smoothing term in the returned
#'   probabilities, matching the neighbour coupling the sampler applies
#'   when it draws the latent states during MCMC. If omitted, the
#'   function returns per-cell-likelihood-only
#'   probabilities and emits a warning (see Details).
#' @param potts_iterations Number of mean-field Potts smoothing sweeps to run
#'   when \code{N} is supplied. Default 5.
#'
#' @details
#' When \code{N} is supplied and \code{nrow(data) == N^2}, the per-cell
#' log-densities are reshaped into \eqn{N \times N} matrices (the same
#' column-major row-order convention as \code{process_data()}) and a small
#' number of mean-field Potts smoothing sweeps are run: each cell's smoothed
#' log-probability is its own per-cell log-density plus the fitted
#' \code{gamma} (posterior mean, extracted the same way as theta/size) times
#' the current-estimate same-label probability mass among its four lattice
#' neighbours, iterated \code{potts_iterations} times so the spatial term
#' propagates a short distance -- the mean-field (deterministic, iterated)
#' analogue of the neighbour-agreement term the sampler applies by
#' drawing latent states during MCMC. If \code{N} is not supplied (or
#' does not satisfy
#' \code{N^2 == nrow(data)}), the function returns the spatial-free
#' probabilities and emits an explicit warning rather than silently guessing.
#' The function:
#' \enumerate{
#'   \item checks the required input columns,
#'   \item optionally relabels the fitted MCMC output,
#'   \item discards the first half of iterations as burn-in,
#'   \item extracts posterior means for each component,
#'   \item evaluates the component-wise count density,
#'   \item normalizes these values into posterior probabilities.
#' }
#' No hard probability threshold is applied. The returned \code{prob1},
#' \code{prob2}, and \code{prob3} columns are intended for users to threshold or
#' otherwise interpret according to their scientific application.
#'
#' @return A \code{data.frame} containing the original columns plus:
#' \itemize{
#'   \item \code{prob1}
#'   \item \code{prob2}
#'   \item \code{prob3}
#' }
#'
#'
#' @section Working with sub-blocks of a large map:
#' This function summarises one contact map. It cannot span sub-blocks:
#' chains fitted to different blocks have different data and different
#' frozen priors, so they are refused rather than pooled, and passing a
#' whole-map \code{data} frame to a single block's fit fails the
#' one-row-per-lattice-cell check.
#'
#' Run it per block and stitch the results, which lets every block keep
#' the chains its own mode screen accepted -- block 1 may be summarised
#' from chains 1 and 3 while block 2 uses chains 1 and 2:
#' \preformatted{parts <- list()
#' for (j in seq_along(fits)) {
#'     parts[[j]] <- compute_HMRFHiC_probabilities(
#'         blocks[[j]], fits[[j]], iterations)
#'     parts[[j]]$block <- names(fits)[j]
#' }
#' combined <- do.call(rbind, parts)}
#' \code{\link{combine_hicpotts_blocks}()} does this with the checks that
#' matter (matching lengths, one row per cell, non-overlapping blocks) and
#' records the chains each block used. Note that component labels are only
#' aligned within a block, not across blocks.
#' @examples
#' set.seed(4921)
#'
#' sim_data <- data.frame(
#'     start = c(1e6, 2e6, 3e6, 4e6),
#'     end = c(2e6, 3e6, 4e6, 5e6),
#'     interactions = c(0, 5, 10, 20),
#'     GC = c(0.40, 0.45, 0.50, 0.55),
#'     TES = c(0.10, 0.20, 0.30, 0.40),
#'     ACC = c(0.30, 0.40, 0.50, 0.60)
#' )
#'
#' iterations <- 10L
#' n_draw <- iterations + 1L
#'
#' make_chain <- function(intercept) {
#'     cbind(
#'         rnorm(n_draw, intercept, 0.05),
#'         rnorm(n_draw, -0.01, 0.01),
#'         rnorm(n_draw, 0.01, 0.01),
#'         rnorm(n_draw, 0.01, 0.01),
#'         rnorm(n_draw, 0.01, 0.01)
#'     )
#' }
#'
#' chain_betas <- list(
#'     list(
#'         chains = list(
#'             make_chain(-1),
#'             make_chain(0),
#'             make_chain(1)
#'         ),
#'         theta = runif(n_draw, 0.05, 0.20),
#'         size = matrix(10, nrow = 3, ncol = n_draw),
#'         gamma = runif(n_draw, 0.5, 0.8)
#'     )
#' )
#'
#' prob_res <- compute_HMRFHiC_probabilities(
#'     data = sim_data,
#'     chain_betas = chain_betas,
#'     iterations = iterations,
#'     dist = "ZINB",
#'     relabel = TRUE
#' )
#'
#' head(prob_res)
#'
#' @seealso \code{\link{relabel_hicpotts}}, \code{\link{dpois}},
#'   \code{\link{dnbinom}}
#'
#' @export
compute_HMRFHiC_probabilities <- function(
    data = NULL,
    chain_betas = NULL,
    iterations = NULL,
    dist = "ZINB",
    max_interactions = NA_integer_,
    consistent_dist = FALSE,
    relabel = FALSE,
    N = NULL,
    potts_iterations = 5L,
    method = c("integrated", "plugin"),
    n_draws = 200L,
    component_definition = hicpotts_component_definition()
) {
    method <- match.arg(method)
    if (inherits(chain_betas, "hicpotts_block_fit")) {
        return(.hicpotts_probability_blocks(
            fit = chain_betas, data = data, iterations = iterations,
            dist = dist, max_interactions = max_interactions,
            consistent_dist = consistent_dist, relabel = relabel,
            N = N, potts_iterations = potts_iterations, method = method,
            n_draws = n_draws,
            component_definition = component_definition
        ))
    }
    if (is.null(data)) stop("'data' must be supplied.", call. = FALSE)
    .check_required_columns(data)

    required_definition <- c("component", "label", "definition")
    if (!is.data.frame(component_definition) ||
        !all(required_definition %in% names(component_definition)) ||
        nrow(component_definition) != 3L ||
        !setequal(component_definition$component, seq_len(3L)) ||
        anyDuplicated(component_definition$component) ||
        anyNA(component_definition[, required_definition, drop = FALSE]) ||
        any(!nzchar(as.character(component_definition$label))) ||
        anyDuplicated(as.character(component_definition$label))) {
        message <- paste0(
            "'component_definition' must contain one unique, non-missing ",
            "component, label and definition for each of components 1:3."
        )
        stop(message, call. = FALSE)
    }
    component_definition <- component_definition[
        match(seq_len(3L), component_definition$component),
        required_definition,
        drop = FALSE
    ]

    dist <- match.arg(dist, c("Poisson", "NB", "ZIP", "ZINB"))

    if (!is.logical(relabel) || length(relabel) != 1L || is.na(relabel)) {
        stop("'relabel' must be TRUE or FALSE.")
    }

    if (inherits(chain_betas, "hicpotts_robust_fit")) {
        chain_betas <- chain_betas$fits
    }
    is_single_fit <- is.list(chain_betas) && is.list(chain_betas$chains) &&
        length(chain_betas$chains) == 3L
    if (is_single_fit) chain_betas <- list(chain_betas)
    valid_fits <- is.list(chain_betas) && length(chain_betas) > 0L &&
        all(vapply(
            chain_betas, function(x) {
                is.list(x) && is.list(x$chains) && length(x$chains) == 3L
            },
            logical(1)
        ))
    if (!valid_fits) {
        stop("'chain_betas' must be a HiCPotts fit, chain list, or robust fit.",
            call. = FALSE
        )
    }

    if (relabel) {
        chain_betas <- relabel_hicpotts(chain_betas)
    }

    mydata <- as.data.frame(data)

    interactions <- mydata$interactions
    distance <- log1p(.hicpotts_genomic_distance(mydata))
    GC <- log1p(mydata$GC)
    TES <- log1p(mydata$TES)
    ACC <- log1p(mydata$ACC)

    if (!is.numeric(iterations) || length(iterations) != 1L || iterations < 2) {
        stop("iterations must be a positive scalar of at least 2.")
    }

    ## Retained draws are derived from each object's ACTUAL length, never from
    ## the requested `iterations`. MCSE early stopping is enabled by default, so
    ## a valid fit routinely holds fewer rows than were requested; indexing by
    ## `iterations` then raised "subscript out of bounds" on a perfectly good
    ## fit. `iterations` is now only validated, not used for indexing, and is
    ## retained purely for backwards compatibility of the call signature.
    retained <- function(n) {
        if (!is.finite(n) || n < 2L) {
            stop(
                "A fitted chain holds fewer than two draws; cannot summarise ",
                "it."
            )
        }
        (floor(n / 2L) + 1L):n
    }

    extract_mean <- function(chain_list, col) {
        vapply(
            chain_list,
            function(chain) {
                vapply(
                    seq_len(3L),
                    function(j) {
                        m <- chain[["chains"]][[j]]
                        mean(m[retained(nrow(m)), col])
                    },
                    numeric(1)
                )
            },
            numeric(3)
        )
    }

    intercept_mat <- extract_mean(chain_betas, 1L)
    distance_mat <- extract_mean(chain_betas, 2L)
    GC_mat <- extract_mean(chain_betas, 3L)
    TES_mat <- extract_mean(chain_betas, 4L)
    ACC_mat <- extract_mean(chain_betas, 5L)

    comp_mean <- function(mat, comp) mean(mat[comp, ])

    intercept <- vapply(
        seq_len(3), function(c) comp_mean(intercept_mat, c),
        numeric(1)
    )
    distance_b <- vapply(
        seq_len(3), function(c) comp_mean(distance_mat, c),
        numeric(1)
    )
    GC_b <- vapply(seq_len(3), function(c) comp_mean(GC_mat, c), numeric(1))
    TES_b <- vapply(seq_len(3), function(c) comp_mean(TES_mat, c), numeric(1))
    ACC_b <- vapply(seq_len(3), function(c) comp_mean(ACC_mat, c), numeric(1))

    theta <- if (dist %in% c("ZIP", "ZINB")) {
        mean(vapply(
            chain_betas,
            function(ch) mean(ch[["theta"]][retained(length(ch[["theta"]]))]),
            numeric(1)
        ))
    } else {
        NULL
    }

    overdisp <- if (dist %in% c("NB", "ZINB")) {
        vapply(
            seq_len(3),
            function(c) {
                mean(vapply(
                    chain_betas,
                    function(ch) {
                        mean(ch[["size"]][c, retained(ncol(ch[[
                            "size"
                        ]]))])
                    },
                    numeric(1)
                ))
            },
            numeric(1)
        )
    } else {
        rep(NA_real_, 3)
    }

    interactions_eval <- if (is.na(max_interactions)) {
        interactions
    } else {
        pmin(interactions, as.integer(max_interactions))
    }

    log_density <- function(y, lambda, theta, overdisp, use_dist) {
        switch(use_dist,
            "Poisson" = dpois(y, lambda = lambda, log = TRUE),
            "NB" = dnbinom(y, size = overdisp, mu = lambda, log = TRUE),
            "ZIP" = ifelse(
                y == 0,
                log(theta + (1 - theta) * exp(-lambda)),
                log1p(-theta) + dpois(y, lambda = lambda, log = TRUE)
            ),
            "ZINB" = ifelse(
                y == 0,
                log(theta + (1 - theta) * dnbinom(0,
                    size = overdisp,
                    mu = lambda
                )),
                log1p(-theta) + dnbinom(y,
                    size = overdisp, mu = lambda,
                    log = TRUE
                )
            ),
            stop("Invalid distribution specified.")
        )
    }

    dist_23 <- if (consistent_dist) {
        dist
    } else if (dist %in% c("Poisson", "ZIP")) {
        "Poisson"
    } else {
        "NB"
    }

    ## ---- incorporate the Potts/neighbour spatial term -----------------------
    n_obs <- nrow(mydata)
    use_spatial <- !is.null(N)
    if (use_spatial) {
        if (!is.numeric(N) || length(N) != 1L || N < 1L || N * as.integer(
            N
        ) != n_obs) {
            warning(
                "compute_HMRFHiC_probabilities: 'N' was supplied but N^2 != ",
                "nrow(data); ",
                "falling back to spatial-free (per-cell-likelihood-only) ",
                "probabilities."
            )
            use_spatial <- FALSE
        }
    } else {
        warning(
            "compute_HMRFHiC_probabilities: 'N' not supplied, so no lattice ",
            "adjacency is available -- returning per-cell-likelihood-only ",
            "probabilities with no spatial/Potts term. Supply N (matching ",
            "process_data()'s row-order convention) for spatially-smoothed ",
            "probabilities consistent with the neighbour coupling the ",
            "sampler applies during MCMC."
        )
    }

    softmax3 <- function(a, b, c) {
        m <- pmax(a, pmax(b, c))
        e1 <- exp(a - m)
        e2 <- exp(b - m)
        e3 <- exp(c - m)
        s <- e1 + e2 + e3
        list(p1 = e1 / s, p2 = e2 / s, p3 = e3 / s)
    }

    ## Zero-padded neighbour shift (matches Neighbours_combined()'s border
    ## convention: missing neighbours contribute 0).
    shift0 <- function(m, di, dj, Nn) {
        out <- matrix(0, Nn, Nn)
        i_src <- seq_len(Nn) - di
        j_src <- seq_len(Nn) - dj
        i_ok <- i_src >= 1 & i_src <= Nn
        j_ok <- j_src >= 1 & j_src <= Nn
        out[i_ok, j_ok] <- m[i_src[i_ok], j_src[j_ok]]
        out
    }

    ## ---- probabilities for ONE parameter draw
    ## -------------------------------- `b` is 5 x 3 (coefficients x component);
    ## theta/size/gamma are that draw's values. Factored out so the plug-in path
    ## and the draw-by-draw path share exactly one implementation.
    probs_for <- function(b, theta_d, overdisp_d, gamma_d) {
        eta_c <- function(c) {
            b[1L, c] + b[2L, c] * distance + b[3L, c] * GC +
                b[4L, c] * TES + b[5L, c] * ACC
        }
        lp1 <- log_density(
            interactions_eval, exp(eta_c(1L)), theta_d,
            overdisp_d[1], dist
        )
        lp2 <- log_density(
            interactions_eval, exp(eta_c(2L)), theta_d,
            overdisp_d[2], dist_23
        )
        lp3 <- log_density(
            interactions_eval, exp(eta_c(3L)), theta_d,
            overdisp_d[3], dist_23
        )

        if (!use_spatial) {
            sm <- softmax3(lp1, lp2, lp3)
            return(list(
                p1 = as.vector(sm$p1), p2 = as.vector(sm$p2),
                p3 = as.vector(sm$p3)
            ))
        }

        Ni <- as.integer(N)
        lp1m <- matrix(lp1, Ni, Ni)
        lp2m <- matrix(lp2, Ni, Ni)
        lp3m <- matrix(lp3, Ni, Ni)
        init <- softmax3(lp1m, lp2m, lp3m)
        p1m <- init$p1
        p2m <- init$p2
        p3m <- init$p3
        for (mf_iter in seq_len(potts_iterations)) {
            neigh1 <- shift0(p1m, -1, 0, Ni) + shift0(p1m, 1, 0, Ni) + shift0(
                p1m, 0, -1, Ni
            ) + shift0(p1m, 0, 1, Ni)
            neigh2 <- shift0(p2m, -1, 0, Ni) + shift0(p2m, 1, 0, Ni) + shift0(
                p2m, 0, -1, Ni
            ) + shift0(p2m, 0, 1, Ni)
            neigh3 <- shift0(p3m, -1, 0, Ni) + shift0(p3m, 1, 0, Ni) + shift0(
                p3m, 0, -1, Ni
            ) + shift0(p3m, 0, 1, Ni)
            sm <- softmax3(
                lp1m + gamma_d * neigh1, lp2m + gamma_d * neigh2,
                lp3m + gamma_d * neigh3
            )
            p1m <- sm$p1
            p2m <- sm$p2
            p3m <- sm$p3
        }
        list(p1 = as.vector(p1m), p2 = as.vector(p2m), p3 = as.vector(p3m))
    }

    gamma_mean <- if (use_spatial) {
        mean(vapply(
            chain_betas,
            function(ch) mean(ch[["gamma"]][retained(length(ch[["gamma"]]))]),
            numeric(1)
        ))
    } else {
        NA_real_
    }

    if (identical(method, "plugin")) {
        ## Legacy behaviour: average the PARAMETERS, then evaluate one density.
        b_mean <- rbind(intercept, distance_b, GC_b, TES_b, ACC_b)
        out <- probs_for(b_mean, theta, overdisp, gamma_mean)
    } else {
        ## ---- M1: integrate over the posterior, draw by draw
        ## -------------------- p(z | y, E[params]) != E[p(z | y, params)].
        ## Averaging parameters first and evaluating a single density
        ## understates uncertainty, and the gap matters precisely where this
        ## model lives: exponentiated linear predictors, negative-binomial
        ## dispersion, zero inflation and label uncertainty. Probabilities are
        ## now computed for each retained draw and averaged, which is the
        ## integrated posterior membership probability the output always claimed
        ## to be.
        ##
        ## Draws are thinned to `n_draws` across all chains to bound the cost;
        ## evaluation is O(n_draws) times the plug-in cost.
        specs <- do.call(rbind, lapply(seq_along(chain_betas), function(ci) {
            ch <- chain_betas[[ci]]
            idx <- retained(nrow(ch[["chains"]][[1L]]))
            data.frame(chain = ci, row = idx)
        }))
        n_draws <- as.integer(n_draws)
        if (length(n_draws) != 1L || is.na(n_draws) || n_draws < 1L) {
            stop("'n_draws' must be a positive integer.")
        }
        if (nrow(specs) > n_draws) {
            specs <- specs[round(seq(1, nrow(specs), length.out = n_draws)), ,
                drop = FALSE
            ]
        }

        acc1 <- acc2 <- acc3 <- numeric(nrow(mydata))
        for (r in seq_len(nrow(specs))) {
            ch <- chain_betas[[specs$chain[r]]]
            k <- specs$row[r]
            b_draw <- vapply(
                seq_len(3L),
                function(j) as.numeric(ch[["chains"]][[j]][k, seq_len(5L)]),
                numeric(5)
            )
            theta_d <- if (dist %in% c("ZIP", "ZINB")) {
                ch[["theta"]][
                    k
                ]
            } else {
                theta
            }
            overdisp_d <- if (dist %in% c("NB", "ZINB")) {
                as.numeric(ch[[
                    "size"
                ]][, k])
            } else {
                overdisp
            }
            gamma_d <- if (use_spatial) ch[["gamma"]][k] else NA_real_
            p <- probs_for(b_draw, theta_d, overdisp_d, gamma_d)
            acc1 <- acc1 + p$p1
            acc2 <- acc2 + p$p2
            acc3 <- acc3 + p$p3
        }
        nd <- nrow(specs)
        out <- list(p1 = acc1 / nd, p2 = acc2 / nd, p3 = acc3 / nd)
        attr(mydata, "n_posterior_draws") <- nd
    }

    mydata$prob1 <- out$p1
    mydata$prob2 <- out$p2
    mydata$prob3 <- out$p3
    attr(mydata, "probability_method") <- method
    attr(mydata, "threshold_policy") <-
        "none in package; user chooses any downstream probability threshold"
    attr(mydata, "classification_role") <-
        paste0("secondary parameter-based sensitivity; official labels come ",
            "from classify_hicpotts()")
    attr(mydata, "component_definition") <- component_definition
    attr(mydata, "probability_components") <- stats::setNames(
        component_definition$label, paste0(
            "prob",
            component_definition$component
        )
    )

    mydata
}

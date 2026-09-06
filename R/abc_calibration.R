.hicpotts_prepare_shared_abc_epsilon <- function(
    epsilon, N, gamma_prior_shape1 = 1, gamma_prior_shape2 = 1,
    abc_epsilon_quantile = 0.10, abc_potts_sweeps = 0L,
    abc_sim_reps = 4L
) {
    if (!is.null(epsilon)) {
        if (!is.numeric(epsilon) || length(epsilon) != 1L ||
            !is.finite(epsilon) || epsilon <= 0) {
            stop("epsilon must be NULL or a finite positive ABC bandwidth.",
                call. = FALSE
            )
        }
        return(list(
            epsilon = as.numeric(epsilon),
            user_supplied = TRUE,
            shared_calibration = FALSE,
            calibration_reps = 0L,
            source = "user"
        ))
    }

    calibration <- withr::with_seed(
        20260721L,
        .hicpotts_calibrate_abc_epsilon_cpp(
            N = as.integer(N),
            gamma_prior_shape1 = gamma_prior_shape1,
            gamma_prior_shape2 = gamma_prior_shape2,
            abc_epsilon_quantile = abc_epsilon_quantile,
            abc_potts_sweeps_arg = as.integer(abc_potts_sweeps),
            abc_sim_reps = as.integer(abc_sim_reps),
            abc_calibration_reps = 60L
        )
    )
    calibration$user_supplied <- FALSE
    calibration$shared_calibration <- TRUE
    calibration$calibration_reps <- calibration$abc_calibration_reps
    calibration$source <- "shared_prior_predictive"
    calibration
}

.hicpotts_record_shared_abc_calibration <- function(fits, calibration) {
    lapply(fits, function(fit) {
        attr(fit$gamma, "abc_epsilon_user_supplied") <-
            isTRUE(calibration$user_supplied)
        attr(fit$gamma, "abc_epsilon_shared_calibration") <-
            isTRUE(calibration$shared_calibration)
        attr(fit$gamma, "abc_calibration_reps_executed") <- 0L
        attr(fit$gamma, "abc_shared_calibration_reps") <-
            as.integer(calibration$calibration_reps)
        attr(fit$gamma, "abc_epsilon_source") <- calibration$source
        fit$sampler_settings$abc_epsilon <- calibration$epsilon
        fit$sampler_settings$abc_epsilon_user_supplied <-
            isTRUE(calibration$user_supplied)
        fit$sampler_settings$abc_epsilon_shared_calibration <-
            isTRUE(calibration$shared_calibration)
        fit$sampler_settings$abc_calibration_reps_executed <- 0L
        fit$sampler_settings$abc_shared_calibration_reps <-
            as.integer(calibration$calibration_reps)
        fit$sampler_settings$abc_epsilon_source <- calibration$source
        fit
    })
}

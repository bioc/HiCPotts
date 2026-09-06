#' Canonical HiCPotts component identities
#'
#' Returns the single package-level definition used by fitting metadata,
#' relabelling and classification. Keeping these labels in one place prevents
#' user-facing names from drifting away from the labelled posterior model.
#'
#' @return A three-row data frame containing component number, short label and
#'   biological definition.
#' @examples
#' hicpotts_component_definition()
#' @noRd
hicpotts_component_definition <- function() {
    data.frame(
        component = seq_len(3L),
        label = c("noise", "signal", "false signal"),
        definition = c(
            "low-mean noise",
            "true interaction with an unrestricted covariate-response pattern",
            paste(
                "elevated noise with a covariate-response pattern resembling",
                "component 1"
            )
        ),
        stringsAsFactors = FALSE
    )
}

.hicpotts_validate_fit_inputs <- function(
    N, y, x_vars, iterations, gamma_prior, theta_start = NULL,
    size_start = NULL, dist = "ZIP", use_data_priors = TRUE,
    user_fixed_priors = NULL
) {
    if (!is.numeric(N) || length(N) != 1L || !is.finite(N) ||
        N < 2 || N != as.integer(N)) {
        stop("'N' must be one integer of at least two.", call. = FALSE)
    }
    if (!is.numeric(iterations) || length(iterations) != 1L ||
        !is.finite(iterations) || iterations < 1 ||
        iterations != as.integer(iterations)) {
        stop("'iterations' must be one positive integer.", call. = FALSE)
    }
    .hicpotts_validate_counts(y)
    if (!identical(dim(y), c(as.integer(N), as.integer(N)))) {
        stop("'y' must have dimensions N by N.", call. = FALSE)
    }
    .hicpotts_validate_covariates(x_vars, as.integer(N))
    for (nm in c("distance", "GC", "TES", "ACC")) {
        entry <- x_vars[[nm]]
        value <- if (is.list(entry)) entry[[1L]] else entry
        if (any(value <= -1)) {
            template <- paste0(
                "x_vars$%s contains values <= -1, which are invalid ",
                "for log1p()."
            )
            stop(sprintf(template, nm), call. = FALSE)
        }
    }
    .hicpotts_validate_starts(
        gamma_prior = gamma_prior, theta_start = theta_start,
        size_start = size_start, dist = dist
    )
    if (dist %in% c("NB", "ZINB") &&
        (is.null(size_start) || length(size_start) != 3L)) {
        stop("NB/ZINB fits require three positive 'size_start' values.",
            call. = FALSE
        )
    }
    if (!is.logical(use_data_priors) || length(use_data_priors) != 1L ||
        is.na(use_data_priors)) {
        stop("'use_data_priors' must be TRUE or FALSE.", call. = FALSE)
    }
    if (!use_data_priors && is.null(user_fixed_priors)) {
        stop("Supply 'user_fixed_priors' when use_data_priors = FALSE.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

.hicpotts_fit_provenance <- function(y, x_vars, seed = NULL,
                                    initialization = NULL) {
    list(
        package = "HiCPotts",
        package_version = as.character(utils::packageVersion("HiCPotts")),
        package_library = normalizePath(
            getNamespaceInfo(asNamespace("HiCPotts"), "path"),
            winslash = "/", mustWork = TRUE
        ),
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        seed = if (is.null(seed)) NULL else as.integer(seed),
        initialization = initialization,
        rng_kind = RNGkind(),
        session_info = utils::sessionInfo(),
        data_fingerprint = hicpotts_fingerprint(y, x_vars)
    )
}

.hicpotts_attach_metadata <- function(fit, y, x_vars, seed = NULL,
                                        initialization = NULL) {
    provenance <- .hicpotts_fit_provenance(
        y, x_vars,
        seed = seed, initialization = initialization
    )
    fit$seed <- provenance$seed
    fit$initialization_method <- initialization
    fit$data_fingerprint <- provenance$data_fingerprint
    fit$component_definition <- hicpotts_component_definition()
    fit$provenance <- provenance
    fit
}

.hicpotts_execute_task <- function(task) {
    run <- function() {
        arguments <- task$arguments
        if (is.null(arguments$z_start) && !is.null(task$initialization)) {
            arguments$z_start <- make_hicpotts_initial_z(
                task$y,
                x_vars = task$x_vars, method = task$initialization,
                seed = task$seed, dist = arguments$dist,
                size_start = arguments$size_start
            )
        }
        do.call(run_metropolis_MCMC_betas, arguments)
    }
    fit <- if (is.null(task$seed)) {
        run()
    } else {
        withr::with_seed(as.integer(task$seed), run())
    }
    .hicpotts_attach_metadata(
        fit, task$y, task$x_vars,
        seed = task$seed,
        initialization = task$initialization
    )
}

.hicpotts_parallel_tasks <- function(tasks, mc_cores = 1L) {
    mc_cores <- as.integer(mc_cores)
    if (length(mc_cores) != 1L || is.na(mc_cores) || mc_cores < 1L) {
        stop("'mc_cores' must be one positive integer.", call. = FALSE)
    }
    if (length(tasks) < 2L || mc_cores == 1L) {
        return(lapply(tasks, .hicpotts_execute_task))
    }
    workers <- min(mc_cores, length(tasks))
    if (.Platform$OS.type != "windows") {
        return(parallel::mclapply(
            tasks, .hicpotts_execute_task,
            mc.cores = workers,
            mc.preschedule = FALSE, mc.set.seed = FALSE
        ))
    }

    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    library_paths <- .libPaths()
    package_library <- dirname(normalizePath(
        find.package("HiCPotts"), winslash = "/", mustWork = TRUE
    ))
    worker_packages <- parallel::clusterCall(cluster, function(paths, lib) {
        .libPaths(paths)
        if ("package:HiCPotts" %in% search()) {
            detach(
                "package:HiCPotts", unload = TRUE,
                character.only = TRUE
            )
        }
        if ("HiCPotts" %in% loadedNamespaces()) {
            unloadNamespace("HiCPotts")
        }
        namespace <- loadNamespace("HiCPotts", lib.loc = lib)
        normalizePath(
            getNamespaceInfo(namespace, "path"),
            winslash = "/", mustWork = TRUE
        )
    }, library_paths, package_library)
    expected_package <- normalizePath(
        file.path(package_library, "HiCPotts"),
        winslash = "/", mustWork = TRUE
    )
    if (!all(vapply(
        worker_packages, identical, logical(1), expected_package
    ))) {
        stop(
            "At least one PSOCK worker loaded HiCPotts from the wrong library.",
            call. = FALSE
        )
    }
    parallel::parLapply(
        cluster, tasks,
        function(task) {
            execute <- get(".hicpotts_execute_task", asNamespace("HiCPotts"))
            execute(task)
        }
    )
}

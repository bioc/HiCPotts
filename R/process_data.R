#' @title CSV file Data Processing for Hi-C Interaction Matrices and Covariates
#'
#' @description The \code{process_data} function takes a data frame with
#' interaction information and associated covariates (genomic bins, GC content,
#' transposable elements(TES), and accessibility) and converts it into a
#' structured list of \eqn{N \times N} matrices suitable for modeling. It also
#' provides options for scaling interaction counts and checking for required
#' columns and missing values.
#'
#' @usage
#' process_data(data, N, scale_max = NA_real_,
#'      standardization_y = FALSE, pad_with_zero = FALSE,
#'      mirror = FALSE)
#'
#' @param data A \code{data.frame} containing the following required columns:
#'   \itemize{
#'     \item \code{start}: Numeric start coordinates for locus i.
#'     \item \code{end}: Numeric end coordinates for locus j.
#'     \item \code{interactions}: Numeric vector of observed interaction counts
#'     between loci i and j.
#'     \item \code{GC}: GC content measure for the given loci.
#'     \item \code{ACC}: Accessibility score for the given loci.
#'     \item \code{TES}: A measure related to transposable elements for the
#'     given loci.
#'   }
#'
#' @param N An integer specifying the dimension of the resulting \eqn{N \times
#' N} matrices. The data provided should correspond to \eqn{N^2} interactions or
#' multiples thereof, as it will be reshaped into one or more \eqn{N \times N}
#' matrices.
#'
#' @param scale_max A numeric value indicating the maximum scaling factor for
#' interaction counts when \code{standardization_y = TRUE}.
#'
#' @param standardization_y A logical value. If \code{TRUE}, interaction counts
#' are scaled by applying a min-max normalization and rounding. If \code{FALSE},
#' the raw interaction counts are used as provided.
#'
#' @param pad_with_zero Logical; if TRUE, zero-pads the last block when
#' nrow(data) is not a multiple of N^2.
#'
#' @param mirror Logical; if \code{TRUE}, treat \code{data} as an upper- or
#' lower-triangle half matrix and fill the opposite triangle by reflection
#' before reshaping. Hi-C contacts are symmetric, so cell \eqn{(i,j)} and
#' cell \eqn{(j,i)} are two records of one measured contact; reflecting
#' therefore adds no fabricated observations, unlike
#' \code{pad_with_zero}. Rows are also returned in the column-major order
#' the rest of the package assumes. A half matrix over \eqn{N} bins holds
#' \eqn{N(N+1)/2} rows with the diagonal or \eqn{N(N-1)/2} without it,
#' neither of which is a multiple of \eqn{N^2}, which is why such input
#' otherwise fails the block check.
#'
#' @details This function processes a long-format data frame where each row
#' represents an interaction between two loci (denoted by \code{start} and
#' \code{end}) and associated covariate values. Key steps include:
#' \enumerate{
#'   \item Validating that required columns are present and checking for missing
#'   values.
#'   \item Optionally scaling the interaction counts to a fixed range to reduce
#'   skew or prepare data for models sensitive to the magnitude of counts.
#'   \item Computing the genomic distance between the two bins of each pair.
#'   When the \code{start.j.} column is present (as in \code{get_data()} output)
#'   the distance is \eqn{|start.j. - start|}, the separation between the two
#'   bin STARTS, which is zero on the diagonal and equal for both mirrored
#'   copies of a contact. Designs carrying only \code{start}/\code{end} fall
#'   back to \eqn{|end - start|}; note that \code{end} is bin j's END, so that
#'   legacy form spans both bins. The definition actually used is reported in a
#'   message.
#'   \item Reshaping the vectorized interaction and covariate data into \eqn{N
#'   \times N} matrices. If the number of rows in \code{data} is greater than
#'   \eqn{N^2}, it splits them into multiple \eqn{N \times N} matrices (one for
#'   each segment of length \eqn{N^2}).
#'   \item Storing the resulting data in a list structure suitable for
#'   downstream analyses.
#' }
#'
#' The returned \code{x_vars} list contains multiple \eqn{N \times N} matrices
#' for each covariate:
#' \itemize{
#'   \item \code{x_vars$distance}: A list of one or more \eqn{N \times N}
#'   matrices containing genomic distances.
#'   \item \code{x_vars$GC}: A list of \eqn{N \times N} matrices for GC content.
#'   \item \code{x_vars$TES}: A list of \eqn{N \times N} matrices for TES
#'   values.
#'   \item \code{x_vars$ACC}: A list of \eqn{N \times N} matrices for
#'   accessibility scores.
#' }
#'
#' The \code{y} element in the returned list contains the scaled (or raw)
#' interaction counts structured as a list of \eqn{N \times N} matrices.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{x_vars}}{A named list containing lists of \eqn{N \times N}
#'   matrices for each covariate: \code{distance}, \code{GC}, \code{TES},
#'   \code{ACC}.}
#'
#'   \item{\code{y}}{A list of \eqn{N \times N} matrices representing the
#'   (scaled) interaction counts corresponding to the covariates in
#'   \code{x_vars}.}
#' }
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'     start = rep(1:10, each = 10),
#'     end = rep(11:20, times = 10),
#'     interactions = rpois(100, 5),
#'     GC = runif(100, 0, 1),
#'     TES = runif(100, 0, 1),
#'     ACC = runif(100, 0, 1)
#' )
#' processed <- process_data(df,
#'     N = 10, scale_max = NA_real_,
#'     standardization_y = FALSE, pad_with_zero = FALSE
#' )
#' # x_vars <- processed$x_vars
#' # y_matrices <- processed$y
#' # str(x_vars) # Show structure of covariates
#' # str(y_matrices) # Show structure of interaction matrices
#'
#'
#' # Extended example with larger dataset
#' # Suppose we have a data frame 'large_df' corresponding to a
#' # 20x20 interaction matrix
#' # large_df <- data.frame(
#' #  start = rep(1:20, each = 20),
#' #  end = rep(21:40, times = 20),
#' #  interactions = rpois(400, 5),
#' #  GC = runif(400, 0, 1),
#' #  TES = runif(400, 0, 1),
#' #  ACC = runif(400, 0, 1)
#' # )
#' # processed <- process_data(large_df, N = 20, scale_max = NA_real_,
#' #   standardization_y = FALSE, pad_with_zero = FALSE)
#' # x_vars <- processed[[1]]
#' # y_matrices <- processed[[2]]
#' # str(x_vars)
#' # str(y_matrices)
#' # See vignette("HiCPotts_vignette") for detailed examples with
#' # real Hi-C data.
#' #
#'
#' @seealso
#' \code{\link{run_chain_betas}} for downstream MCMC inference
#' and modeling steps once data is processed into matrix form.
#'
#' @export
#
#
process_data <- function(
    data, N, scale_max = NA_real_,
    standardization_y = FALSE, pad_with_zero = FALSE, mirror = FALSE
) {
    # Check required columns in the data
    .check_required_columns(data)

    # Check for NA values in the data
    if (anyNA(data)) {
        stop(
            "The data contains NA values. Please handle missing data before ",
            "processing."
        )
    }

    if (!is.numeric(N) || length(N) != 1L || N < 1L) {
        stop("N must be a positive integer scalar.")
    }

    if (!is.logical(mirror) || length(mirror) != 1L || is.na(mirror)) {
        stop("'mirror' must be TRUE or FALSE.", call. = FALSE)
    }
    if (isTRUE(mirror)) {
        data <- .hicpotts_mirror_half_matrix(data)
    }

    block <- N * N
    nr <- nrow(data)
    if (nr %% block != 0L && !pad_with_zero) {
        ## A triangular row count is the usual cause: the caller supplied
        ## half of a symmetric contact map. Say so, because padding such
        ## input with zeros would invent roughly as many observations as
        ## it was given.
        n_bins <- length(unique(c(data$start, data$start.j.)))
        hint <- ""
        if (n_bins * (n_bins + 1L) / 2L == nr ||
            n_bins * (n_bins - 1L) / 2L == nr) {
            hint <- sprintf(
                paste0(
                    " 'data' holds one triangle of a %d-bin map, so it ",
                    "looks like a half matrix: call process_data() with ",
                    "mirror = TRUE and N = %d to fill the other triangle ",
                    "from the measured contacts."
                ),
                n_bins, n_bins
            )
        } else if (n_bins * n_bins == nr) {
            hint <- sprintf(
                " 'data' spans %d distinct bins, so N = %d.",
                n_bins, n_bins
            )
        }
        template <- paste0(
            "nrow(data) = %d is not a multiple of N^2 = %d. Either ",
            "choose a different N, or pass pad_with_zero = TRUE to ",
            "zero-pad the final (partial) matrix (not recommended for ",
            "inference).%s"
        )
        stop(sprintf(template, nr, block, hint))
    }
    ## ---- interactions -------------------------------------------------------
    y_dat <- data$interactions

    ## HiCPotts models interaction COUNTS with Poisson/NB/ZIP/ZINB
    ## distributions, which are only defined for non-negative integers.
    ## Normalised / balanced Hi-C values (KR-, ICE-, VC-normalised, etc.) are
    ## continuous: dpois()/ dnbinom() then return -Inf for every non-integer
    ## cell (invalidating the likelihood) and emit one warning per element
    ## (making the sampler dozens of times slower). Fail fast here with an
    ## actionable message instead.
    if (!is.numeric(y_dat)) {
        stop("'interactions' must be a numeric vector of counts.")
    }
    if (any(y_dat < 0, na.rm = TRUE)) {
        stop(
            "'interactions' contains negative values; counts must be ",
            "non-negative."
        )
    }
    if (!isTRUE(standardization_y)) {
        n_noninteger <- sum(abs(y_dat - round(y_dat)) > 1e-8, na.rm = TRUE)
        if (n_noninteger > 0L) {
            template <- paste0(
                "'interactions' has %d of %d non-integer values. ",
                "HiCPotts is a count model (Poisson/NB/ZIP/ZINB) and ",
                "requires RAW integer contact counts, not ",
                "normalised/balanced Hi-C values. Supply the raw counts. ",
                "(If you deliberately want to discretise continuous ",
                "values, round them yourself first: data$interactions <- ",
                "round(data$interactions).)"
            )
            stop(sprintf(template, n_noninteger, length(y_dat)))
        }
    }

    if (isTRUE(standardization_y)) {
        if (is.na(scale_max) || !is.numeric(scale_max) || scale_max <= 1) {
            stop(
                "scale_max must be a numeric > 1 when standardization_y = TRUE."
            )
        }
        warning(
            "standardization_y = TRUE min-max rescales count data. This is ",
            "not recommended for Poisson/NB/ZIP/ZINB regression because it ",
            "breaks the count assumption. Use only for backwards compatibility."
        )
        min_val <- min(y_dat)
        max_val <- max(y_dat)
        if (isTRUE(all.equal(min_val, max_val))) {
            warning("All interactions are equal; scaling disabled.")
            scaled_data <- y_dat
        } else {
            scaled_data <- round((y_dat - min_val) / (max_val - min_val) *
                (scale_max - 1) + 1)
        }
    } else {
        scaled_data <- y_dat
    }

    ## covariates##
    # Genomic distance between the two bins of each pair (see
    # .hicpotts_genomic_distance(): prefers start.j., falls back to end/start).
    x_dist <- .hicpotts_genomic_distance(data)
    x_gc <- data$GC
    x_te <- data$TES
    x_acc <- data$ACC

    split_into_blocks <- function(v) {
        split(v, ceiling(seq_along(v) / block))
    }

    to_matrices <- function(lst) {
        lapply(lst, function(x) {
            pad <- max(block - length(x), 0L)
            matrix(c(x, rep(0, pad)), nrow = N, ncol = N)
        })
    }

    x_vars <- list(
        distance = to_matrices(split_into_blocks(x_dist)),
        GC       = to_matrices(split_into_blocks(x_gc)),
        TES      = to_matrices(split_into_blocks(x_te)),
        ACC      = to_matrices(split_into_blocks(x_acc))
    )
    y_out <- to_matrices(split_into_blocks(scaled_data))

    list(x_vars = x_vars, y = y_out)
}

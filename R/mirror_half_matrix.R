## Reflect a half matrix of Hi-C bin pairs into the full square lattice.
##
## Contacts are symmetric: (i,j) and (j,i) are two records of the same
## measurement, so reflecting recovers real observations rather than
## inventing them. Rows come back in the column-major order (bin i
## fastest) that process_data() and classify_hicpotts() assume, because
## probabilities are later attached positionally.
.hicpotts_mirror_half_matrix <- function(data) {
    bins <- sort(unique(c(data$start, data$start.j.)))
    pair_key <- function(a, b) paste(a, b, sep = "\r")
    seen <- pair_key(data$start, data$start.j.)
    if (anyDuplicated(seen)) {
        stop(
            "'data' contains duplicate bin pairs, so it is not a half ",
            "matrix.",
            call. = FALSE
        )
    }
    off_diagonal <- data$start != data$start.j.
    reversed <- pair_key(data$start.j., data$start)
    if (any(reversed[off_diagonal] %in% seen)) {
        stop(
            "'data' already holds both orientations of at least one bin ",
            "pair, so it is a full matrix rather than a half matrix. Call ",
            "process_data() without mirror = TRUE.",
            call. = FALSE
        )
    }
    swapped <- data[off_diagonal, , drop = FALSE]
    bin_i_start <- swapped$start
    bin_i_end <- swapped$end.i.
    swapped$start <- swapped$start.j.
    swapped$end.i. <- swapped$end
    swapped$start.j. <- bin_i_start
    swapped$end <- bin_i_end
    full <- rbind(data, swapped)

    grid <- expand.grid(start = bins, start.j. = bins)
    idx <- match(
        pair_key(grid$start, grid$start.j.),
        pair_key(full$start, full$start.j.)
    )
    if (anyNA(idx)) {
        template <- paste0(
            "Reflection left %d of %d lattice cells empty, so 'data' is ",
            "not a complete triangle over its %d bins."
        )
        stop(
            sprintf(template, sum(is.na(idx)), nrow(grid), length(bins)),
            call. = FALSE
        )
    }
    out <- full[idx, , drop = FALSE]
    rownames(out) <- NULL
    out
}

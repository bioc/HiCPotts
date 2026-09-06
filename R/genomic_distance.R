# no roxygen needed as it is a helper function
#'
#' Genomic distance between the two bins of a bin pair.
#'
#' \code{get_data()} returns \code{start} = bin i's START and \code{end} =
#' bin j's END, so \code{abs(end - start)} spans BOTH bins: it is one bin width
#' on the diagonal instead of zero, and it differs between the two mirrored
#' copies of the same contact. The correct separation is between the two bin
#' starts, which requires the \code{start.j.} column.
#'
#' Older inputs (notably the generative simulation designs) carry only
#' \code{start}/\code{end} and were generated under the legacy expression, so
#' they are still supported: truth and fit stay self-consistent there. Which
#' definition was used is reported so it is never ambiguous after the fact.
#'
#' @param data A data frame with \code{start}/\code{end}, optionally
#' \code{start.j.}.
#' @param quiet Logical; suppress the informational message.
#' @return Numeric vector of genomic distances, one per row.
#' @noRd
.hicpotts_genomic_distance <- function(data, quiet = FALSE) {
    has_startj <- "start.j." %in% names(data) &&
        is.numeric(data[["start.j."]]) &&
        !anyNA(data[["start.j."]])

    if (has_startj) {
        if (!isTRUE(quiet)) {
            message(
                "Genomic distance: using abs(start.j. - start)"
            )
        }
        abs(as.numeric(data[["start.j."]]) - as.numeric(data[["start"]]))
    } else {
        if (!isTRUE(quiet)) message("Genomic distance: using abs(end - start)")
        abs(as.numeric(data[["end"]]) - as.numeric(data[["start"]]))
    }
}

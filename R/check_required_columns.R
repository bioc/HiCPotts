# no roxygen needed as it is a helper function
#' @noRd
.check_required_columns <- function(data,
                                    required_cols = c(
                                        "start", "end",
                                        "interactions", "GC", "TES", "ACC"
                                    )) {
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0L) {
        stop(
            "The following required columns are missing: ",
            paste(missing_cols, collapse = ", ")
        )
    }
}

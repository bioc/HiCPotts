## R/globals.R
if (getRversion() >= "2.15.1") {
    utils::globalVariables(
        c(
            "bin1_chrom", "bin2_chrom",
            "bin1_start", "bin1_end",
            "bin2_start", "bin2_end"
        )
    )
}

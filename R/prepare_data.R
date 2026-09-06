#' @title Extract Hi-C Bin Interactions with GC, Accessibility, and TE Counts
#'
#' @description Reads a Hi-C contact matrix from a \code{.cool} or \code{.h5}
#' (HDF5) file, subsets to a specified genomic region, rebins to a desired
#' resolution (merging or slicing as needed), and returns a complete grid of
#' bin-pair interactions. For each bin-pair it computes:
#' \itemize{
#'   \item GC content from a BSgenome
#'   \item DNase-I accessibility (if provided, lifted over)
#'   \item Transposable element (TE) overlap counts
#'   \item Interaction counts
#' }
#'
#' @usage
#' get_data(
#'   file_path, chr, start, end, resolution,
#'   genome_package = NULL, acc_wig = NULL,
#'   chain_file = NULL, te_granges = NULL
#' )
#'
#' @param file_path Path to a \code{.cool} file or a supported HDF5
#' schema. Binary \code{.hic} files are not read directly; convert them to
#' \code{.cool} or extract counts with \code{strawr::straw()} first.
#' @param chr Chromosome name (e.g. \code{"chr2L"}).
#' @param start 1-based start coordinate of the region.
#' @param end End coordinate of the region.
#' @param resolution Target bin size (bp); smaller native bins are merged,
#' larger ones are sliced.
#' @param genome_package Optional name of a BSgenome package (e.g.
#' \code{"BSgenome.Dmelanogaster.UCSC.dm6"}) for GC calculations. When
#' \code{NULL}, the returned GC column is \code{NA}.
#' @param acc_wig (Optional) Path to a DNase-I wig (bedGraph) file. Requires
#' \code{chain_file} to lift over.
#' @param chain_file (Optional) LiftOver chain file for mapping \code{acc_wig}.
#' @param te_granges (Optional) Path to a BED/GTF of TE annotations or a
#' \code{GRanges} object. Only TEs on \code{chr} are counted.
#'
#' @return
#' A \code{data.frame} with one row per bin-pair across \code{[start,end]} at
#' \code{resolution}, containing:
#' \describe{
#'   \item{\code{start}}{Bin1 start.}
#'   \item{\code{end.i.}}{Bin1 end.}
#'   \item{\code{start.j.}}{Bin2 start.}
#'   \item{\code{end}}{Bin2 end.}
#'   \item{\code{chrom}}{Chromosome name.}
#'   \item{\code{GC}}{Combined GC fraction of the two bins.}
#'   \item{\code{ACC}}{Mean DNase-I score per bin or \code{NA}.}
#'   \item{\code{TES}}{Sum of TE overlaps in both bins or \code{NA}.}
#'   \item{\code{interactions}}{Observed contact count (\code{0} if absent).}
#' }
#'
#' @seealso
#' \code{\link[rhdf5]{h5read}}, \code{\link[rtracklayer]{liftOver}}
#'
#' @examples
#' ## A 500 kb window of the bundled chr4 Hi-C matrix at 100 kb resolution.
#' ## GC, accessibility and TE columns stay NA unless the optional
#' ## annotation arguments are supplied.
#' cool_file <- system.file(
#'     "extdata", "BG3_WT_merged_hic_matrix_chr4_100Kb.cool",
#'     package = "HiCPotts"
#' )
#' bins <- get_data(
#'     file_path  = cool_file,
#'     chr        = "chr4",
#'     start      = 1,
#'     end        = 500000,
#'     resolution = 100000
#' )
#' head(bins)
#'
#' @importFrom Biostrings getSeq DNAString letterFrequency
#' @importFrom GenomicRanges GRanges seqnames countOverlaps findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer import.chain import liftOver
#' @importFrom utils read.table
#' @importFrom stats aggregate
#' @importFrom rhdf5 h5ls h5read
#' @importFrom S4Vectors queryHits
#' @export
get_data <- function(file_path, chr, start, end, resolution,
                    genome_package = NULL,
                    acc_wig = NULL, chain_file = NULL,
                    te_granges = NULL) {
    if (grepl("\\.hic$", file_path, ignore.case = TRUE)) {
        stop(
            ".hic files are not HDF5 and are not supported by get_data(). ",
            "Convert to .cool with cooler/hic2cool or use strawr::straw() ",
            "directly."
        )
    }

    required_pkgs <- c("rhdf5", "rtracklayer", "GenomicRanges")
    missing_pkgs <- required_pkgs[!vapply(required_pkgs,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
    )]
    if (length(missing_pkgs)) {
        stop(
            "Missing required packages:\n  - ",
            paste(missing_pkgs, collapse = "\n  - "),
            "\nPlease install them before running get_data()."
        )
    }

    has_genome <- FALSE
    genome <- NULL
    if (!is.null(genome_package)) {
        if (requireNamespace("BSgenome", quietly = TRUE) &&
            requireNamespace("Biostrings", quietly = TRUE) &&
            requireNamespace(genome_package, quietly = TRUE)) {
            genome <- get(genome_package, envir = asNamespace(genome_package))
            has_genome <- TRUE
        } else {
            warning(
                "BSgenome/Biostrings or ", genome_package,
                " not installed; GC will be NA."
            )
        }
    }

    ## DNase-I (ACC)
    acc_gr <- NULL
    if (!is.null(acc_wig) && !is.null(chain_file)) {
        dnase_df <- utils::read.table(acc_wig,
            header = FALSE, skip = 1,
            col.names = c("chr", "start", "end", "score"),
            stringsAsFactors = FALSE
        )
        if (!grepl("^chr", dnase_df$chr[1])) {
            dnase_df$chr <- paste0("chr", dnase_df$chr)
        }
        acc_gr <- GenomicRanges::GRanges(dnase_df$chr,
            IRanges::IRanges(dnase_df$start + 1, dnase_df$end),
            score = as.numeric(dnase_df$score)
        )
        acc_gr <- unlist(rtracklayer::liftOver(
            acc_gr,
            rtracklayer::import.chain(chain_file)
        ))
        acc_gr <- acc_gr[GenomicRanges::seqnames(acc_gr) == chr]
    }

    ## TEs
    te_in <- NULL
    if (!is.null(te_granges)) {
        if (is.character(te_granges)) {
            te_in <- rtracklayer::import(te_granges)
        } else if (inherits(te_granges, "GRanges")) {
            te_in <- te_granges
        } else {
            stop("`te_granges` must be a file path or a GRanges object")
        }
        te_in <- te_in[GenomicRanges::seqnames(te_in) == chr]
    }

    ## H5 schema dispatch
    hc <- rhdf5::h5ls(file_path)
    has_bins <- any(hc$group == "/" & hc$name == "bins")
    has_intervals <- any(hc$group == "/" & hc$name == "intervals")

    ## Attach bin coordinates to a pixel list and keep only the pairs with
    ## both ends on 'ch'. 'bin_id' is a dense zero-based row index into
    ## 'bins_tbl', so both joins are positional lookups. Doing them with
    ## merge() costs seconds per call because the whole pixel table -- often
    ## millions of rows genome-wide -- is joined twice before anything is
    ## filtered. Selecting the wanted pairs first and then indexing gives the
    ## same rows and the same columns.
    join_bin_coordinates <- function(bin1_id, bin2_id, count, bins_tbl, ch) {
        n_bins <- nrow(bins_tbl)
        i1 <- bin1_id + 1
        i2 <- bin2_id + 1
        ## merge() drops pixels whose bin is absent from bins_tbl; an
        ## out-of-range index must therefore be discarded, not turned into NA.
        keep <- !is.na(i1) & !is.na(i2) &
            i1 >= 1 & i1 <= n_bins & i2 >= 1 & i2 <= n_bins
        i1 <- i1[keep]
        i2 <- i2[keep]
        chrom <- bins_tbl$chrom
        on_chr <- chrom[i1] == ch & chrom[i2] == ch
        on_chr[is.na(on_chr)] <- FALSE
        i1 <- i1[on_chr]
        i2 <- i2[on_chr]
        sel <- which(keep)[on_chr]
        out <- data.frame(
            bin2_id = bin2_id[sel],
            bin1_id = bin1_id[sel],
            count = count[sel],
            bin1_chrom = chrom[i1],
            bin1_start = bins_tbl$start[i1],
            bin1_end = bins_tbl$end[i1],
            bin2_chrom = chrom[i2],
            bin2_start = bins_tbl$start[i2],
            bin2_end = bins_tbl$end[i2],
            stringsAsFactors = FALSE
        )
        row.names(out) <- NULL
        out
    }

    read_cool <- function(fp, ch) {
        codes_raw <- rhdf5::h5read(fp, "bins/chrom")
        names_map <- rhdf5::h5read(fp, "chroms/name")
        add_chr_prefix <- function(x) {
            x <- as.character(x)
            ifelse(grepl("^chr", x), x, paste0("chr", x))
        }

        if (is.factor(codes_raw)) {
            cl <- as.character(codes_raw)
            bin_chrom <- add_chr_prefix(cl)
        } else {
            bin_chrom <- add_chr_prefix(names_map[codes_raw + 1])
        }
        bin_start <- rhdf5::h5read(fp, "bins/start") + 1
        bin_end <- rhdf5::h5read(fp, "bins/end")
        bins_tbl <- data.frame(
            bin_id = seq_along(bin_chrom) - 1,
            chrom = bin_chrom, start = bin_start, end = bin_end,
            stringsAsFactors = FALSE
        )
        join_bin_coordinates(
            rhdf5::h5read(fp, "pixels/bin1_id"),
            rhdf5::h5read(fp, "pixels/bin2_id"),
            rhdf5::h5read(fp, "pixels/count"),
            bins_tbl, ch
        )
    }

    read_interval <- function(fp, ch) {
        ints <- rhdf5::h5read(fp, "intervals")
        ints$start_list <- ints$start_list + 1
        bins_tbl <- data.frame(
            bin_id = seq_along(ints$start_list) - 1,
            chrom = ints$chr_list,
            start = ints$start_list, end = ints$end_list,
            stringsAsFactors = FALSE
        )
        dat <- rhdf5::h5read(fp, "matrix")
        rows <- rep(seq_len(dat$shape[1]) - 1, diff(dat$indptr))
        join_bin_coordinates(rows, dat$indices, dat$data, bins_tbl, ch)
    }

    raw_df <- if (has_bins) {
        read_cool(file_path, chr)
    } else if (has_intervals) {
        read_interval(file_path, chr)
    } else {
        stop("Unrecognised .h5 schema: neither /bins nor /intervals found")
    }

    raw_df <- subset(
        raw_df,
        bin1_start <= end & bin1_end >= start &
            bin2_start <= end & bin2_end >= start
    )

    ## If the cooler file stores only the upper triangle, mirror off-diagonal
    ## entries so the downstream full grid is not filled with artificial
    ## lower-triangle zeros.
    if (nrow(raw_df) > 0 &&
        all(raw_df$bin1_start <= raw_df$bin2_start, na.rm = TRUE)) {
        offdiag <- raw_df$bin1_start != raw_df$bin2_start |
            raw_df$bin1_end != raw_df$bin2_end

        if (any(offdiag)) {
            mirrored <- raw_df[offdiag, , drop = FALSE]

            tmp_chrom <- mirrored$bin1_chrom
            tmp_start <- mirrored$bin1_start
            tmp_end <- mirrored$bin1_end

            mirrored$bin1_chrom <- mirrored$bin2_chrom
            mirrored$bin1_start <- mirrored$bin2_start
            mirrored$bin1_end <- mirrored$bin2_end

            mirrored$bin2_chrom <- tmp_chrom
            mirrored$bin2_start <- tmp_start
            mirrored$bin2_end <- tmp_end

            raw_df <- rbind(raw_df, mirrored)
        }
    }

    ## Native resolution & re-binning
    if (nrow(raw_df) > 0) {
        native_res <- unique(raw_df$bin1_end - raw_df$bin1_start + 1)
        native_res <- if (length(native_res) > 1) {
            min(
                native_res
            )
        } else {
            native_res
        }
        obs_df <- raw_df
    } else {
        native_res <- resolution
        obs_df <- raw_df[FALSE, ]
    }

    if (native_res != resolution && nrow(obs_df) > 0) {
        if (native_res < resolution) {
            ## finer native -> coarser target: aggregate
            obs_df$new_s1 <- floor((
                obs_df$bin1_start - 1) / resolution) * resolution + 1
            obs_df$new_e1 <- obs_df$new_s1 + resolution - 1
            obs_df$new_s2 <- floor((
                obs_df$bin2_start - 1) / resolution) * resolution + 1
            obs_df$new_e2 <- obs_df$new_s2 + resolution - 1
            agg <- stats::aggregate(count ~ new_s1 + new_e1 + new_s2 + new_e2,
                data = obs_df, FUN = sum
            )
            obs_df <- data.frame(
                bin1_start = agg$new_s1, bin1_end = agg$new_e1,
                bin2_start = agg$new_s2, bin2_end = agg$new_e2,
                count = agg$count, stringsAsFactors = FALSE
            )
        } else {
            ## coarser native -> finer target: slice with weights
            if (native_res != resolution && nrow(obs_df) > 0) {
                if (native_res < resolution) {
                    ## finer native -> coarser target: aggregate
                    obs_df$new_s1 <- floor((
                        obs_df$bin1_start - 1) / resolution) * resolution + 1
                    obs_df$new_e1 <- obs_df$new_s1 + resolution - 1
                    obs_df$new_s2 <- floor((
                        obs_df$bin2_start - 1) / resolution) * resolution + 1
                    obs_df$new_e2 <- obs_df$new_s2 + resolution - 1

                    agg <- stats::aggregate(
                        count ~ new_s1 + new_e1 + new_s2 + new_e2,
                        data = obs_df, FUN = sum
                    )

                    obs_df <- data.frame(
                        bin1_start = agg$new_s1,
                        bin1_end = agg$new_e1,
                        bin2_start = agg$new_s2,
                        bin2_end = agg$new_e2,
                        count = agg$count,
                        stringsAsFactors = FALSE
                    )
                } else {
                    stop(
                        "Requested resolution (", resolution,
                        ") is finer than the native Hi-C matrix resolution (",
                        native_res,
                        "). Upsampling coarse Hi-C bins into finer bins is ",
                        "not statistically meaningful. ",
                        "Please use resolution >= native_res."
                    )
                }
            }
        }
    }

    ## ---- ONE shared bin grid, used by both the observed and full-grid bins
    ## ---- Bins are anchored to chromosome coordinate 1: the bin covering
    ## coordinate p starts at floor((p - 1) / resolution) * resolution + 1.
    ## read_cool() and read_interval() convert the cooler's 0-based starts to
    ## this same 1-based convention, and the rebinning above uses the identical
    ## rule, so observed and generated bins now share one origin. Previously the
    ## observed side was anchored at coordinate 1 while the output grid was
    ## anchored at the user-supplied `start`, so an arbitrary genomic window
    ## could produce two incompatible grids and return only zeros.
    ##
    ## The requested region is snapped OUTWARD to whole bins, so a request that
    ## does not land on bin boundaries returns the bins covering it rather than
    ## a shifted grid or a truncated partial bin.
    bin_of <- function(pos) floor((as.numeric(pos) - 1) / resolution)
    start_of <- function(idx) idx * resolution + 1

    first_bin <- bin_of(start)
    last_bin <- bin_of(end)
    bin_idx <- seq.int(first_bin, last_bin)
    n_bins_g <- length(bin_idx)
    bs <- start_of(bin_idx)

    ## Keep observed bins lying inside the snapped target region.
    if (nrow(obs_df) > 0) {
        ob1 <- bin_of(obs_df$bin1_start)
        ob2 <- bin_of(obs_df$bin2_start)
        obs_df <- obs_df[ob1 >= first_bin & ob1 <= last_bin &
            ob2 >= first_bin & ob2 <= last_bin, , drop = FALSE]
    }

    full <- expand.grid(
        bin1_start = bs, bin2_start = bs,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    full$bin1_end <- full$bin1_start + resolution - 1
    full$bin2_end <- full$bin2_start + resolution - 1
    full$interactions <- 0
    if (nrow(obs_df) > 0) {
        ## Join on NUMERIC bin indices, never on default-formatted strings.
        ## Observed coordinates arrive from HDF5 as integers while the generated
        ## grid is double, so paste() keys like "100000" and "1e+05" compare
        ## unequal even though the coordinates are identical: every match() was
        ## NA and the pre-filled zeros survived, silently destroying all counts.
        ## The failure depended on options(scipen), so results moved with a
        ## global print option. Bin indices are small integers held exactly in
        ## double.
        key_full <- (bin_of(full$bin1_start) - first_bin) * n_bins_g +
            (bin_of(full$bin2_start) - first_bin)
        key_obs <- (bin_of(obs_df$bin1_start) - first_bin) * n_bins_g +
            (bin_of(obs_df$bin2_start) - first_bin)
        idx <- match(key_obs, key_full)
        keep <- !is.na(idx)
        full$interactions[idx[keep]] <- obs_df$count[keep]
    }

    ## GC content
    if (has_genome) {
        full$GC <- mapply(function(s1, e1, s2, e2) {
            seq1 <- Biostrings::getSeq(genome, chr, s1, e1)
            seq2 <- Biostrings::getSeq(genome, chr, s2, e2)
            dna <- Biostrings::DNAString(paste0(
                as.character(seq1),
                as.character(seq2)
            ))
            ## letterFrequency with c("G","C") and default OR="|" returns GC
            ## fraction.
            sum(Biostrings::letterFrequency(dna, c("G", "C"), as.prob = TRUE))
        }, full$bin1_start, full$bin1_end, full$bin2_start, full$bin2_end)
    } else {
        full$GC <- NA_real_
    }

    ## ACC
    full$ACC <- if (is.null(acc_gr)) {
        NA_real_
    } else {
        ## Accessibility is averaged over the bedGraph records overlapping a
        ## bin, WEIGHTED by the width of each record's overlap with that bin. An
        ## unweighted mean lets a record covering a handful of bases count as
        ## much as one spanning the whole bin, which biases the score whenever
        ## record widths are unequal (they usually are).
        bin_score <- function(s, e) {
            g <- GenomicRanges::GRanges(chr, IRanges::IRanges(s, e))
            hits <- GenomicRanges::findOverlaps(acc_gr, g)
            if (length(hits) == 0L) {
                return(NA_real_)
            }
            recs <- acc_gr[S4Vectors::queryHits(hits)]
            ov <- pmin(GenomicRanges::end(recs), e) - pmax(GenomicRanges::start(
                recs
            ), s) + 1
            ov <- pmax(as.numeric(ov), 0)
            sc <- recs$score
            ok <- !is.na(sc) & ov > 0
            if (!any(ok)) {
                return(NA_real_)
            }
            sum(sc[ok] * ov[ok]) / sum(ov[ok])
        }
        mapply(function(s1, e1, s2, e2) {
            m1 <- bin_score(s1, e1)
            m2 <- bin_score(s2, e2)
            acc_val <- mean(c(m1, m2), na.rm = TRUE)
            if (is.nan(acc_val)) NA_real_ else acc_val
        }, full$bin1_start, full$bin1_end, full$bin2_start, full$bin2_end)
    }

    ## TES
    full$TES <- if (is.null(te_in)) {
        NA_integer_
    } else {
        mapply(function(s1, e1, s2, e2) {
            g1 <- GenomicRanges::GRanges(chr, IRanges::IRanges(s1, e1))
            g2 <- GenomicRanges::GRanges(chr, IRanges::IRanges(s2, e2))
            sum(GenomicRanges::countOverlaps(g1, te_in)) +
                sum(GenomicRanges::countOverlaps(g2, te_in))
        }, full$bin1_start, full$bin1_end, full$bin2_start, full$bin2_end)
    }

    data.frame(
        start = full$bin1_start,
        `end.i.` = full$bin1_end,
        `start.j.` = full$bin2_start,
        end = full$bin2_end,
        chrom = chr,
        GC = full$GC,
        ACC = full$ACC,
        TES = full$TES,
        interactions = as.double(full$interactions), # FIX: no destructive round
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

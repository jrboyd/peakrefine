
#' remove positional read duplicates from data.table
#'
#' @param reads_dt a data.table containing read span information.  required
#' columns are which_label, start, end, and strand.
#' @param max_dupes integer >= 1, the maximum allowed duplicates per position.
#'
#' @return data.table with reads in excess of max_dupes removed
#'
#' @examples
#' library(data.table)
#' n = 4
#' dt = data.table(
#'   which_label = seq_len(n),
#'   seqnames = "chr1",
#'   strand = c("+", "-"),
#'   start = seq_len(n),
#'   end = seq_len(n) + 10
#' )
#' dt[, width := end - start + 1]
#' make_dupes = seq_len(nrow(dt))
#' make_dupes = rep(make_dupes, make_dupes)
#' dt = dt[make_dupes]
#' lapply(1: n, function(x)peakrefine:::.remove_duplicates(dt, max_dupes = x))
.remove_duplicates = function(reads_dt, max_dupes){
    stopifnot(max_dupes >= 1)
    stopifnot(class(reads_dt) == "data.table")
    stopifnot(all(c("which_label", "start", "end", "strand") %in%
                      colnames(reads_dt)))
    # browser()
    ndupe = which_label = NULL
    reads_dt[, ndupe := 1L]
    reads_dt[strand == "+",
             ndupe := seq_len(.N)[order(width, decreasing = TRUE)],
             by = list(which_label, start)]
    reads_dt[strand == "-",
             dupe := seq_len(.N)[order(width, decreasing = TRUE)],
             by = list(which_label, end)]
    reads_dt = reads_dt[ndupe <= max_dupes]
    reads_dt$ndupe = NULL
    reads_dt
}

#' .fetch_bam_stranded
#'
#' @param bam_f character. a .bam file, must have index file in same directory:
#'   .bam.bai
#' @param qgr GRanges, regions to retrieve reads for.
#' @param max_dupes numeric, maximum positional duplicates allowed.  See
#'   \code{\link{.remove_duplicates}}
#' @param ... additional arguments passed to
#'   \code{\link[Rsamtools]{ScanBamParam}}
#'
#' @return a tidy data.table of read pileup coverage from bam_f for regions in
#'   qgr
#' @export
#' @rawNamespace import(data.table, except = c(shift, first, last, second))
#' @import Rsamtools
#' @examples
#' bam_file = system.file("extdata",
#'   "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' bam_input = system.file("extdata",
#'   "MCF10A_input.random5.bam", package = "peakrefine")
#' np = system.file("extdata",
#'   "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#' .fetch_bam_stranded(bam_file, qgr)
#'
.fetch_bam_stranded = function(bam_f,
                               qgr,
                               max_dupes = Inf,
                               ...){
    stopifnot(is.numeric(max_dupes))
    stopifnot(max_dupes >= 1)
    sbgr = qgr
    strand(sbgr) = "*"
    sbParam = Rsamtools::ScanBamParam(
        which = sbgr,
        what = c("rname", "strand", "pos", "qwidth", "cigar"),
        ...
    )
    bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
    bam_dt = lapply(bam_raw, function(x){
        data.table::data.table(seqnames = x$rname,
                               strand = x$strand,
                               start = x$pos,
                               width = x$qwidth,
                               cigar = x$cigar)
    })
    bam_dt = data.table::rbindlist(bam_dt,
                                   use.names = TRUE,
                                   idcol = "which_label")

    bam_dt[, end := start + width - 1L]

    if(max_dupes < Inf){
        bam_dt = .remove_duplicates(bam_dt, max_dupes)
    }
    bam_dt
}

#' Title
#'
#' @param reads_dt
#' @param fragLen
#'
#' @return
#' @export
#'
#' @examples
strandsExtend = function(reads_dt, fragLen){
    reads_dt[strand == "+", end := start + as.integer(fragLen) - 1L]
    reads_dt[strand == "-", start := end - as.integer(fragLen) + 1L]
    return(reads_dt)
}

#' Title
#'
#' @param score_dt
#' @param window_size
#' @param anchor
#'
#' @return
#'
#' @examples
.shift_anchor = function(score_dt, window_size, anchor){
    x = id = NULL
    shift = round(window_size/2)
    switch(anchor,
           center = {
               score_dt[, x := start - min(start) + shift, by = id]
               score_dt[, x := x - round(mean(x)), by = id]
               score_dt[strand == "-", x := -1*x]
           },
           center_unstranded = {
               score_dt[, x := start - min(start) + shift, by = id]
               score_dt[, x := x - round(mean(x)), by = id]
           },
           left = {
               score_dt[, x := -1]
               score_dt[strand != "-", x := start - min(start) + shift,
                        by = id]
               #flip negative
               score_dt[strand == "-", x := -1*(end - max(end) - shift),
                        by = id]
           },
           left_unstranded = {
               score_dt[, x := start - min(start) + shift, by = id]
           }
    )

    score_dt[, start := start - shift + 1]
    score_dt[, end := end + window_size - shift]
    return(score_dt)
}

#' Title
#'
#' @param score_gr
#' @param qgr
#' @param window_size
#' @param anchor
#'
#' @return
#' @export
#'
#' @examples
viewGRangesWinSample_dt = function(score_gr, qgr, window_size,
                                   anchor = c("center", "center_unstranded",
                                              "left", "left_unstranded")[1]){
    #reserve bindings for data.table
    x = id = queryHits = NULL
    stopifnot(class(score_gr) == "GRanges")
    stopifnot(!is.null(score_gr$score))
    stopifnot(class(qgr) == "GRanges")
    stopifnot(is.numeric(window_size))
    stopifnot(window_size >= 1)
    stopifnot(window_size %% 1 == 0)
    stopifnot(anchor %in% c("center", "center_unstranded",
                            "left", "left_unstranded"))
    if (is.null(qgr$id)) {
        if (is.null(names(qgr))) {
            names(qgr) = paste0("region_", seq_along(qgr))
        }
        qgr$id = names(qgr)
    }
    windows = slidingWindows(qgr, width = window_size, step = window_size)

    # names(windows) = qgr$id
    # windows's handling of names seems to have changed and now every nest
    # GRanges has parent's name
    names(windows) = NULL
    windows = unlist(windows)
    windows$id = names(windows)
    windows = resize(windows, width = 1, fix = "center")
    olaps = suppressWarnings(data.table::as.data.table(
        findOverlaps(query = windows,
                     subject = score_gr,
                     ignore.strand = TRUE)
    ))
    # patch up missing/out of bound data with 0
    missing_idx = setdiff(seq_along(windows), olaps$queryHits)
    if (length(missing_idx) > 0) {
        olaps = rbind(olaps, data.table::data.table(
            queryHits = missing_idx,
            subjectHits = length(score_gr) + 1))[order(queryHits)]
        suppressWarnings({
            score_gr = c(score_gr,
                         GRanges("chrPatchZero",
                                 IRanges::IRanges(1, 1), score = 0))
        })
    }
    # set y and output windows = windows[olaps$queryHits]
    windows$y = score_gr[olaps$subjectHits]$score
    score_dt = data.table::as.data.table(windows)

    return(.shift_anchor(score_dt, window_size, anchor))
}

#' Title
#'
#' @param bam_dt
#' @param qgr
#' @param win_size
#'
#' @return
#' @export
#' @import GenomicRanges
#' @examples
strandsCoverage = function(bam_dt, qgr, win_size = 1){
    ext_cov = GenomicRanges::coverage(GenomicRanges::split(GenomicRanges::GRanges(bam_dt[strand == "+"]), bam_dt$which_label))
    score_gr = GenomicRanges::GRanges(ext_cov)
    pos_dt = viewGRangesWinSample_dt(score_gr, qgr, win_size)
    pos_dt$strand = "+"

    ext_cov = GenomicRanges::coverage(GenomicRanges::split(GenomicRanges::GRanges(bam_dt[strand == "-"]), bam_dt$which_label))
    score_gr = GenomicRanges::GRanges(ext_cov)
    neg_dt = viewGRangesWinSample_dt(score_gr, qgr, win_size)
    neg_dt$strand = "-"

    rbind(pos_dt, neg_dt)
}

#' Title
#'
#' @param reads_dt
#' @param qgr
#' @param fragLen
#' @param win_size
#'
#' @return
#' @export
#' @rawNamespace import(data.table, except = c(shift, first, last, second))
#' @importFrom stats cor
#' @examples
calcStrandCorr = function(reads_dt, qgr, fragLen = NA, win_size = 1){
    if(!is.na(fragLen)){
        reads_dt = strandsExtend(reads_dt, fragLen)
    }
    pile_dt = strandsCoverage(reads_dt, qgr, win_size)
    # ggplot(pile_dt) + geom_path(aes(x = x, y = y, color = strand)) + facet_wrap("id")
    dc_dt = data.table::dcast(pile_dt, id + x ~ strand, value.var = "y")
    dc_dt = dc_dt[, .(corr = cor(`+`, `-`)) , by = .(id)]
    dc_dt
}





#
# tgr = GRanges(c("chr1", "chr6"), IRanges(c(1e6, 2e6), c(10e6, 16e6)), qValue = 1:2, id = paste("reg", 1:2))
# names(tgr) = tgr$id
# bigcorr_res = crossCorrByExtension(bam_file, tgr, frag_min = 1, frag_max = 600)
# bigcorr_res
#
# qgr = qgr[1:1600]
#
# st_times = list()
# st_res = list()
#
# options("mc.cores" = 16)
#
# # for(nc in c(1, 2, 4, 8, 16)){
# nc = 16
# print(nc)
# nam = paste("nc", nc)
# st_times[[nam]] = system.time({
#     st_res[[nam]] = crossCorrByExtensionFull(bam_file, qgr, fragLen = corr_res$frag_length, ncores = nc)
# })
# # }

#' Title
#'
#' @param bam_file
#' @param qgr
#' @param fl
#' @param rl
#' @param roi
#'
#' @return
#' @export
#'
#' @examples
eval_fragvsread = function(bam_file, qgr, fl, rl, roi){
    bdt = .fetch_bam_stranded(bam_file, qgr[roi], max_dupes = 1)
    p1 = ggplot(bdt, aes(x = x, y = y, color = strand)) + geom_path() +
        facet_wrap("id")

    bdt = .fetch_bam_stranded(bam_file, qgr[roi], max_dupes = 1)
    strandsExtend(bdt, fl)
    p2 = ggplot(bdt, aes(x = x, y = y, color = strand)) + geom_path() +
        facet_wrap("id")

    return(list(
        read_ext = p1 + labs(title = paste0("reads (", rl, ")")),
        fragment_ext = p2 + labs(title = paste0("fragments (", fl, ")"))
    ))
}

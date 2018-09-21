
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
    stopifnot("data.table" %in% class(reads_dt))
    stopifnot(all(c("which_label", "start", "end", "strand") %in%
                      colnames(reads_dt)))
    # browser()
    ndupe = which_label = NULL
    reads_dt[, ndupe := 1L]
    reads_dt[strand == "+",
             ndupe := seq_len(.N)[order(width, decreasing = TRUE)],
             by = list(which_label, start)]
    reads_dt[strand == "-",
             ndupe := seq_len(.N)[order(width, decreasing = TRUE)],
             by = list(which_label, end)]
    reads_dt = reads_dt[ndupe <= max_dupes]
    reads_dt$ndupe = NULL
    reads_dt
}

#' .fetch_bam_stranded
#'
#' @param bam_f character. a .bam file, must have index file in same directory:
#'   .bam.bai
#' @param query_gr GRanges, regions to retrieve reads for.
#' @param max_dupes numeric, maximum positional duplicates allowed.  See
#'   \code{\link{.remove_duplicates}}
#' @param ... additional arguments passed to
#'   \code{\link[Rsamtools]{ScanBamParam}}
#'
#' @return a tidy data.table of read pileup coverage from bam_f for regions in
#'   query_gr
#' @rawNamespace import(data.table, except = c(shift, first, last, second))
#' @import Rsamtools
#' @examples
#' bam_file = system.file("extdata",
#'   "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' bam_input = system.file("extdata",
#'   "MCF10A_input.random5.bam", package = "peakrefine")
#' np = system.file("extdata",
#'   "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' query_gr = rtracklayer::import(np, format = "narrowPeak")
#' peakrefine:::.fetch_bam_stranded(bam_file, query_gr)
#'
.fetch_bam_stranded = function(bam_f,
                               query_gr,
                               max_dupes = Inf,
                               ...){
    stopifnot(is.numeric(max_dupes))
    stopifnot(max_dupes >= 1)
    sbgr = query_gr
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

#' extend reads to specified length in a strand sensitve way
#'
#' @param reads_dt data.table. Contains at least 3 mandatory columns: strand, start,
#'   and end.
#' @param  frag_len integer. the fragment length to extend to
#'
#' @return data.table after extending to fragment length
#'
#' @examples
#' bam_file = system.file("extdata",
#'   "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' bam_input = system.file("extdata",
#'   "MCF10A_input.random5.bam", package = "peakrefine")
#' np = system.file("extdata",
#'   "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' query_gr = rtracklayer::import(np, format = "narrowPeak")
#' reads_dt = peakrefine:::.fetch_bam_stranded(bam_file, query_gr)
#' peakrefine:::.extend_reads(reads_dt, 150)
.extend_reads = function(reads_dt,  frag_len){
    reads_dt[strand == "+", end := start + as.integer( frag_len) - 1L]
    reads_dt[strand == "-", start := end - as.integer( frag_len) + 1L]
    return(reads_dt)
}

#' applies an anchor which controls x coordinate system 1) should 0 on the
#' x-axis be on the left or in the center? 2) should - strand signal be flipped?
#'
#' @param score_dt data.table, mandatory columns are start, end, and id.  strand
#'   is also mandatory if not anchoring *_unstranded.
#' @param window_size the window size at which score was calculated/retrieved
#' @param anchor character, one of center, center_unstranded, left,
#'   left_unstranded.
#'
#' @return data.table with x added according to anchor
#'
#' @examples
#' library(GenomicRanges)
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' bam_input = system.file("extdata", "MCF10A_input.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#' strand(qgr) = c("+", "-", "-", "+", "-")
#'
#' reads_dt = peakrefine:::.fetch_bam_stranded(bam_file, qgr)
#' ext_cov = GenomicRanges::coverage(GenomicRanges::split(
#'     GenomicRanges::GRanges(reads_dt[strand == "+"]), reads_dt$which_label))
#' score_gr = GenomicRanges::GRanges(ext_cov)
#' bam_dt = peakrefine::.view_gr_by_window_sample(score_gr, qgr, 1,
#'                                  anchor = "center_unstranded")
#' peakrefine:::.shift_anchor(bam_dt, 1, "center_unstranded")
#'
.shift_anchor = function(score_dt, window_size, anchor){
    stopifnot(length(anchor) == 1)
    stopifnot(class(anchor) == "character")
    stopifnot(anchor %in% c("center", "center_unstranded",
                            "left", "left_unstranded"))
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

    score_dt[, start := start - shift ]
    score_dt[, end := end + window_size - shift - 1]
    return(score_dt)
}

#' view GRanges with score at specified window size by sampling values
#'
#' @param score_gr GRanges with score mcol
#' @param query_gr GRanges that was used to fecth score_gr
#' @param window_size ineger. desired window size.
#' @param anchor character, anchoring method. one of "center", "center_unstranded",
#' "left", "left_unstranded"
#'
#' @return data.table of scores stored in y value sample every window_size
#' @import GenomicRanges rtracklayer
#' @examples
#' library(GenomicRanges)
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' bam_input = system.file("extdata", "MCF10A_input.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#' strand(qgr) = c("+", "-", "-", "+", "-")
#'
#' reads_dt = peakrefine:::.fetch_bam_stranded(bam_file, qgr)
#' ext_cov = coverage(split(
#'     GRanges(reads_dt[strand == "+"]), reads_dt$which_label))
#' score_gr = GRanges(ext_cov)
#' peakrefine::.view_gr_by_window_sample(score_gr, qgr, 1,
#'                                  anchor = "center_unstranded")
.view_gr_by_window_sample = function(score_gr, query_gr, window_size,
                                   anchor = c("center", "center_unstranded",
                                              "left", "left_unstranded")[1]){
    #reserve bindings for data.table
    x = id = queryHits = NULL
    stopifnot(class(score_gr) == "GRanges")
    stopifnot(!is.null(score_gr$score))
    stopifnot(class(query_gr) == "GRanges")
    stopifnot(is.numeric(window_size))
    stopifnot(window_size >= 1)
    stopifnot(window_size %% 1 == 0)
    stopifnot(anchor %in% c("center", "center_unstranded",
                            "left", "left_unstranded"))
    if (is.null(query_gr$id)) {
        if (is.null(names(query_gr))) {
            names(query_gr) = paste0("region_", seq_along(query_gr))
        }
        query_gr$id = names(query_gr)
    }
    windows = GenomicRanges::slidingWindows(query_gr, width = window_size, step = window_size)

    # names(windows) = query_gr$id
    # windows's handling of names seems to have changed and now every nest
    # GRanges has parent's name
    names(windows) = NULL
    windows = unlist(windows)
    windows$id = names(windows)
    windows = resize(windows, width = 1, fix = "center")
    olaps = suppressWarnings(data.table::as.data.table(
        GenomicRanges::findOverlaps(query = windows,
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
    # apply shift and return
    .shift_anchor(score_dt, window_size, anchor)
}

#' Calculate coverage per strand
#'
#' @param reads_dt data.table containing read spans. mandatory columns are
#'   which_label, strand, start, end
#' @param query_gr GRanges of regions of interest
#' @param window_size resolution used to represent coverage.  window_size of 10
#'   retrieves coverage every 10 bp.
#'
#' @return data.table containing stranded read coverage
#' @import GenomicRanges
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#' bam_dt = .fetch_bam_stranded(bam_file, qgr)
#' cov_dt = .calc_stranded_coverage(bam_dt, qgr)
#'
#'
.calc_stranded_coverage = function(reads_dt, query_gr, window_size = 1){
    ext_cov = GenomicRanges::coverage(GenomicRanges::split(
        GenomicRanges::GRanges(reads_dt[strand == "+"]), reads_dt$which_label))
    score_gr = GenomicRanges::GRanges(ext_cov)
    pos_dt = .view_gr_by_window_sample(score_gr, query_gr, window_size,
                                     anchor = "center_unstranded")
    pos_dt$strand = "+"

    ext_cov = GenomicRanges::coverage(GenomicRanges::split(
        GenomicRanges::GRanges(reads_dt[strand == "-"]), reads_dt$which_label))
    score_gr = GenomicRanges::GRanges(ext_cov)
    neg_dt = .view_gr_by_window_sample(score_gr, query_gr, window_size,
                                     anchor = "center_unstranded")
    neg_dt$strand = "-"

    rbind(pos_dt, neg_dt)
}

#' Calculate correlation between signal on DNA strands
#'
#' @param reads_dt data.table containing read spans. mandatory columns are
#'   which_label, strand, start, end
#' @param query_gr GRanges of regions of interest
#' @param  frag_len fragment size to extend reads to.  frag_len of NA causes read
#' extension to be skipped, ie. raw read pileup.
#' @param window_size integer, resolution used to represent coverage.  window_size of 10
#'   retrieves coverage every 10 bp.
#'
#' @return data.table of correlation per region in query_gr
#' @export
#' @rawNamespace import(data.table, except = c(shift, first, last, second))
#' @importFrom stats cor
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#' bam_dt = .fetch_bam_stranded(bam_file, qgr)
#' calcStrandCorr(bam_dt, qgr) #read correlation
#' calcStrandCorr(bam_dt, qgr, frag_len = 150) #fragment 150 bp correlation
calcStrandCorr = function(reads_dt, query_gr,  frag_len = NA, window_size = 1){
    if(!is.na( frag_len)){
        reads_dt = .extend_reads(reads_dt,  frag_len)
    }
    pile_dt = .calc_stranded_coverage(reads_dt, query_gr, window_size)
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
# query_gr = query_gr[1:1600]
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
#     st_res[[nam]] = crossCorrByExtensionFull(bam_file, query_gr,  frag_len = corr_res$frag_length, ncores = nc)
# })
# # }

#' #' Title
#' #'
#' #' @param bam_file
#' #' @param query_gr
#' #' @param fl
#' #' @param rl
#' #' @param roi
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' eval_fragvsread = function(bam_file, query_gr, fl, rl, roi){
#'     bdt = .fetch_bam_stranded(bam_file, query_gr[roi], max_dupes = 1)
#'     p1 = ggplot(bdt, aes(x = x, y = y, color = strand)) + geom_path() +
#'         facet_wrap("id")
#'
#'     bdt = .fetch_bam_stranded(bam_file, query_gr[roi], max_dupes = 1)
#'     .extend_reads(bdt, fl)
#'     p2 = ggplot(bdt, aes(x = x, y = y, color = strand)) + geom_path() +
#'         facet_wrap("id")
#'
#'     return(list(
#'         read_ext = p1 + labs(title = paste0("reads (", rl, ")")),
#'         fragment_ext = p2 + labs(title = paste0("fragments (", fl, ")"))
#'     ))
#' }

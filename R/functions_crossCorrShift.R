#' Calculate cross strand correlation using the read shift method.
#'
#' @param bam_file character. Path to .bam file, must have index at .bam.bai.
#' @param query_gr GRanges.  Regions to calculate cross correlation for.
#' @param n_regions integer.  query_gr will be downsampled to this many regions
#'   for speed. Use NA to skip downsampling.
#' @param max_dupes integer.  Duplicate reads above this value will be removed.
#' @param frag_min integer.  extension value to start at.
#' @param frag_max integer. extension value to end at.
#' @param step integer.  proceed from frag_min measuring correlation every step.
#' @param small_step integer.  after measuring correlation every step, a second
#'   round of fragment size refinement is done using small_step within +/- step
#'   of maximum.
#' @param include_plots logical. Should plots be included in output?
#'
#' @return data.table of metrics
#' @export
#' @import pbapply
#' @import ggplot2
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#' crossCorrByShift(bam_file, qgr[1:2], frag_min = 50,
#' frag_max = 250, step = 50, small_step = 10)
crossCorrByShift = function(bam_file,
                            query_gr,
                            n_regions = 20,
                            max_dupes = 1,
                            frag_min = 0, frag_max = 250,
                            step = 10, small_step = 1, include_plots = TRUE){
    which_label = x = y = shiftn = NULL #reserve for data.table
    stopifnot(is.numeric(n_regions))
    stopifnot(n_regions >= 1)
    if(is.na(n_regions) || n_regions >= length(query_gr)){
        test_gr = query_gr
    }else{
        test_gr = sample(query_gr, n_regions)
    }
    if(is.null(test_gr$id)){
        test_gr$id = paste0("peak_", seq_along(test_gr))
    }
    if(is.null(names(test_gr))){
        names(test_gr) = test_gr$id
    }

    message("fetch reads...")
    # browser()
    reads_dt = .fetch_bam_stranded(bam_file, test_gr, max_dupes = max_dupes)

    tab = table(reads_dt$width)
    read_length = as.numeric(names(tab[which(tab == max(tab))]))
    # reads_dt = .extend_reads(reads_dt, 50)
    message("correlate coarse...")

    reads_dt$shiftn = 0
    reads_dt[ strand == "-", c("start", "end") := list(start + width, end + width)]
    corrVals = pbapply::pblapply(
        seq(from = frag_min, to = frag_max, by = step),
        function(shiftLen){
            ###TODO HERE
            # if(shiftLen == 100) browser()
            reads_dt[strand == "+", shiftn := shiftLen ]
            dc_dt = .calc_cross_corr(reads_dt[, list(which_label,
                                                seqnames,
                                                start = start + shiftn, end = end + shiftn,
                                                strand)],
                                   test_gr, frag_len = NA, window_size = 1)
            cov_dt = .calc_stranded_coverage(reads_dt[, list(which_label,
                                                          seqnames,
                                                          start = start + shiftn, end = end + shiftn,
                                                          strand)], test_gr, 1)
            ggplot(cov_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id")
            dc_dt$shiftLen = shiftLen
            dc_dt
        })
    rbindlist(corrVals)
}

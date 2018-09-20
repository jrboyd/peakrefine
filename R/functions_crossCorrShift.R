#' Title
#'
#' @param bam_file
#' @param qgr
#' @param nbest
#' @param revbest
#' @param qual_metric
#' @param max_dupes
#' @param shift_min
#' @param shift_max
#' @param step
#' @param include_plots
#'
#' @return
#' @export
#' @import pbapply
#' @import ggplot2
#' @examples
crossCorrByShift = function(bam_file, qgr, nbest = 20, revbest = FALSE,
                        qual_metric = "qValue", max_dupes = 1,
                        shift_min = 0, shift_max = 250,
                        step = 10, small_step = 1, include_plots = TRUE){
    # browser()
    if(is.na(nbest)){
        test_gr = qgr
    }else{
        test_gr = qgr[order(mcols(qgr)[[qual_metric]], decreasing = !revbest)][seq_len(min(nbest, length(qgr)))]
    }
    if(is.null(names(test_gr))){
        names(test_gr) = test_gr$id
    }

    message("fetch reads...")
    # browser()
    reads_dt = .fetch_bam_stranded(bam_file, test_gr, max_dupes = max_dupes)

    tab = table(reads_dt$width)
    read_length = as.numeric(names(tab[which(tab == max(tab))]))
    # reads_dt = strandsExtend(reads_dt, 50)
    message("correlate coarse...")

    reads_dt$shiftn = 0
    reads_dt[ strand == "-", c("start", "end") := .(start + width, end + width)]
    corrVals = pbapply::pblapply(
        seq(from = shift_min, to = shift_max, by = step),
        function(shiftLen){
            ###TODO HERE
            # if(shiftLen == 100) browser()
            reads_dt[strand == "+", shiftn := shiftLen ]
            dc_dt = calcStrandCorr(reads_dt[, .(which_label,
                                                seqnames,
                                                start = start + shiftn, end = end + shiftn,
                                                strand)],
                                   test_gr, fragLen = NA, win_size = 1)
            cov_dt = strandsCoverage(reads_dt[, .(which_label,
                                                  seqnames,
                                                  start = start + shiftn, end = end + shiftn,
                                                  strand)], test_gr, 1)
            ggplot(cov_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("id")
            dc_dt$shiftLen = shiftLen
            dc_dt
        })
    rbindlist(corrVals)
}

#' Title
#'
#' @param bam_file
#' @param qgr
#' @param nbest
#' @param qual_metric
#' @param max_dupes
#' @param frag_min
#' @param frag_max
#' @param step
#' @param include_plots
#'
#' @return
#' @export
#' @import pbapply
#' @depends
#' @examples
crossCorrByExtension = function(bam_file, qgr, nbest = 20,
                         qual_metric = "qValue", max_dupes = 1,
                         frag_min = 50, frag_max = 250,
                         step = 10, small_step = 1, include_plots = TRUE){
    if(is.na(nbest)){
        test_gr = qgr
    }else{
        test_gr = qgr[order(mcols(qgr)[[qual_metric]], decreasing = TRUE)][seq_len(min(nbest, length(qgr)))]
    }
    if(is.null(test_gr$id)){
        test_gr$id = paste0("peak_", seq_along(test_gr))
    }
    if(is.null(names(test_gr))){
        names(test_gr) = test_gr$id
    }
    # browser()
    message("fetch reads...")
    reads_dt = .fetch_bam_stranded(bam_file, test_gr, max_dupes = max_dupes)

    cnt_dt = reads_dt[, .N, by = .(which_label)]
    test_dt = data.table(which_label = as.character(test_gr), id = test_gr$id)
    cnt_dt = merge(cnt_dt, test_dt, all = TRUE)
    cnt_dt[is.na(N), N := 0]
    cnt_dt = cnt_dt[, .(id, count = N)]

    read_corr = calcStrandCorr(reads_dt, test_gr)
    read_coverage = strandsCoverage(reads_dt, test_gr)

    tab = table(reads_dt$width)
    read_length = as.numeric(names(tab[which(tab == max(tab))]))
    message("correlate coarse...")
    corrVals = pbapply::pblapply(seq(from = frag_min, to = frag_max, by = step), function(fragLen){
        dc_dt = calcStrandCorr(reads_dt, test_gr, fragLen)
        dc_dt$fragLen = fragLen
        dc_dt
    })

    corrVals = data.table::rbindlist(corrVals)

    # corrVals$crank = NULL
    corrVals[, crank := rank(-corr), by = .(id)]
    center = round(mean(corrVals[crank < 2 & !is.na(corr)]$fragLen))
    message("correlate fine...")
    corrValsDetail = pbapply::pblapply(seq(from = center-step, to = center+step, by = small_step), function(fragLen){
        dc_dt = calcStrandCorr(reads_dt, test_gr, fragLen)
        dc_dt$fragLen = fragLen
        dc_dt
    })
    corrValsDetail = rbindlist(corrValsDetail)
    corrValsDetail[, crank := rank(-corr), by = .(id)]
    corrValsDetail[crank == 1]
    bestFragLen = round(mean(corrValsDetail[crank < 2]$fragLen))

    frag_corr = corrValsDetail[fragLen == bestFragLen, 1:2]

    if(include_plots){
        tp = sample(unique(corrVals$id), min(12, length(test_gr)))
        message("plot sampled regions...")
        p = ggplot(corrVals[id %in% tp], aes(x = fragLen, y = corr, group = id)) + geom_path() +
            geom_path(data = corrValsDetail[id %in% tp], color = "red") + facet_wrap("id") +
            geom_point(data = corrValsDetail[id %in% tp][crank == 1], color = "red")
        out = list(
            read_length = read_length,
            frag_length = bestFragLen,
            read_corr = read_corr,
            frag_corr = frag_corr,
            count = cnt_dt,
            sample_plot = p
        )
    }else{
        out = list(
            read_length = read_length,
            frag_length = bestFragLen,
            read_corr = read_corr,
            frag_corr = frag_corr,
            count = cnt_dt
        )
    }
    return(out)
}

#' Title
#'
#' @param bam_file
#' @param qgr
#' @param fragLen
#' @param qual_metric
#' @param max_dupes
#' @param ncores
#'
#' @return
#' @export
#'
#' @examples
crossCorrByExtensionFull = function(bam_file, qgr, fragLen,
                                    qual_metric = "id",
                                    max_dupes = 1, ncores = 1, output_withGRanges = TRUE){
    # browser()
    options(mc.cores = ncores)
    assignments = ceiling(seq_along(qgr) / (length(qgr)/ncores))
    cres = mclapply(seq_len(ncores), function(i){
        crossCorrByExtension(bam_file,
                             qgr[assignments == i],
                             nbest = NA,
                             step = 0,
                             qual_metric = qual_metric,
                             max_dupes = max_dupes,
                             frag_min = fragLen,
                             frag_max = fragLen,
                             include_plots = FALSE
        )
    })
    out = list(
        read_corr =
            rbindlist(lapply(cres, function(x)x$read_corr)),
        frag_corr =
            rbindlist(lapply(cres, function(x)x$frag_corr)),
        count =
            rbindlist(lapply(cres, function(x)x$count))
    )
    colnames(out$read_corr)[2] = "read_corr"
    colnames(out$frag_corr)[2] = "frag_corr"
    out = merge(merge(out$read_corr, out$frag_corr), out$count)
    if(output_withGRanges){
        out = GRanges(merge(out, qgr, by = "id"))
    }
    return(out)

}

#' Calculate cross correlation by extending reads
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
#' @return named list of results
#' @export
#' @import pbapply
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#' crossCorrByExtension(bam_file, qgr[1:2], frag_min = 50,
#' frag_max = 250, step = 50, small_step = 10)
crossCorrByExtension = function(bam_file,
                                query_gr,
                                n_regions = 20,
                                max_dupes = 1,
                                frag_min = 50,
                                frag_max = 250,
                                step = 10,
                                small_step = 1,
                                include_plots = TRUE){
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

    # browser()
    message("fetch reads...")
    reads_dt = .fetch_bam_stranded(bam_file, test_gr, max_dupes = max_dupes)

    cnt_dt = reads_dt[, .N, by = .(which_label)]
    test_dt = data.table(which_label = as.character(test_gr), id = test_gr$id)
    cnt_dt = merge(cnt_dt, test_dt, all = TRUE)
    cnt_dt[is.na(N), N := 0]
    cnt_dt = cnt_dt[, .(id, count = N)]

    read_corr = calcStrandCorr(reads_dt, test_gr)
    read_coverage = .calc_stranded_coverage(reads_dt, test_gr)

    tab = table(reads_dt$width)
    read_length = as.numeric(names(tab[which(tab == max(tab))]))
    message("correlate coarse...")
    corrVals = pbapply::pblapply(seq(from = frag_min, to = frag_max, by = step), function(frag_len){
        dc_dt = calcStrandCorr(reads_dt, test_gr, frag_len)
        dc_dt$frag_len = fragLen
        dc_dt
    })

    corrVals = data.table::rbindlist(corrVals)

    # corrVals$crank = NULL
    corrVals[, crank := rank(-corr), by = .(id)]
    center = round(mean(corrVals[crank < 2 & !is.na(corr)]$frag_len))
    message("correlate fine...")
    corrValsDetail = pbapply::pblapply(seq(from = center-step, to = center+step, by = small_step), function(frag_len){
        dc_dt = calcStrandCorr(reads_dt, test_gr, frag_len)
        dc_dt$frag_len = fragLen
        dc_dt
    })
    corrValsDetail = rbindlist(corrValsDetail)
    corrValsDetail[, crank := rank(-corr), by = .(id)]
    corrValsDetail[crank == 1]
    bestFragLen = round(mean(corrValsDetail[crank < 2]$frag_len))

    frag_corr = corrValsDetail[frag_len == bestFragLen, 1:2]

    if(include_plots){
        tp = sample(unique(corrVals$id), min(12, length(test_gr)))
        message("plot sampled regions...")
        p = ggplot(corrVals[id %in% tp], aes(x = frag_len, y = corr, group = id)) + geom_path() +
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

#' Measure cross correlation using specified frag_len for all regions
#'
#' @param bam_file character. Path to .bam file, must have index at .bam.bai.
#' @param query_gr GRanges.  Regions to calculate cross correlation for.
#' @param frag_len integer. Fragment length to calculate cross correlation for.
#' @param max_dupes integer.  Duplicate reads above this value will be removed.
#' @param ncores integer.  ncores to use to split up the cross correlation
#'   calculation.
#' @param output_withGRanges logical.  Should results be merged back into
#'   query_gr? If TRUE output is GRanges. If FALSE output is data.table.
#'
#' @return Either a GRanges equivalent to query_gr with added columns for
#'   correlation metics or a data.table of metrics.
#' @export
#'
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#' crossCorrByExtensionFull(bam_file, qgr[1:2], frag_len = 150, ncores = 2)
crossCorrByExtensionFull = function(bam_file, query_gr, frag_len,
                                    max_dupes = 1,
                                    ncores = 1,
                                    output_withGRanges = TRUE){
    # browser()
    options(mc.cores = ncores)
    assignments = ceiling(seq_along(query_gr) / (length(query_gr)/ncores))
    cres = mclapply(seq_len(ncores), function(i){
        crossCorrByExtension(bam_file,
                             query_gr[assignments == i],
                             n_regions = NA,
                             step = 0,
                             max_dupes = max_dupes,
                             frag_min = frag_len,
                             frag_max = frag_len,
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
        out = GRanges(merge(out, query_gr, by = "id"))
    }
    return(out)

}

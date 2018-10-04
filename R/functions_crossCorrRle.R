#' Calculate cross correlation by extending reads
#'
#' @param bam_file character. Path to .bam file, must have index at .bam.bai.
#' @param query_gr GRanges.  Regions to calculate cross correlation for.
#' @param max_dupes integer.  Duplicate reads above this value will be removed.
#' @param fragment_range integer.  fragment size range to search for maximum
#'   correlation.
#' @param read_length integer. Any values outside fragment_range that must be
#'   searched.  If not supplied will be determined from bam_file.  Set as NA
#'   to disable this behavior.
#' @param include_plots logical. Should plots be included in output?
#'
#' @return named list of results
#' @export
#' @import GenomicRanges GenomicAlignments pbapply Rsamtools S4Vectors
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' query_gr = rtracklayer::import(np, format = "narrowPeak")
#' crossCorrByExtension(bam_file, query_gr[1:2], fragment_range = 50:250,
#' step = 50, small_step = 10)
crossCorrByRle = function(bam_file,
                          query_gr,
                          max_dupes = 1,
                          fragment_range = 50:300,
                          read_length = NULL,
                          include_plots = TRUE){
    rn = NULL # reserve for data.table
    if(is.null(query_gr$name)){
        if(is.null(names(query_gr))){
            query_gr$name = paste0("peak_", seq_along(query_gr))
        }else{
            query_gr$name = names(query_gr)
        }
    }else{
        if(is.null(names(query_gr))){
            names(query_gr) = query_gr$name
        }else{
            #both names() and $name are set, leave it alone
        }

    }
    names(query_gr) = query_gr$name
    # query_gr = resize(query_gr, 500, fix = "center")

    Param <- ScanBamParam(which=query_gr,
                          what=c("flag","mapq"))
    temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
    dt = as.data.table(temp)
    # browser()
    if(is.null(read_length)){
        read_length = getReadLength(bam_file, query_gr)
    }
    if(is.na(read_length)){
        read_length = numeric()
    }
    fragment_range = range(fragment_range)
    fragment_sizes = union(read_length, seq(fragment_range[1], fragment_range[2]))

    PosCoverage <- coverage(GenomicRanges::shift(GRanges(temp[strand(temp)=="+"])), -read_length)
    PosCoverage = PosCoverage[query_gr]
    names(PosCoverage) = query_gr$name

    NegCoverage <- coverage(GRanges(temp[strand(temp)=="-"]))
    NegCoverage = NegCoverage[query_gr]
    names(NegCoverage) = query_gr$name
    ShiftMatCor = pbapply::pbsapply(seq_along(query_gr), simplify = FALSE, function(i){
        ShiftsCorTemp <- S4Vectors::shiftApply(fragment_sizes,
                                               PosCoverage[[i]],
                                               NegCoverage[[i]],
                                               cor, simplify = FALSE,
                                               verbose = FALSE)
    })
    #necessary due to singleton query_gr or shift not resulting in matrix
    ShiftMatCor = matrix(unlist(ShiftMatCor),
                         byrow = FALSE,
                         nrow = length(fragment_sizes),
                         ncol = length(query_gr))
    colnames(ShiftMatCor) = query_gr$name
    rownames(ShiftMatCor) = fragment_sizes
    shift_dt = as.data.table(ShiftMatCor, keep.rownames = TRUE)
    shift_dt[, shift := as.numeric(rn)]
    shift_dt$rn = NULL
    shift_dt = melt(shift_dt, id.vars = "shift",
                    variable.name = "id", value.name = "correlation")
    return(shift_dt)
}



#' get the read length from a bam file
#'
#' @param bam_file character. Path to .bam file, must have index at .bam.bai.
#' @param query_gr GRanges.  Regions to calculate cross correlation for.  Only first 10 actually used.
#'
#' @return numeric.  the most common read length observed in bam.
#' @export
#'
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' query_gr = rtracklayer::import(np, format = "narrowPeak")
#' getReadLength(bam_file, query_gr)
#'
getReadLength = function(bam_file,
                         query_gr){
    Param <- Rsamtools::ScanBamParam(which=sample(query_gr, min(10, length(query_gr))),
                                     what=c("flag","mapq"))
    temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
    readlength=as.numeric(names(sort(table(width(temp)), decreasing = TRUE))[1])
    readlength
}

gather_metrics = function(peak_strand_corr, read_length = NULL){
    max_dt = peak_strand_corr[, .(shift = shift[which.max(correlation)], correlation = max(correlation)), by = .(id)]
    fl = round(mean(max_dt$shift, na.rm = TRUE))
    flex_frag_corrs = max_dt[, .(shift, id, correlation)]
    mean_frag_corrs = peak_strand_corr[shift == fl]

    if(!is.null(read_length)){
        read_corrs = peak_strand_corr[shift == read_length]
        out = list(read_length = read_length,
                   fragment_length = fl,
                   read_correlation = read_corrs,
                   flex_fragment_correlation = flex_frag_corrs,
                   mean_fragment_correlation = mean_frag_corrs,
                   full_correlation_results = peak_strand_corr)
    }else{
        out = list(fragment_length = fl,
                   flex_fragment_correlation = flex_frag_corrs,
                   mean_fragment_correlation = mean_frag_corrs,
                   full_correlation_results = peak_strand_corr)
    }
    out
}

#' Title
#'
#' @param bam_file
#' @param qgr
#' @param frag_min
#' @param frag_max
#' @param bam_md5
#' @param cache_path
#' @param cach_version
#' @param n_splits
#'
#' @return
#' @export
#'
#' @examples
calcCorrMetrics = function(bam_file, qgr, frag_min, frag_max,
                           bam_md5 = NULL,
                           cache_path = "~/.cache_peakrefine",
                           cach_version = "v1",
                           n_splits = getOption("mc.cores", 1L)){
    if(is.null(bam_md5)){
        bam_md5 = tools::md5sum(bam_file)
    }
    bfc_corr = BiocFileCache::BiocFileCache(cache_path, ask = FALSE)
    corr_key = paste(digest::digest(base_gr), bam_md5, frag_min, frag_max, cach_version, sep = "_")
    corr_res = bfcif(bfc_corr, corr_key, function(){
        message("cached results not found, gathering correlation info.")
        nper = ceiling(length(qgr) / n_splits)
        grps = ceiling(seq_along(qgr)/ nper)
        table(grps)

        rl = getReadLength(bam_file, qgr)
        lres = parallel::mclapply(unique(grps), function(g){
            k = grps == g
            crossCorrByRle(bam_file, qgr[k], fragment_range = c(frag_min, frag_max), read_length = rl)
        })
        peak_strand_corr = rbindlist(lres)
        gather_metrics(peak_strand_corr, rl)
        #
        # read_corrs = peak_strand_corr[shift == rl]
        # max_dt = peak_strand_corr[, .(shift = shift[which.max(correlation)], correlation = max(correlation)), by = .(id)]
        # fl = round(mean(max_dt$shift, na.rm = TRUE))
        # flex_frag_corrs = max_dt[, .(shift, id, correlation)]
        # mean_frag_corrs = peak_strand_corr[shift == fl]
        # list(read_length = rl,
        #      fragment_length = fl,
        #      read_correlation = read_corrs,
        #      flex_fragment_correlation = flex_frag_corrs,
        #      mean_fragment_correlation = mean_frag_corrs,
        #      full_correlation_results = peak_strand_corr)
    })
}

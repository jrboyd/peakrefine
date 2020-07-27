#' Calculate cross correlation by extending reads
#'
#' @param bam_file character. Path to .bam file, must have index at .bam.bai.
#' @param query_gr GRanges.  Regions to calculate cross correlation for.
#' @param max_dupes integer.  Duplicate reads above this value will be removed.
#' @param fragment_sizes integer.  fragment size range to search for maximum
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
                          fragment_sizes = 50:300,
                          read_length = NULL,
                          include_plots = TRUE,
                          ...){
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

    query_gr = harmonize_seqlengths(query_gr, bam_file)

    Param <- ScanBamParam(which=query_gr,
                          what=c("flag","mapq"),
                          ...)
    temp <- GenomicAlignments::readGAlignments(bam_file,param=Param)
    dt = as.data.table(temp)
    # browser()
    if(is.null(read_length)){
        read_length = getReadLength(bam_file, query_gr)
    }
    if(is.na(read_length)){
        read_length = numeric()
    }
    fragment_sizes = sort(union(read_length, fragment_sizes))

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
    ShiftMatCor[is.nan(ShiftMatCor)] = 0

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

#calculate mode
.my_mode = function(x){
    names(sort(-table(x)))[1]
}

gather_metrics = function(peak_strand_corr, read_length = NULL){
    max_dt = peak_strand_corr[, .(shift = shift[which.max(correlation)], correlation = max(correlation)), by = .(id)]
    # if(is.null(read_length)){
    #     fl = round(.my_mode(max_dt[shift != min(shift, na.rm = TRUE) & shift != max(shift, na.rm = TRUE)]$shift, na.rm = TRUE))
    # }else{
    #     fl = round(.my_mode(max_dt[shift != read_length][shift != min(shift, na.rm = TRUE) & shift != max(shift, na.rm = TRUE)]$shift, na.rm = TRUE))
    # }
    flex_frag_corrs = max_dt[, .(shift, id, correlation)]


    average_corr = peak_strand_corr[, .(correlation = mean(correlation)), .(shift)]

    fl = average_corr[, shift[which.max(correlation)[1]]]


    stable_frag_corrs = peak_strand_corr[shift == fl]


    if(!is.null(read_length)){
        read_corrs = peak_strand_corr[shift == read_length]
        out = list(read_length = read_length,
                   fragment_length = fl,
                   read_correlation = read_corrs,
                   flex_fragment_correlation = flex_frag_corrs,
                   stable_fragment_correlation = stable_frag_corrs,
                   full_correlation_results = peak_strand_corr,
                   average_correlation = average_corr)
    }else{
        out = list(fragment_length = fl,
                   flex_fragment_correlation = flex_frag_corrs,
                   stable_fragment_correlation = stable_frag_corrs,
                   full_correlation_results = peak_strand_corr,
                   average_correlation = average_corr)
    }
    out
}

#' Title
#'
#' @param bam_file
#' @param qgr
#' @param frag_sizes
#' @param fetch_size
#' @param bam_md5
#' @param qgr_md5
#' @param cache_path
#' @param cach_version
#' @param force_overwrite
#' @param n_splits
#'
#' @return
#' @export
#'
#' @examples
calcSCCMetrics = function(bam_file, qgr, frag_sizes, fetch_size = 3*max(frag_sizes),
                          bam_md5 = NULL, qgr_md5 = NULL,
                          cache_path = "~/.cache_peakrefine",
                          cach_version = "v1", force_overwrite = FALSE,
                          n_splits = getOption("mc.cores", 1L),
                          ...){
    if(is.null(bam_md5)){
        bam_md5 = tools::md5sum(bam_file)
    }
    if(is.null(qgr_md5)){
        qgr_md5 = digest::digest(qgr)
    }
    if(fetch_size <= max(frag_sizes)){
        stop("fetch_size (", fetch_size, ") must exceed max of frag_sizes (", max(frag_sizes), ").")
    }
    stopifnot(file.exists(bam_file))
    stopifnot(class(qgr) == "GRanges")
    if(!file.exists(paste0(bam_file, ".bai"))){
        stop("bam_file index not found. expected at ", paste0(bam_file, ".bai"),
             "\ntry running:\nsamtools index ", bam_file)
    }
    stopifnot(n_splits >= 1)

    qgr = resize(qgr, fetch_size, fix = 'center')


    bfc_corr = BiocFileCache::BiocFileCache(cache_path, ask = FALSE)
    corr_key = paste(qgr_md5, bam_md5, digest::digest(frag_sizes), fetch_size, cach_version, sep = "_")
    corr_res = bfcif(bfc_corr, corr_key, function(){
        message("cached results not found, gathering correlation info.")
        nper = ceiling(length(qgr) / n_splits)
        grps = ceiling(seq_along(qgr)/ nper)
        table(grps)
        # browser()
        rl = getReadLength(bam_file, qgr)
        lres = parallel::mclapply(unique(grps), function(g){
            k = grps == g
            crossCorrByRle(bam_file, qgr[k], fragment_sizes = frag_sizes, read_length = rl, ...)
        })
        peak_strand_corr = rbindlist(lres)
        gather_metrics(peak_strand_corr, rl)
    }, force_overwrite = force_overwrite)
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
calcCorrMetrics = function(bam_file, qgr, frag_min, frag_max, fetch_size = 3*frag_max,
                           bam_md5 = NULL, qgr_md5 = NULL,
                           cache_path = "~/.cache_peakrefine",
                           cach_version = "v1", force_overwrite = FALSE,
                           n_splits = getOption("mc.cores", 1L)){
    calcSCCMetrics(bam_file = bam_file,
                   qgr = qgr,
                   frag_sizes = seq(frag_min, frag_max),
                   fetch_size = fetch_size,
                   bam_md5 = bam_md5,
                   qgr_md5 = qgr_md5,
                   cache_path = cache_path,
                   cach_version = cach_version,
                   force_overwrite = force_overwrite,
                   n_splits = n_splits)
}

#' Calculate cross correlation by extending reads
#'
#' @param bam_file character. Path to .bam file, must have index at .bam.bai.
#' @param query_gr GRanges.  Regions to calculate cross correlation for.
#' @param max_dupes integer.  Duplicate reads above this value will be removed.
#' @param frag_min integer.  extension value to start at.
#' @param frag_max integer. extension value to end at.
#' @param include_plots logical. Should plots be included in output?
#'
#' @return named list of results
#' @export
#' @import GenomicRanges GenomicAlignments pbapply Rsamtools S4Vectors
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' query_gr = rtracklayer::import(np, format = "narrowPeak")
#' crossCorrByExtension(bam_file, query_gr[1:2], frag_min = 50,
#' frag_max = 250, step = 50, small_step = 10)
crossCorrByRle = function(bam_file,
                                query_gr,
                                max_dupes = 1,
                                frag_min = 50,
                                frag_max = 300,
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
    readlength=as.numeric(names(sort(table(width(temp)), decreasing = TRUE))[1])
    PosCoverage <- coverage(GenomicRanges::shift(GRanges(temp[strand(temp)=="+"])), -readlength)
    PosCoverage = PosCoverage[query_gr]
    names(PosCoverage) = query_gr$name

    NegCoverage <- coverage(GRanges(temp[strand(temp)=="-"]))
    NegCoverage = NegCoverage[query_gr]
    names(NegCoverage) = query_gr$name
    ShiftMatCor = pbapply::pbsapply(seq_along(query_gr), simplify = FALSE, function(i){
        ShiftsCorTemp <- S4Vectors::shiftApply(seq(frag_min,frag_max),
                                    PosCoverage[[i]],
                                    NegCoverage[[i]],
                                    cor, simplify = FALSE,
                                    verbose = FALSE)
    })
    #necessary due to singleton query_gr or shift not resulting in matrix
    ShiftMatCor = matrix(unlist(ShiftMatCor),
                         byrow = FALSE,
                         nrow = length(seq(frag_min,frag_max)),
                         ncol = length(query_gr))
    colnames(ShiftMatCor) = query_gr$name
    rownames(ShiftMatCor) = seq(frag_min,frag_max)
    shift_dt = as.data.table(ShiftMatCor, keep.rownames = TRUE)
    shift_dt[, shift := as.numeric(rn)]
    shift_dt$rn = NULL
    shift_dt = melt(shift_dt, id.vars = "shift",
                    variable.name = "id", value.name = "correlation")
return(shift_dt)
    # if(include_plots){
    #     tp = sample(unique(corrVals$id), min(12, length(test_gr)))
    #     message("plot sampled regions...")
    #     p = ggplot(corrVals[id %in% tp], aes(x = frag_len, y = corr, group = id)) + geom_path() +
    #         geom_path(data = corrValsDetail[id %in% tp], color = "red") + facet_wrap("id") +
    #         geom_point(data = corrValsDetail[id %in% tp][crank == 1], color = "red")
    #     out = list(
    #         read_length = read_length,
    #         frag_length = bestFragLen,
    #         read_corr = read_corr,
    #         frag_corr = frag_corr,
    #         corr_vals = corrVals,
    #         count = cnt_dt,
    #         sample_plot = p
    #     )
    # }else{
    #     out = list(
    #         read_length = read_length,
    #         frag_length = bestFragLen,
    #         read_corr = read_corr,
    #         frag_corr = frag_corr,
    #         corr_vals = corrVals,
    #         count = cnt_dt
    #     )
    # }
    # return(out)
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

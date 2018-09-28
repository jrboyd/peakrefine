eval_newFUN = T
eval_equal = FALSE
eval_speed = F


if(eval_newFUN){
    library(Rsamtools)
    library(data.table)
    library(ggplot2)
    library(GenomicRanges)
    library(peakrefine)

    qgr = seqsetvis::easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled_peaks_passIDR.05.narrowPeak")[[1]]
    qgr = dropSeqlevels(qgr, "chrU13369.1")
    qgr = sample(qgr, 3200)
    qgr = resize(qgr, 800, fix = "center")
    bam_file = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled.bam"

    njobs = 32
    nper = ceiling(length(qgr) / njobs)
    grps = ceiling(seq_along(qgr)/ nper)
    table(grps)
    options(mc.cores = 32)
    getOption("mc.cores", 2L)
    st = system.time({
        lres = parallel::mclapply(unique(grps), function(g){
            k = grps == g
            crossCorrByRle(bam_file, qgr[k])
        })
    })

    rl = get_readLength(bam_file, qgr)

    cce = rbindlist(lres)
    read_corrs = cce[shift == rl]
    ggplot(cce[id %in% unique(id)[1:50]],
           aes(x = shift, y = correlation, group = id)) +
        geom_path(alpha = .1, size = 3)
}
if(eval_equal){ #correlation of Rle look equivalent
    dat = cbind(c(0,0, 1, 2, 3, 2, 1, 0, 0, 2, 2, 1),
                c(1,0,1,2,4,2,1,0,0, 1, 1, 2))
    cor(dat)[1,2]
    cor(Rle(dat[,1]), Rle(dat[,2]))
}

if(eval_speed){ #correlation by Rle is significantly faster
    library(Rsamtools)
    library(data.table)
    library(ggplot2)
    library(GenomicRanges)

    qgr = seqsetvis::easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled_peaks_passIDR.05.narrowPeak")[[1]]
    qgr = dropSeqlevels(qgr, "chrU13369.1")
    qgr = sample(qgr, 5)
    names(qgr) = qgr$name
    qgr$id = qgr$name


    bam_file = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled.bam"

    ###ChIPQC param setup

    bamFile = bam_file
    shiftWindowStart=50
    shiftWindowEnd=300

    qgr = resize(qgr, shiftWindowEnd*2, fix = "center")

    tp = qgr$id[1:20]

    ### dt cross corr
    library(peakrefine)
    shift_corr = crossCorrByShift(
        bam_file, qgr, length(qgr),
        step = 20, frag_min = shiftWindowStart, frag_max = shiftWindowEnd,
        max_dupes = Inf)
    shift_corr$id = factor(shift_corr$id, levels = qgr$id)
    # ggplot(shift_corr[id %in% tp], aes(x = shiftLen, y = corr, color = id)) + facet_wrap("id") + geom_path() + guides(color = "none")
    ext_corr_res = crossCorrByExtension(
        bam_file, qgr, length(qgr),
        frag_min = shiftWindowStart, frag_max = shiftWindowEnd, step = 20,
        small_step = 10, include_plots = FALSE, max_dupes = Inf)
    ext_corr = ext_corr_res$corr_vals
    # ggplot(ext_corr[id %in% tp], aes(x = frag_len, y = corr, color = id)) + facet_wrap("id") + geom_path() + guides(color = "none")

    ### ChIPQC corr method
    # ChrLengths <- scanBamHeader(bamFile)[[1]]$targets

    # ShiftMatCor <- NULL

    Param <- ScanBamParam(which=qgr,
                          what=c("flag","mapq"))
    temp <- GenomicAlignments::readGAlignments(bamFile,param=Param)
    readlength=as.numeric(names(sort(table(width(temp)), decreasing = TRUE))[1])
    PosCoverage <- coverage(GenomicRanges::shift(GRanges(temp[strand(temp)=="+"])), -readlength)
    PosCoverage = PosCoverage[qgr]
    names(PosCoverage) = qgr$name

    NegCoverage <- coverage(GRanges(temp[strand(temp)=="-"]))
    NegCoverage = NegCoverage[qgr]
    names(NegCoverage) = qgr$name
    ShiftMatCor = pbapply::pbsapply(seq_along(qgr), function(i){
        ShiftsCorTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),
                                    PosCoverage[[i]],NegCoverage[[i]],cor, verbose = FALSE)
    })

    colnames(ShiftMatCor) = qgr$name
    rownames(ShiftMatCor) = seq(shiftWindowStart,shiftWindowEnd)
    shift_dt = as.data.table(ShiftMatCor, keep.rownames = TRUE)
    shift_dt[, shift := as.numeric(rn)]
    shift_dt$rn = NULL
    shift_dt = melt(shift_dt, id.vars = "shift",
                    variable.name = "id", value.name = "correlation")
    # ggplot(shift_dt[id %in% tp][order(shift)], aes(x = shift, color = id, y = correlation, group = id)) + geom_path() + facet_wrap("id") + guides(color = "none")
    pdt = rbindlist(list(
        shift_corr[, .(id, corr, fl = shiftLen, method = "dt_shift")],
        ext_corr[, .(id, corr, fl = frag_len, method = "dt_extension")],
        shift_dt[, .(id, corr = correlation, fl = shift, method = "gr_shift")]
    ))

    ggplot(pdt[id %in% tp],
           aes(x = fl, color = method, y = corr, group = method)) +
        geom_path() + facet_wrap("id")# + guides(color = "none")

    max_dt = pdt[, .(fl = fl[which.max(corr)], corr = max(corr)), by = .(id, method)]
    max_dt = dcast(max_dt, id ~ method, value.var = "fl")
    plot(max_dt[, 2:4], xlim = c(50, 300), ylim = c(50, 300), pch = 16, col = rgb(0,0,0,.2))
}


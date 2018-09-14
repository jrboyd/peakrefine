# source("https://bioconductor.org/biocLite.R")
# biocLite("PWMEnrich")
# biocLite("PWMEnrich.Hsapiens.background")
# biocLite("BSgenome.Hsapiens.UCSC.hg38")
setwd("~/R/peakrefine/")
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)

library(seqsetvis)
library(data.table)
library(peakrefine)
data("PWMLogn.hg19.MotifDb.Hsap")

ncores = 32

registerCoresPWMEnrich(ncores)
useBigMemoryPWMEnrich(TRUE)

#
# motifs.denovo = readMotifs(file = "MA0511.1.jaspar")
# genomic.acgt = getBackgroundFrequencies("hg19")
#
# pwms.denovo = toPWM(motifs.denovo, prior=genomic.acgt)
# bg.denovo = makeBackground(pwms.denovo, organism="hg19", type="logn", quick=TRUE)
#
#
# PWMLogn.hg19.MotifDb.Hsap@pwms



print(load("../KZ_Runx1_BRG1_ESR1_overlap/corr_MCF7_runx1_and_esr1.save"))
setwd("/slipstream/galaxy/uploads/working/qc_framework/output_JR_bookmarking_blocked_RUNX1_U13369masked/")

qgr = easyLoad_narrowPeak("MCF10A-blocked_Runx1-4336BF_pooled/MCF10A-blocked_Runx1-4336BF_pooled_peaks.narrowPeak")[[1]]
bam_file = "MCF10A-blocked_Runx1-4336BF_pooled/MCF10A-blocked_Runx1-4336BF_pooled.bam"

bam_file = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled.bam"
bam_input = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_input_pooled/MDA231_input_pooled.bam"
qgr = easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled_peaks.narrowPeak")[[1]]

qgr$id = paste0("region_", seq_along(qgr))
names(qgr) = qgr$id
bgr = qgr

# cc = ssvCrossCorr(bam_file, bgr, shift_min = 0, shift_max = 300, step = 10, nbest = 16)
# corr_res = ssvStrandCorr(bam_file, bgr, frag_min = 1, frag_max = 300, nbest = 16)
#
# roi150 = c(1475, 1477, 1478, 1484)
# roi100 = c(1505, 1510, 828)
#
# p1 = ggplot(cc, aes(x = shiftLen, y = corr)) +
#     facet_wrap("id") +
#     geom_path() +
#     labs(title = "cross corr")
# p2 = corr_res$sample_plot +
#     labs(title = "frag corr")


score_motif = function(bam_file, qgr, fl, ncores = 16){
    ncores = 32

    registerCoresPWMEnrich(ncores)
    useBigMemoryPWMEnrich(TRUE)
}

all_res = list()

# fl = corr_res$frag_length
# fl = 180
todo_fl = c(180, 80, 130, 150, 230, 280, 330)
# rl = corr_res$read_length

### top enriched can be cached but score and pval are group calcs
subset_MotifEnrichmentResults = function(me_res, k){
    res = me_res@res
    res$sequences = res$sequences[k]
    res$sequence.nobg = res$sequence.nobg[k, ]
    res$sequence.bg = res$sequence.bg[k, ]
    res$sequence.norm = res$sequence.norm[k, ]

    sequences = res$sequences
    score = res$score
    #TODO fix HARDCODE
    pwmobj = PWMLogn.hg19.MotifDb.Hsap
    pwms = res$pwms
    verbose = FALSE
    cutoff = NULL
    bg = res$bg
    group.only = FALSE

### from motifEnrichment
    if (score == "affinity") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        # res$sequence.nobg = motifScores(sequences, pwms, verbose = verbose)
        res$group.nobg = PWMEnrich:::affinitySequenceSet(res$sequence.nobg,
                                             seq.len, pwm.len)
    } else if (score == "clover") {
        # res$sequence.nobg = motifScores(sequences, pwms, verbose = verbose)
        res$group.nobg = PWMEnrich:::cloverScore(res$sequence.nobg, verbose = verbose)
    } else if (score == "cutoff") {
        res$params = list(cutoff = cutoff)
        # res$sequence.nobg = motifScores(sequences, pwms, cutoff = cutoff,
                                        # verbose = verbose)
        res$group.nobg = colSums(res$sequence.nobg)
    } else {
        stop(paste("Unknown scoring algorithm: '", score, "'. Please select one of: 'affinity', 'cutoff', 'clover'",
                   sep = ""))
    }

    if (bg == "none") {
        # res$sequence.bg = NULL
        res$group.bg = NULL
    } else if (bg == "logn") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        # res$sequence.bg = logNormPval(res$sequence.nobg, seq.len,
        #                               pwm.len, pwmobj@bg.mean, pwmobj@bg.sd, pwmobj@bg.len)
        # colnames(res$sequence.bg) = names(pwms)
        res$sequence.norm = apply(res$sequence.bg, 1:2, qlnorm,
                                  lower.tail = FALSE)
        if (score == "affinity") {
            if (is.matrix(pwmobj@bg.mean)) {
                res$group.bg = apply(res$sequence.bg, 2, function(x) {
                    pchisq(-2 * sum(log(x)), 2 * length(x), lower.tail = FALSE)
                })
                res$group.norm = sapply(res$group.bg, qlnorm,
                                        lower.tail = FALSE)
            } else {
                res$group.bg = logNormPvalSequenceSet(res$sequence.nobg,
                                                      seq.len, pwm.len, pwmobj@bg.mean, pwmobj@bg.sd,
                                                      pwmobj@bg.len)
                res$group.norm = sapply(res$group.bg, qlnorm,
                                        lower.tail = FALSE)
            }
        } else if (score == "clover") {
            res$group.bg = cloverScore(res$sequence.norm, verbose = verbose)
        }
    } else if (bg == "z") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        # res$sequence.bg = cutoffZscore(res$sequence.nobg, seq.len,
        #                                pwm.len, pwmobj@bg.P)
        res$group.bg = cutoffZscoreSequenceSet(res$sequence.nobg,
                                               seq.len, pwm.len, pwmobj@bg.P)
    } else if (bg == "pval") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        usecutoff = NULL
        if (score == "cutoff")
            usecutoff = cutoff
        if (!group.only) {
            # res$sequence.bg = t(sapply(1:length(seq.len), function(i) empiricalPvalue(res$sequence.nobg[i,
            #                                                                                             ], seq.len[i], pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev,
            #                                                                           cutoff = usecutoff, B = B, verbose = verbose)))
        }
        if (score == "clover")
            res$group.bg = cloverPvalue1seq(res$sequence.nobg,
                                            seq.len, pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev,
                                            B = B, verbose = verbose, clover = res$group.nobg)
        else res$group.bg = empiricalPvalueSequenceSet(res$sequence.nobg,
                                                       seq.len, pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev, cutoff = usecutoff,
                                                       B = B, verbose = verbose)
    } else if (bg == "ms") {
        usecutoff = NULL
        if (score == "cutoff")
            usecutoff = cutoff
        # res$sequence.bg = matrixShuffleZscorePerSequence(res$sequence.nobg,
        #                                                  sequences, pwms, cutoff = usecutoff, B = motif.shuffles)
    } else if (bg == "gev") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        # res$sequence.bg = gevPerSequence(res$sequence.nobg, seq.len,
        #                                  pwm.len, pwmobj@bg.loc, pwmobj@bg.scale, pwmobj@bg.shape)
        res$group.bg = NULL
    } else {
        stop(paste("Uknown background correction algorithm: '",
                   bg, "', Please select one of: 'none', 'logn', 'z', 'pval', 'ms'",
                   sep = ""))
    }

    if (group.only) {
        seq.n = grep("^sequence[.]", names(res))
        if (length(seq.n) > 0) {
            res = res[-seq.n]
        }
    }
    if ("sequence.nobg" %in% names(res)) {
        rownames(res$sequence.nobg) = names(sequences)
    }
    if ("sequence.bg" %in% names(res)) {
        rownames(res$sequence.bg) = names(sequences)
    }

    me_res@res = res
    return(me_res)
}

for(fl in todo_fl){
    message("fragLen ", fl)
    odir = file.path("~/R/peakrefine/cache", sub(".bam", "", basename(bam_file)))
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)
    setwd(odir)

    # cowplot::plot_grid(p1 +
    #                        annotate("line", x = rep(fl, 2), y = c(-.5, 1), color = "green") +
    #                        annotate("line", x = rep(rl, 2), y = c(-.5, 1), color = "red"),
    #                    p2 +
    #                        annotate("line", x = rep(fl, 2), y = c(-.5, 1), color = "green") +
    #                        annotate("line", x = rep(rl, 2), y = c(-.5, 1), color = "red")
    # )
    cfile = paste0("cache_strandRes_", fl, ".save")
    if(file.exists(cfile)){
        load(cfile)
    }else{
        message("calculating strandRes...")
        strandRes = ssvStrandCorrFull(bam_file, qgr, fragLen = fl, ncores = ncores, output_withGRanges = TRUE)
        save(strandRes, file = cfile)
    }

    cfile = paste0("cache_strandResInput_", fl, ".save")
    if(file.exists(cfile)){
        load(cfile)
    }else{
        message("calculating strandResInput...")
        strandResInput = ssvStrandCorrFull(bam_input, qgr, fragLen = fl, ncores = ncores, output_withGRanges = TRUE)
        save(strandResInput, file = cfile)
    }

    corrInput_dt = as.data.table(strandResInput)
    corrInput_dt = corrInput_dt[, .(id,
                     input_read_corr = read_corr,
                     input_frag_corr = frag_corr,
                     input_count = count)]

    corr_dt = as.data.table(strandRes)
    corr_dt = merge(corr_dt, corrInput_dt)
    corr_dt = corr_dt[seqnames != "chrU13369.1"]
    corr_dt[, diff_corr := frag_corr - read_corr]
    #stratify by various metrics
    ngroup = 8
    g = factor(paste0("g", seq_len(ngroup)), levels = paste0("g", seq_len(ngroup)))
    ##
    corr_dt[, diff_group := g[ceiling(rank(-diff_corr) / .N * ngroup)]]
    corr_dt[, frag_group := g[ceiling(rank(-frag_corr) / .N * ngroup)]]
    ## fragment length independent metrics
    if(fl == todo_fl[1]){
        corr_dt[, read_group := g[ceiling(rank(-read_corr) / .N * ngroup)]]
        corr_dt[, signalValue_group := g[ceiling(rank(-signalValue) / .N * ngroup)]]
        corr_dt[, pValue_group := g[ceiling(rank(-pValue) / .N * ngroup)]]
    }
    ## input metrics
    corr_dt[, input_read_group := g[ceiling(rank(-input_read_corr) / .N * ngroup)]]
    corr_dt[, input_frag_group := g[ceiling(rank(-input_frag_corr) / .N * ngroup)]]
    corr_dt[, input_count_group := g[ceiling(rank(-input_count) / .N * ngroup)]]

    # plot(sort(corr_dt$frag_corr - corr_dt$read_corr))
    # plot(sort(corr_dt$frag_corr))
    # plot(corr_dt$frag_corr, corr_dt$pValue)
    # plot(corr_dt$diff_corr, corr_dt$pValue)

    todo_groups = colnames(corr_dt)[grepl("group", colnames(corr_dt))]

    myMotif = function(qgr){
        qgr = resize(qgr, 200, fix = "center")
        # qgr = GenomicRanges::shift(qgr, 600)

        sequences = getSeq(Hsapiens, qgr)

        res19 = motifEnrichment(sequences, PWMLogn.hg19.MotifDb.Hsap)
        groupReport(res19)
    }


    all_motif = lapply(todo_groups, function(grp){
        lapply(levels(corr_dt[[grp]]), function(sel){
            cfile = paste0("cache_", grp, "_", sel, "_", ngroup, "_", fl, ".save")
            if(file.exists(cfile)){
                load(cfile)
            }else{
                message("calculating motifs for ",
                        grp, " ", sel)
                qgr = GRanges(corr_dt[corr_dt[[grp]] == sel])
                myMotifRes = myMotif(qgr)
                save(myMotifRes, file = cfile)
            }
            myMotifRes
        })
    })
    setwd("~/R/peakrefine/")

    names(all_motif) = todo_groups

    all_motif$signalValue_group

    dt_motif = rbindlist(use.names = TRUE, idcol = "metric",
                         lapply(all_motif, function(x){
                             names(x) = g
                             x
                             rbindlist(use.names = TRUE, idcol = "group",
                                       lapply(x, function(y){
                                           dt = as.data.table(as.data.frame(y))
                                           colnames(dt) = gsub("\\.", "_", colnames(dt))
                                           dt
                                       }))
                         }))

    dt_motif$group = factor(dt_motif$group, levels = rev(levels(g)))
    dt_motif = dt_motif[order(group)][order(metric)][order(id)]
    # ggplot(dt_motif[grepl("RUNX", id)], aes(x = group, y = top_motif_prop, color = metric, group =  paste(id, metric))) +
    #     geom_path() +
    #     facet_grid("id~metric") +
    #     theme(panel.grid.major.y = element_line())

    all_res[[paste("fl_", fl)]] = dt_motif

    # ggplot(dt_motif[id %in% c("Hsapiens-jolma2013-RUNX2-3")], aes(x = group, y = top_motif_prop, color = metric, group =  paste(id, metric))) +
    #     geom_path() +
    #     facet_wrap("id") +
    #     # facet_grid("id~metric") +
    #     theme(panel.grid.major.y = element_line())
}

# all_res
#
# ggplot(all_res$`fl_ 180`[id %in% c("Hsapiens-jolma2013-RUNX2-3")], aes(x = group, y = top_motif_prop, color = metric, group =  paste(id, metric))) +
#     geom_path() +
#     facet_wrap("id") +
#     # facet_grid("id~metric") +
#     theme(panel.grid.major.y = element_line())
#
# ggplot(all_res$`fl_ 280`[id %in% c("Hsapiens-jolma2013-RUNX2-3")], aes(x = group, y = top_motif_prop, color = metric, group =  paste(id, metric))) +
#     geom_path() +
#     facet_wrap("id") +
#     # facet_grid("id~metric") +
#     theme(panel.grid.major.y = element_line())

dt = rbindlist(all_res, use.names = TRUE, idcol = "fragLen")
dt$fragLen = factor(dt$fragLen, levels = paste("fl_", sort(todo_fl)))
mycol = c("gray", "black", "red")[c(1, 1, 3, 2, 2, 2, 1)]
names(mycol) = paste("fl_", sort(todo_fl))
ggplot(dt[id %in% c("Hsapiens-jolma2013-RUNX2-3")], aes(x = group, y = top_motif_prop, color = fragLen, group =  paste(id, metric, fragLen))) +
    geom_path() +
    # facet_wrap("id") +
    facet_grid("id~metric") +
    theme(panel.grid.major.y = element_line()) +
    scale_color_manual(values = mycol) + theme_classic()

# ggplot(dt[id %in% c("Hsapiens-jolma2013-RUNX2-3")], aes(x = group, y = -log10(p_value), color = fragLen, group =  paste(id, metric, fragLen))) +
#     geom_path() +
#     # facet_wrap("id") +
#     facet_grid("id~metric") +
#     theme(panel.grid.major.y = element_line())

ggplot(dt[id %in% c("Hsapiens-jolma2013-RUNX2-3")], aes(x = group, y = raw_score, color = fragLen, group =  paste(id, metric, fragLen))) +
    geom_path() +
    # facet_wrap("id") +
    facet_grid("id~metric") +
    theme(panel.grid.major.y = element_line()) +
    scale_color_manual(values = mycol) + theme_classic()

# qgr = GRanges(runx_corr[signalValue > 10 & qValue > 10 & fragment > .8 & diff > .1])
# qgr = qgr[order(qgr$pValue, decreasing = TRUE)][1:2000]


#
# rep19 = as.data.table(as.data.frame(rep19raw))
# gsub("\\.", "_", colnames(rep19))
# rep19[order(`top.motif.prop`, decreasing = TRUE)]
# rep19[grepl("RUN",rep19$target)]
#
# plot(rep19raw[1:10], fontsize=7, id.fontsize=5)
#
#
# res.denovo = motifEnrichment(sequences, bg.denovo)
# rep = as.data.frame(sequenceReport(res.denovo, seq_along(qgr)))
# head(rep)
# ggplot(rep, aes(x = raw.score)) + geom_histogram()
# ggplot(rep, aes(x = p.value)) + geom_histogram()
# groupReport(res.denovo)
#
# score = motifScores(sequences, motifs.denovo, raw.scores = TRUE)
# # score
#
#
# plotMotifScores(score[1:10])

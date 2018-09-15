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

pre_motif = function(qgr, pwm = PWMLogn.hg19.MotifDb.Hsap, nbases = 200, ncores = 32){
    #TODO biocfilecache and digest args to auto cache, helper get cache file function
    registerCoresPWMEnrich(ncores)
    useBigMemoryPWMEnrich(TRUE)
    sequences = getSeq(Hsapiens,
                       resize(qgr, nbases, fix = "center"))
    motif_res = motifEnrichment(sequences, pwm)
    return(motif_res)
}

motif_res = pre_motif(qgr)

score_motif = function(bam_file, qgr, fl, motif_res = NULL,
                       pwm = PWMLogn.hg19.MotifDb.Hsap,
                       include_fl_independent = FALSE,
                       out_dir = file.path("~/R/peakrefine/cache", sub(".bam", "", basename(bam_file))),
                       nbases = 200, ncores = 8){
    group = metric = id = NULL #declare data.table bindings
    registerCoresPWMEnrich(ncores)
    useBigMemoryPWMEnrich(TRUE)

    message("fragLen ", fl)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(out_dir)

    cfile = paste0("cache_strandRes_", fl, ".save")
    if(file.exists(cfile)){
        load(cfile)
    }else{
        message("calculating strandRes...")
        strandRes = ssvStrandCorrFull(bam_file,
                                      qgr,
                                      fragLen = fl,
                                      ncores = ncores,
                                      output_withGRanges = TRUE)
        save(strandRes, file = cfile)
    }

    cfile = paste0("cache_strandResInput_", fl, ".save")
    if(file.exists(cfile)){
        load(cfile)
    }else{
        message("calculating strandResInput...")
        strandResInput = ssvStrandCorrFull(bam_input,
                                           qgr,
                                           fragLen = fl,
                                           ncores = ncores,
                                           output_withGRanges = TRUE)
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
    g = factor(paste0("g", seq_len(ngroup)),
               levels = paste0("g", seq_len(ngroup)))
    ##
    corr_dt[, diff_group := g[ceiling(rank(-diff_corr) / .N * ngroup)]]
    corr_dt[, frag_group := g[ceiling(rank(-frag_corr) / .N * ngroup)]]
    ## fragment length independent metrics
    if(include_fl_independent){
        corr_dt[, read_group := g[ceiling(rank(-read_corr) / .N * ngroup)]]
        corr_dt[, signalValue_group := g[ceiling(rank(-signalValue) / .N * ngroup)]]
        corr_dt[, pValue_group := g[ceiling(rank(-pValue) / .N * ngroup)]]
    }
    ## input metrics
    corr_dt[, input_read_group := g[ceiling(rank(-input_read_corr) / .N * ngroup)]]
    corr_dt[, input_frag_group := g[ceiling(rank(-input_frag_corr) / .N * ngroup)]]
    corr_dt[, input_count_group := g[ceiling(rank(-input_count) / .N * ngroup)]]

    todo_groups = colnames(corr_dt)[grepl("group", colnames(corr_dt))]

    ## loading full motifs.
    cfile = paste0("cache_motif_res_", nbases, "nbases_", length(qgr), "seq.save")
    if(is.null(motif_res)){
        if(file.exists(cfile)){
            message("loading cached full motif enrichment.")
            load(cfile)
        }else{
            message("calculating full motif enrichment...")
            sequences = getSeq(Hsapiens,
                               resize(qgr, nbases, fix = "center"))
            motif_res = motifEnrichment(sequences, pwm)
            save(motif_res, file = cfile)
            message("saved cached full motif enrichment.")
        }
    }else if(is.character(motif_res)){
        yes_cache = motif_res != cfile
        if(yes_cache){
            message("loading externally cached full motif enrichment.")
        }else{
            message("loading cached full motif enrichment.")
        }

        load(motif_res)
        if(!class(motif_res) == "MotifEnrichmentResults"){
            stop("cached file was not MotifEnrichmentResults!")
        }
        if(yes_cache){
            save(motif_res, file = cfile)
            message("caching full motif enrichment.")
        }
    }else{
        if(!class(motif_res) == "MotifEnrichmentResults"){
            stop("motif_res must be one of:
                 NULL - MotifEnrichmentResults will be calculated - SLOW
                 character - file path to cached MotifEnrichmentResults
                 MotifEnrichmentResults - object of class MotifEnrichmentResults")
        }
        save(motif_res, file = cfile)
        message("caching full motif enrichment.")
    }


    all_motif = lapply(todo_groups, function(grp){
        lapply(levels(corr_dt[[grp]]), function(sel){
            cfile = paste0("cache_", grp, "_", sel, "_", ngroup, "_", fl, ".save")
            if(file.exists(cfile)){
                load(cfile)
            }else{
                message("calculating motifs for ",
                        grp, " ", sel)
                # qgr = GRanges(corr_dt[corr_dt[[grp]] == sel])
                k = which(corr_dt[[grp]] == sel)
                myMotifRes = groupReport(subset_MotifEnrichmentResults(motif_res, k))
                # myMotifRes = myMotif(qgr, nbases)
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
                             rbindlist(use.names = TRUE, idcol = "group",
                                       lapply(x, function(y){
                                           dt = as.data.table(as.data.frame(y))
                                           colnames(dt) = gsub("\\.", "_", colnames(dt))
                                           dt
                                       }))
                         }))

    dt_motif$group = factor(dt_motif$group, levels = rev(levels(g)))
    dt_motif = dt_motif[order(group)][order(metric)][order(id)]

    dt_motif
}

### top enriched can be cached but score and pval are group calcs
subset_MotifEnrichmentResults = function(me_res, k){
    # browser()
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
        # res$sequence.bg = PWMEnrich:::logNormPval(res$sequence.nobg, seq.len,
        #                               pwm.len, pwmobj@bg.mean, pwmobj@bg.sd, pwmobj@bg.len)
        # colnames(res$sequence.bg) = names(pwms)
        # res$sequence.norm = apply(res$sequence.bg, 1:2, qlnorm,
        #                           lower.tail = FALSE)
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

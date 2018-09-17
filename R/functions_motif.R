# # source("https://bioconductor.org/biocLite.R")
# # biocLite("PWMEnrich")
# # biocLite("PWMEnrich.Hsapiens.background")
# # biocLite("BSgenome.Hsapiens.UCSC.hg38")
# setwd("~/R/peakrefine/")
# library(PWMEnrich)
# library(PWMEnrich.Hsapiens.background)
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(BSgenome)
#
# library(seqsetvis)
# library(data.table)
# library(peakrefine)
# data("PWMLogn.hg19.MotifDb.Hsap")

#' Calculate a large motif enrichment once for later rounds of
#' subsetting and analysis
#'
#' @param qgr
#' @param pwm
#' @param nbases
#' @param ncores
#'
#' @return
#' @export
#' @importClassesFrom PWMEnrich MotifEnrichmentResults
#'
#' @examples
pre_motif = function(qgr, pwm = PWMLogn.hg19.MotifDb.Hsap, nbases = 200, ncores = 32){
    #TODO biocfilecache and digest args to auto cache, helper get cache file function
    registerCoresPWMEnrich(ncores)
    useBigMemoryPWMEnrich(TRUE)
    sequences = getSeq(Hsapiens,
                       resize(qgr, nbases, fix = "center"))
    motif_res = motifEnrichment(sequences, pwm)
    return(motif_res)
}

#' Assembles several metrics for a set of peaks from:
#' 1) the original peak call
#' 2) the ChIP-seq bam file
#' 3) an optional input (no pull-down) bam file
#'
#' Stratify by each metric and conduct motif enrichment for each.
#'
#' @param bam_file
#' @param bam_input
#' @param qgr
#' @param fl
#' @param motif_res
#' @param pwm
#' @param include_fl_independent
#' @param out_dir
#' @param nbases
#' @param ncores
#'
#' @return
#' @export
#'
#' @examples
score_motif = function(bam_file, bam_input, qgr, fl, motif_res = NULL,
                       pwm = PWMLogn.hg19.MotifDb.Hsap,
                       include_fl_independent = FALSE,
                       out_dir = file.path("~/R/peakrefine/cache", sub(".bam", "", basename(bam_file))),
                       nbases = 200, ncores = 32){
    group = metric = id = NULL #declare data.table bindings
    registerCoresPWMEnrich(ncores)
    useBigMemoryPWMEnrich(TRUE)
    options(mc.cores = ncores)

    message("fragLen ", fl)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(out_dir)
    # browser()
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
        if(!file.exists(cfile)){
            save(motif_res, file = cfile)
            message("caching full motif enrichment.")
        }

    }

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
    set.seed(0)
    corr_dt[, diff_group := g[ceiling(rank(-diff_corr, ties = "random") / .N * ngroup)]]
    corr_dt[, frag_group := g[ceiling(rank(-frag_corr, ties = "random") / .N * ngroup)]]

    ## fragment length independent metrics
    if(include_fl_independent){
        corr_dt[, input_count_group := g[ceiling(rank(-input_count, ties = "random") / .N * ngroup)]]
        corr_dt[, count_group := g[ceiling(rank(-count, ties = "random") / .N * ngroup)]]
        corr_dt[, read_group := g[ceiling(rank(-read_corr, ties = "random") / .N * ngroup)]]
        corr_dt[, signalValue_group := g[ceiling(rank(-signalValue, ties = "random") / .N * ngroup)]]
        corr_dt[, pValue_group := g[ceiling(rank(-pValue, ties = "random") / .N * ngroup)]]
        corr_dt[, input_read_group := g[ceiling(rank(-input_read_corr, ties = "random") / .N * ngroup)]]
    }
    ## input metrics

    corr_dt[, input_frag_group := g[ceiling(rank(-input_frag_corr, ties = "random") / .N * ngroup)]]


    corr_dt$id = factor(corr_dt$id, levels = qgr$id)
    corr_dt = corr_dt[order(id)]

    todo_groups = colnames(corr_dt)[grepl("group", colnames(corr_dt))]
# browser()

    all_motif = lapply(todo_groups, function(grp){
        parallel::mclapply(levels(corr_dt[[grp]]), function(sel){
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


#' Title
#'
#' @param me_res
#' @param k
#'
#' @return
#' @export
#'
#' @examples
subset_MotifEnrichmentResults = function(me_res, k){
    ### top enriched can be cached but score and pval are group calcs
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


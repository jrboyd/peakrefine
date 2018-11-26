#' Assembles several metrics for a set of peaks from:
#' 1) the original peak call
#' 2) the ChIP-seq bam file
#' 3) an optional input (no pull-down) bam file
#'
#' Stratify by each metric and conduct motif enrichment for each.
#'
#' @param bam_file character. Path to treatment/ChIP-seq .bam file, must have index at .bam.bai.
#' @param bam_input character. Path to input/non-ChIP control .bam file, must have index at .bam.bai.
#' @param query_gr GRanges.  Regions to calculate cross correlation for.
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
score_motif = function(bam_file, bam_input, query_gr, fl, motif_res = NULL,
                       pwm = PWMLogn.hg19.MotifDb.Hsap,
                       include_fl_independent = FALSE,
                       out_dir = file.path("~/R/peakrefine/cache", sub(".bam", "", basename(bam_file))),
                       nbases = 200, ncores = 32){
    group = metric = id = NULL #declare data.table bindings
    registerCoresPWMEnrich(ncores)
    useBigMemoryPWMEnrich(TRUE)
    options(mc.cores = ncores)

    message("frag_len ", fl)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(out_dir)
    ## loading full motifs.
    cfile = paste0("cache_motif_res_", nbases, "nbases_", length(query_gr), "seq.save")
    if(is.null(motif_res)){
        if(file.exists(cfile)){
            message("loading cached full motif enrichment.")
            load(cfile)
        }else{
            message("calculating full motif enrichment...")
            sequences = getSeq(Hsapiens,
                               resize(query_gr, nbases, fix = "center"))
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
                                      query_gr,
                                      frag_len = fl,
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
                                           query_gr,
                                           frag_len = fl,
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
    corr_dt[, diff_group := g[ceiling(rank(-diff_corr, ties.method = "random") / .N * ngroup)]]
    corr_dt[, frag_group := g[ceiling(rank(-frag_corr, ties.method = "random") / .N * ngroup)]]

    ## fragment length independent metrics
    if(include_fl_independent){
        corr_dt[, input_count_group := g[ceiling(rank(-input_count, ties.method = "random") / .N * ngroup)]]
        corr_dt[, count_group := g[ceiling(rank(-count, ties.method = "random") / .N * ngroup)]]
        corr_dt[, read_group := g[ceiling(rank(-read_corr, ties.method = "random") / .N * ngroup)]]
        corr_dt[, signalValue_group := g[ceiling(rank(-signalValue, ties.method = "random") / .N * ngroup)]]
        corr_dt[, pValue_group := g[ceiling(rank(-pValue, ties.method = "random") / .N * ngroup)]]
        corr_dt[, input_read_group := g[ceiling(rank(-input_read_corr, ties.method = "random") / .N * ngroup)]]
    }
    ## input metrics

    corr_dt[, input_frag_group := g[ceiling(rank(-input_frag_corr, ties.method = "random") / .N * ngroup)]]


    corr_dt$id = factor(corr_dt$id, levels = query_gr$id)
    corr_dt = corr_dt[order(id)]

    todo_groups = colnames(corr_dt)[grepl("group", colnames(corr_dt))]

    all_motif = lapply(todo_groups, function(grp){
        parallel::mclapply(levels(corr_dt[[grp]]), function(sel){
            cfile = paste0("cache_", grp, "_", sel, "_", ngroup, "_", fl, ".save")
            if(file.exists(cfile)){
                load(cfile)
            }else{
                message("calculating motifs for ",
                        grp, " ", sel)
                # query_gr = GRanges(corr_dt[corr_dt[[grp]] == sel])
                k = which(corr_dt[[grp]] == sel)
                myMotifRes = groupReport(subset_MotifEnrichmentResults(motif_res, k))
                # myMotifRes = myMotif(query_gr, nbases)
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
#' @param corr_res
#' @param base_gr
#'
#' @return
#' @export
#'
#' @examples
scoreMetrics = function(corr_res, base_gr, ngroup = 8){
    rl = corr_res$read_length
    fl = corr_res$fragment_length
    read_corrs = corr_res$read_correlation

    if(is.null(corr_res$read_correlation$treatment)){
        flex_frag_corrs = corr_res$flex_fragment_correlation
        stable_frag_corrs = corr_res$stable_fragment_correlation
        peak_strand_corr = corr_res$full_correlation_results

        qgr = base_gr
        stopifnot(all(read_corrs$id == qgr$name))
        stopifnot(all(flex_frag_corrs$id == qgr$name))
        stopifnot(all(stable_frag_corrs$id == qgr$name))
        qgr$read_corr = read_corrs$correlation
        qgr$flex_frag_corr = flex_frag_corrs$correlation
        qgr$stable_frag_corr = stable_frag_corrs$correlation
        qgr$flex_frag_len = flex_frag_corrs$shift

        qdt = as.data.table(qgr)
        to_score = c("signalValue", "qValue", "stable_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len")

    }else{
        read_corrs = corr_res$read_correlation[treatment == "pulldown"]
        flex_frag_corrs = corr_res$flex_fragment_correlation[treatment == "pulldown"]
        stable_frag_corrs = corr_res$stable_fragment_correlation[treatment == "pulldown"]

        read_corrs.input = corr_res$read_correlation[treatment == "input"]
        flex_frag_corrs.input = corr_res$flex_fragment_correlation[treatment == "input"]
        stable_frag_corrs.input = corr_res$stable_fragment_correlation[treatment == "input"]

        qgr = base_gr
        stopifnot(all(read_corrs$id == qgr$name))
        stopifnot(all(flex_frag_corrs$id == qgr$name))
        stopifnot(all(stable_frag_corrs$id == qgr$name))

        qgr$read_corr = read_corrs$correlation
        qgr$flex_frag_corr = flex_frag_corrs$correlation
        qgr$stable_frag_corr = stable_frag_corrs$correlation
        qgr$flex_frag_len = flex_frag_corrs$shift

        qgr$read_corr_input = read_corrs.input$correlation
        qgr$flex_frag_corr_input = flex_frag_corrs.input$correlation
        qgr$stable_frag_corr_input = stable_frag_corrs.input$correlation
        qgr$flex_frag_len_input = flex_frag_corrs.input$shift

        qdt = as.data.table(qgr)
        to_score = c("signalValue", "qValue", "stable_frag_corr",
                     "flex_frag_corr", "read_corr", "flex_frag_len",
                     "read_corr_input", "flex_frag_corr_input",
                     "stable_frag_corr_input", "flex_frag_len_input")
    }

    stopifnot(all(to_score %in% colnames(qdt)))
    # for(ts in to_score){
    #     if(any(is.nan(qdt[[ts]])))
    #         qdt[is.nan(get(ts)), ][[ts]] = 0
    # }
    #
    # g = factor(paste0("g", seq_len(ngroup)),
    #            levels = rev(paste0("g", seq_len(ngroup))))
    #
    # for(ts in to_score){
    #     set.seed(seed)
    #     qdt[, paste0(ts, "_group") := g[ceiling(rank(-get(ts), ties.method = "random") / .N * ngroup)]]
    # }
    for(ts in to_score){
        group_var(qdt, ts, 0, 8)
    }

    qdt
}


group_var = function(qdt, ts, nan_val = 0, ngroup = 8, seed = 0){
    # for(ts in to_score){
    if(any(is.nan(qdt[[ts]])))
        qdt[is.nan(get(ts)), ][[ts]] = nan_val
    # }

    g = factor(paste0("g", seq_len(ngroup)),
               levels = rev(paste0("g", seq_len(ngroup))))

    # for(ts in to_score){
    set.seed(seed)
    qdt[, paste0(ts, "_group") := g[ceiling(rank(-get(ts), ties.method = "random") / .N * ngroup)]]
    # }
    qdt
}

#' Title
#'
#' @param corr_res
#' @param base_gr
#' @param pwm
#' @param seq
#' @param nbases
#' @param cache_path
#' @param cach_version
#'
#' @return
#' @export
#'
#' @examples
calcMotifEnrichment = function(corr_res, base_gr, qdt, todo_groups, pwm, seq, bam_md5, qgr_md5 = NULL, nbases = 200, ncores = 8,
                               cache_path = "~/.cache_peakrefine", force_overwrite = FALSE, force_overwrite_pre = FALSE,
                               cach_version = "v1"){
    base_gr$name = factor(base_gr$name, levels = qdt$name)
    base_gr = base_gr[order(base_gr$name)]
    stopifnot(all(base_gr$name == qdt$name))
    if(is.null(qgr_md5)){
        qgr_md5 = digest::digest(base_gr)
    }
    options(mc.cores = ncores)
    bfc_motif = BiocFileCache::BiocFileCache(cache_path, ask = FALSE)
    motif_key = paste(qgr_md5, digest::digest(pwm),
                      digest::digest(seq), nbases, cach_version, sep = "_")
    motif_res = bfcif(bfc_motif, motif_key, function(){
        message("cached results not found, gathering motif info.")
        pre_motif(base_gr, pwm = pwm, seq_reference = seq, nbases = nbases, ncores = ncores)
    }, force_overwrite = force_overwrite_pre)
    motif_md5 = digest::digest(motif_res)
    all_motif = lapply(todo_groups, function(grp){
        motif_group_key = paste(bam_md5, digest::digest(base_gr),
                                motif_md5,
                                digest::digest(qdt[[grp]]), grp,
                                cach_version, nbases, sep = "_")
        bfcif(bfc_motif, motif_group_key, function(){
            message("calculating motifs for ", grp)
            grp_motif = parallel::mclapply(levels(qdt[[grp]]), mc.preschedule = FALSE, function(sel){
                k = which(as.character(qdt[[grp]]) == sel)
                groupReport(subset_MotifEnrichmentResults(motif_res, k, pwm))
            })
            names(grp_motif) = levels(qdt[[grp]])
            rbindlist(lapply(grp_motif, function(y){
                dt = as.data.table(as.data.frame(y))
                colnames(dt) = gsub("\\.", "_", colnames(dt))
                dt
            }), use.names = TRUE, idcol = "group")
        }, force_overwrite = force_overwrite)


    })
    names(all_motif) = todo_groups
    dt_motif = rbindlist(all_motif, use.names = TRUE, idcol = "metric")
    lev = unique(dt_motif$group)
    lev = lev[order(as.numeric(gsub("[A-Za-z]", "", lev)), decreasing = TRUE)]
    dt_motif$group = factor(dt_motif$group, levels = lev)
    dt_motif = dt_motif[order(group)][order(metric)][order(id)]
    dt_motif
}

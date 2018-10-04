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
    # browser()
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
# browser()

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

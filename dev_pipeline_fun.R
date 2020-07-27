library(Rsamtools)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(peakrefine)
library(digest)
library(BiocFileCache)

pipeline = function(bam_file, qgr, cach_version, inputs_file = NULL,
                    output_dir = "results_default",
                    frag_min = 50, ncores = 16, gen = "hg38",
                    to_score = c("signalValue", "qValue", "stable_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len"),
                    frag_max = 250,
                    pdfname = file.path(ouptut_dir, paste0("motif_res_", sub(".bam", "", basename(bam_file)), ".pdf")),
                    skip_motif = FALSE,
                    force_overwrite_motif = FALSE,
                    force_overwrite_corr = FALSE){
    # if(file.exists(pdfname)) return()
    # prep query granges

    if(is.null(qgr$signalValue)){
        qgr$signalValue = 1
    }
    if(is.null(qgr$qValue)){
        qgr$qValue = 1
    }
    subset(qgr, !grepl("_", seqnames(qgr)))
    qgr = qgr[order(qgr$signalValue, decreasing = TRUE)][1:min(30000, length(qgr)) ]
    qgr$name = paste0("peak_", seq_along(qgr))
    qgr = dropSeqlevels(qgr, "chrU13369.1", pruning.mode="coarse")
    qgr = sort(qgr)
    qgr = resize(qgr, 800, fix = "center")
    qgr = harmonize_seqlengths(qgr, bam_file)
    base_gr = qgr
    # precalc md5
    bam_md5 = tools::md5sum(path.expand(bam_file))
    qgr_md5 = digest::digest(base_gr)
    bfc_motif = BiocFileCache::BiocFileCache("~/.cache_motif", ask = FALSE)
    options(mc.cores = ncores)
    getOption("mc.cores", 1L)
    if(is.null(inputs_file)){
        # calc correlation
        corr_res = calcCorrMetrics(bam_file, qgr, frag_min, frag_max,
                                   cache_path = "~/.cache_corr/",
                                   cach_version = cach_version,
                                   bam_md5 = bam_md5, qgr_md5 = qgr_md5, force_overwrite = force_overwrite_corr)

        uid = unique(corr_res$full_correlation_results$id)
        tp = sample(uid, min(500, length(uid)))

        p1 = ggplot(corr_res$full_correlation_results[id %in% tp],
                    aes(x = shift, y = correlation, group = id)) +
            geom_path(alpha = .05, size = 3)

        p2 = ggplot(corr_res$full_correlation_results[id %in% tp],
                    aes(x = shift, y = correlation)) +
            stat_density_2d(geom = "raster", aes(fill = stat(density)), contour = FALSE)
    }else{
        input_md5 = tools::md5sum(path.expand(inputs_file))
        corr_res = calcCorrMetrics(bam_file, qgr, frag_min, frag_max,
                                   cache_path = "~/.cache_corr/",
                                   cach_version = cach_version,
                                   bam_md5 = bam_md5, qgr_md5 = qgr_md5, force_overwrite = force_overwrite_corr)
        input_res = calcCorrMetrics(inputs_file, qgr, frag_min, frag_max,
                                    cache_path = "~/.cache_corr/",
                                    cach_version = cach_version,
                                    bam_md5 = input_md5, qgr_md5 = qgr_md5, force_overwrite = force_overwrite_corr)

        corr_res$full_correlation_results$treatment = "pulldown"
        input_res$full_correlation_results$treatment = "input"
        corr_res$full_correlation_results = rbind(corr_res$full_correlation_results, input_res$full_correlation_results)

        corr_res$stable_fragment_correlation$treatment = "pulldown"
        input_res$stable_fragment_correlation = input_res$full_correlation_results[shift == corr_res$fragment_length]
        corr_res$stable_fragment_correlation = rbind(corr_res$stable_fragment_correlation, input_res$stable_fragment_correlation)

        corr_res$flex_fragment_correlation$treatment = "pulldown"
        input_res$flex_fragment_correlation$treatment = "input"
        corr_res$flex_fragment_correlation = rbind(corr_res$flex_fragment_correlation, input_res$flex_fragment_correlation)

        corr_res$read_correlation$treatment = "pulldown"
        input_res$read_correlation$treatment = "input"
        corr_res$read_correlation = rbind(corr_res$read_correlation, input_res$read_correlation)

        remove("input_res")

        uid = unique(corr_res$full_correlation_results$id)
        tp = sample(uid, min(500, length(uid)))

        p1 = ggplot(corr_res$full_correlation_results[id %in% tp],
                    aes(x = shift, y = correlation, group = id, color = treatment)) +
            geom_path(alpha = .05, size = 3) +
            facet_wrap("treatment", ncol = 1)

        p2 = ggplot(corr_res$full_correlation_results[id %in% tp],
                    aes(x = shift, y = correlation)) +
            stat_density_2d(geom = "raster", aes(fill = stat(density)), contour = FALSE) +
            facet_wrap("treatment", ncol = 1)
    }




    format(object.size(corr_res$full_correlation_results), units = "GB")

    qdt = scoreMetrics(corr_res, base_gr)

    tmp = qgr
    seqlengths(tmp) = NA
    names(tmp) = tmp$name
    bdt = seqsetvis::ssvFetchBam(bam_file, tmp, fragLens = 1, return_data.table = TRUE, target_strand = "both")
    bdt = bdt[, .(total_reads = sum(y)), .(id, strand)]
    bdt = dcast(bdt, id ~ strand, value.var = "total_reads")
    colnames(bdt)[colnames(bdt) == "+"] = "pos"
    colnames(bdt)[colnames(bdt) == "-"] = "neg"
    bdt[, total := pos + neg ]
    bdt[, pos_fraction := (pos ) / (total ) ]
    bdt[is.nan(pos_fraction) & pos > 0, pos_fraction := 1]
    bdt[is.nan(pos_fraction) & pos == 0, pos_fraction := 0]

    pulldown_dt = bdt

    bdt = seqsetvis::ssvFetchBam(inputs_file, tmp, fragLens = 1, return_data.table = TRUE, target_strand = "both")
    bdt = bdt[, .(total_reads = sum(y)), .(id, strand)]
    bdt = dcast(bdt, id ~ strand, value.var = "total_reads")
    colnames(bdt)[colnames(bdt) == "+"] = "pos"
    colnames(bdt)[colnames(bdt) == "-"] = "neg"
    bdt[, total := pos + neg ]
    bdt[, pos_fraction := (pos ) / (total ) ]
    bdt[is.nan(pos_fraction) & pos > 0, pos_fraction := 1]
    bdt[is.nan(pos_fraction) & pos == 0, pos_fraction := 0]

    input_dt = bdt
    colnames(input_dt)[-1] = paste0(colnames(input_dt)[-1], "_input")

    qdt = merge(qdt, pulldown_dt, by.x = "name", by.y = "id")
    qdt = merge(qdt, input_dt, by.x = "name", by.y = "id")

    if(skip_motif) return(list(qdt = qdt, corr_res = corr_res))

    if(gen == "hg38"){
        library(PWMEnrich.Hsapiens.background)
        data(PWMLogn.hg19.MotifDb.Hsap)
        library(BSgenome.Hsapiens.UCSC.hg38)
        pwm = PWMLogn.hg19.MotifDb.Hsap
        seq = Hsapiens
    }else if(gen == "mm10"){
        library(PWMEnrich.Mmusculus.background)
        data(PWMLogn.mm9.MotifDb.Mmus)
        library(BSgenome.Mmusculus.UCSC.mm10)
        pwm = PWMLogn.mm9.MotifDb.Mmus
        seq = Mmusculus
    }
    nbases = 200
    qdt$name = factor(qdt$name, levels = base_gr$name)
    qdt = qdt[order(name)]
    todo_groups = colnames(qdt)[grepl("group", colnames(qdt))]
    dt_motif = calcMotifEnrichment(corr_res, base_gr, qdt, todo_groups,
                                   pwm, seq, bam_md5 = bam_md5,
                                   cache_path = "~/.cache_motif/",
                                   cach_version = cach_version, ncores = ncores, force_overwrite_pre = force_overwrite_motif, force_overwrite = force_overwrite_motif)

    library("GGally")
    p2a = ggpairs(qdt[sample(.N, min(500, nrow(qdt))), .(signalValue, qValue, stable_frag_corr,
                                                         flex_frag_corr, read_corr, flex_frag_len)],
                  lower = list(continuous = "density"), upper = list(continuous = "density"))


    metric_dt = melt(qdt[, c("name", to_score), with = F], id.vars = "name",
                     value.name = "score", variable.name = "metric")
    metric_dt[, score_rank := rank(-score, ties.method = "first"), by = .(metric)]
    metric_dt[, score_norm := (score - min(score)) / (max(score) - min(score)), by = .(metric) ]
    metric_dt[, score_rank_norm := (score_rank - min(score_rank)) / (max(score_rank) - min(score_rank)), by = .(metric) ]
    name_tp = sample(unique(metric_dt$name), min(500, length(unique(metric_dt$name))))

    p3 = ggplot(metric_dt[name %in% name_tp], aes(x = score, y = score_rank_norm)) +
        geom_point() +
        facet_wrap("metric", scales = "free_x")

    id_oi = dt_motif[, .(score_range = quantile(top_motif_prop, .75)),
                     by = .(id)][order(score_range, decreasing = TRUE)]$id[1:3]
    p4 = ggplot(dt_motif[id %in% id_oi],
                aes(x = group,
                    y = top_motif_prop,
                    color = id,
                    group =  paste(id, metric))) +
        geom_path() +
        geom_point() +
        facet_wrap("metric", ncol = 3) +
        theme(panel.grid.major.y = element_line()) +
        theme_classic() +
        theme(strip.text.y = element_text(angle = 0))

    p5 = ggplot(dt_motif[id %in% id_oi],
                aes(x = group,
                    y = raw_score,
                    color = id,
                    group =  paste(id, metric))) +
        geom_path() +
        geom_point() +
        facet_wrap("metric", ncol = 3) +
        theme(panel.grid.major.y = element_line()) +
        theme_classic() +
        theme(strip.text.y = element_text(angle = 0))


    qdt[, diff_corr := stable_frag_corr - read_corr]
    # qdt$combined_passing = FALSE
    # qdt[signalValue > 5 & flex_frag_corr > .75 & read_corr < .6 &
    #         flex_frag_len > 50 & flex_frag_len <= 150 & diff_corr > .1,
    #     combined_passing := TRUE ]

    #add combined_passing
    qdt$combined_passing = TRUE
    qdt[read_corr > .8, combined_passing := FALSE]
    qdt[stable_frag_corr - read_corr < .1, combined_passing := FALSE]
    qdt[flex_frag_len < 100, combined_passing := FALSE]
    qdt[flex_frag_len > 300, combined_passing := FALSE]
    todo_groups = c("combined_passing", todo_groups)



    npass = sum(qdt$combined_passing)
    npass
    nrow(qdt)

    todo_groups_direction = rep(1, length(todo_groups))
    names(todo_groups_direction) = todo_groups
    todo_groups_direction["read_corr_group"] = -1

    for(group in todo_groups){
        qdt[[sub("_group", "_passing", group)]] =
            factor(qdt[, rank(-1*todo_groups_direction[group]*get(sub("_group", "", group)), ties.method = "first") <= npass])
    }

    todo_passing = colnames(qdt)[grepl("passing", colnames(qdt))]
    dt_motif_pass = calcMotifEnrichment(corr_res, base_gr, qdt, todo_passing,
                                        pwm, seq,
                                        bam_md5 = bam_md5,
                                        cache_path = "~/.cache_motif/",
                                        cach_version = cach_version, force_overwrite = force_overwrite_motif)

    id_oi_pass = dt_motif_pass[, .(score_range = quantile(top_motif_prop, .75)), by = .(id)][order(score_range, decreasing = TRUE)]$id[1]
    p6 = ggplot(dt_motif_pass[id %in% id_oi_pass],
                aes(x = group,
                    y = top_motif_prop,
                    fill = group,
                    group =  paste(id, metric))) +
        geom_bar(stat = "identity") +
        # geom_point() +
        facet_wrap("metric", ncol = 3) +
        theme(panel.grid.major.y = element_line()) +
        theme_classic() +
        theme(strip.text.y = element_text(angle = 0))

    # ggplot(dt_motif_pass[id %in% id_oi_pass],
    #        aes(x = metric,
    #            y = top_motif_prop,
    #            fill = group,
    #            group =  paste(group, id, metric))) +
    #     geom_bar(stat = "identity") +
    #     # geom_point() +
    #     facet_wrap("group", ncol = 3) +
    #     theme(panel.grid.major.y = element_line()) +
    #     theme_classic() +
    #     theme(strip.text.y = element_text(angle = 0))

    p7 = ggplot(dt_motif_pass[id %in% id_oi_pass],
                aes(x = group,
                    y = raw_score,
                    fill = group,
                    group =  paste(id, metric))) +
        geom_bar(stat = "identity") +
        # geom_point() +
        facet_wrap("metric", ncol = 3) +
        theme(panel.grid.major.y = element_line()) +
        theme_classic() +
        theme(strip.text.y = element_text(angle = 0))

    p8 = ggplot(dt_motif_pass[id %in% id_oi_pass],
                aes(x = metric,
                    y = raw_score,
                    fill = group,
                    group =  paste(group, id, metric))) +
        geom_bar(stat = "identity") +
        # geom_point() +
        facet_wrap("group", ncol = 3) +
        theme(panel.grid.major.y = element_line()) +
        theme_classic() +
        theme(strip.text.y = element_text(angle = 0), axis.text.x = element_text(angle = 90))

    save(qdt, corr_res, dt_motif, dt_motif_pass, file = paste0(sub(".pdf", "", pdfname), "_data.save"))
    plist = list(p1, p2, p2a, p3, p4, p5, p6, p7, p8)
    pdf(pdfname)
    lapply(plist, print)
    dev.off()
    return(corr_res)
}

get_corr_res = function(bam_file, qgr, cach_version, force_overwrite = FALSE, frag_min = 50, ncores = 16,
                        frag_max = 250, pdfname = paste0("motif_res_", sub(".bam", "", basename(bam_file)), ".pdf")){
    # if(file.exists(pdfname)) return()
    # prep query granges
    subset(qgr, !grepl("_", seqnames(qgr)))
    qgr = qgr[order(qgr$signalValue, decreasing = TRUE)][1:min(30000, length(qgr)) ]
    qgr$name = paste0("peak_", seq_along(qgr))
    qgr = dropSeqlevels(qgr, "chrU13369.1", pruning.mode="coarse")
    qgr = sort(qgr)
    qgr = resize(qgr, 800, fix = "center")
    qgr = harmonize_seqlengths(qgr, bam_file)
    base_gr = qgr
    # precalc md5
    bam_md5 = tools::md5sum(path.expand(bam_file))
    qgr_md5 = digest::digest(base_gr)
    bfc_motif = BiocFileCache::BiocFileCache("~/.cache_motif", ask = FALSE)
    options(mc.cores = 16)
    getOption("mc.cores", 1L)

    # calc correlation
    corr_res = calcCorrMetrics(bam_file, qgr, frag_min, frag_max,
                               cache_path = "~/.cache_corr/",
                               cach_version = cach_version,
                               bam_md5 = bam_md5, qgr_md5 = qgr_md5, force_overwrite = force_overwrite)

    corr_res
}

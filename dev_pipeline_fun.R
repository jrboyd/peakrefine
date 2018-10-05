library(Rsamtools)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(peakrefine)
library(digest)
library(BiocFileCache)

pipeline = function(bam_file, qgr, cach_version, frag_min = 50,
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
    options(mc.cores = 32)
    getOption("mc.cores", 1L)

    # calc correlation
    corr_res = calcCorrMetrics(bam_file, qgr, frag_min, frag_max,
                               cache_path = "~/.cache_corr/",
                               cach_version = cach_version,
                               bam_md5 = bam_md5, qgr_md5 = qgr_md5)

    p1 = ggplot(corr_res$full_correlation_results[id %in% unique(id)[1:500]],
                aes(x = shift, y = correlation, group = id)) +
        geom_path(alpha = .05, size = 3)

    p2 = ggplot(corr_res$full_correlation_results[id %in% unique(id)[1:500]],
                aes(x = shift, y = correlation)) +
        stat_density_2d(geom = "raster", aes(fill = stat(density)), contour = FALSE)

    format(object.size(corr_res$full_correlation_results), units = "GB")

    library(PWMEnrich.Hsapiens.background)
    data(PWMLogn.hg19.MotifDb.Hsap)
    library(BSgenome.Hsapiens.UCSC.hg38)
    nbases = 200
    pwm = PWMLogn.hg19.MotifDb.Hsap
    seq = Hsapiens
    qdt = scoreMetrics(corr_res, base_gr)
    todo_groups = colnames(qdt)[grepl("group", colnames(qdt))]

    dt_motif = calcMotifEnrichment(corr_res, base_gr, qdt, todo_groups,
                                   pwm, seq, bam_md5 = bam_md5, qgr_md5 = qgr_md5,
                                   cache_path = "~/.cache_motif/",
                                   cach_version = cach_version)

    library("GGally")
    p2a = ggpairs(qdt[sample(.N, min(500, nrow(qdt))), .(signalValue, qValue, mean_frag_corr,
                                   flex_frag_corr, read_corr, flex_frag_len)],
            lower = list(continuous = "density"))

    to_score = c("signalValue", "qValue", "mean_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len")
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


    qdt[, diff_corr := mean_frag_corr - read_corr]
    qdt$combined_passing = FALSE
    qdt[signalValue > 5 & flex_frag_corr > .75 & read_corr < .6 &
            flex_frag_len > 50 & flex_frag_len <= 150 & diff_corr > .1,
        combined_passing := TRUE ]

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
                                        bam_md5 = bam_md5, qgr_md5 = qgr_md5,
                                        cache_path = "~/.cache_motif/",
                                        cach_version = cach_version)

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

    plist = list(p1, p2, p2a, p3, p4, p5, p6, p7, p8)
    pdf(pdfname)
    lapply(plist, print)
    dev.off()
    return(corr_res)
}

get_corr_res = function(bam_file, qgr, cach_version, frag_min = 50,
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
    options(mc.cores = 32)
    getOption("mc.cores", 1L)

    # calc correlation
    corr_res = calcCorrMetrics(bam_file, qgr, frag_min, frag_max,
                               cache_path = "~/.cache_corr/",
                               cach_version = cach_version,
                               bam_md5 = bam_md5, qgr_md5 = qgr_md5)

    corr_res
}

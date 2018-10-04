library(Rsamtools)
library(data.table)
library(ggplot2)
library(GenomicRanges)
library(peakrefine)
library(digest)
library(BiocFileCache)

cach_version = "v5"

# bam_file = "~/ENCODE_EGR1/K562_EGR1_myers_rep1_ENCFF645IYJ.bam"
bam_file = "~/ENCODE_EGR1/K562_EGR1_myers_rep2_ENCFF021CMK.bam"
# qgr = rtracklayer::import("~/ENCODE_EGR1/K562_EGR1_myers_rep1_ENCFF178PTI.bed", format = "narrowPeak")
qgr = rtracklayer::import("~/ENCODE_EGR1/K562_EGR1_myers_rep2_ENCFF913QDG.bed", format = "narrowPeak")
subset(qgr, !grepl("_", seqnames(qgr)))
qgr = qgr[order(qgr$signalValue, decreasing = TRUE)][1:30000]
qgr$name = paste0("peak_", seq_along(qgr))

qgr = dropSeqlevels(qgr, "chrU13369.1", pruning.mode="coarse")
qgr = sort(qgr)
qgr = resize(qgr, 800, fix = "center")

qgr = harmonize_seqlengths(qgr, bam_file)

base_gr = qgr


#temp
if(F){
    bw_dt = ssvRecipes::myFetchStrandedBam(bam_file, GRanges(gsub(",", "", "chr6:36,153,642-36,155,265")))
    ggplot(bw_dt, aes(x = x, y = y, color = strand)) + geom_path()
}

bam_md5 = tools::md5sum(bam_file)

frag_min = 50
frag_max = 150


bfc_motif = BiocFileCache::BiocFileCache("~/.cache_motif", ask = FALSE)

options(mc.cores = 32)



getOption("mc.cores", 1L)



##TO DO force overwrite
force_overwrite = TRUE

corr_res = calcCorrMetrics(bam_file, qgr, frag_min, frag_max,
                           cache_path = "~/.cache_corr/",
                           cach_version = cach_version, bam_md5 = bam_md5)

# rl = corr_res$read_length
# fl = corr_res$fragment_length
# read_corrs = corr_res$read_correlation
# flex_frag_corrs = corr_res$flex_fragment_correlation
# mean_frag_corrs = corr_res$mean_fragment_correlation
# peak_strand_corr = corr_res$full_correlation_results

ggplot(peak_strand_corr[id %in% unique(id)[1:500]],
       aes(x = shift, y = correlation, group = id)) +
    geom_path(alpha = .05, size = 3)

ggplot(peak_strand_corr[id %in% unique(id)[1:500]],
       aes(x = shift, y = correlation)) +
    stat_density_2d(geom = "raster", aes(fill = stat(density)), contour = FALSE)

# df = dcast(peak_strand_corr, id~shift, value.var = "correlation")
# mat = as.matrix(df[, -1])
# rownames(mat) = df[[1]]
# hclust(mat, dist)
# gplots::heatmap.2(mat[1:100,], trace = "n", Colv = FALSE, labRow = "", labCol = "", dendrogram = "none",
#                   col = rgb(colorRamp(c("blue", "white", "red"))(0:50/50)/255))

# read_corrs
# flex_frag_corrs
# mean_frag_corrs
format(object.size(peak_strand_corr), units = "GB")

library(PWMEnrich.Hsapiens.background)
data(PWMLogn.hg19.MotifDb.Hsap)
library(BSgenome.Hsapiens.UCSC.hg38)
nbases = 200
pwm = PWMLogn.hg19.MotifDb.Hsap
seq = Hsapiens

calcMotifEnrichment = function(corr_res, base_gr, pwm, seq, nbases, cache_path, cach_version){
    motif_key = paste(digest::digest(base_gr), digest::digest(pwm),
                      digest::digest(seq), nbases, cach_version, sep = "_")
    motif_res = bfcif(bfc_motif, motif_key, function(){
        message("cached results not found, gathering motif info.")
        pre_motif(base_gr, pwm = pwm, seq_reference = Hsapiens, nbases = nbases)
    })
    motif_md5 = digest::digest(motif_res)

    rl = corr_res$read_length
    fl = corr_res$fragment_length
    read_corrs = corr_res$read_correlation
    flex_frag_corrs = corr_res$flex_fragment_correlation
    mean_frag_corrs = corr_res$mean_fragment_correlation
    peak_strand_corr = corr_res$full_correlation_results

    qgr = base_gr
    stopifnot(all(read_corrs$id == qgr$name))
    stopifnot(all(flex_frag_corrs$id == qgr$name))
    stopifnot(all(mean_frag_corrs$id == qgr$name))
    qgr$read_corr = read_corrs$correlation
    qgr$flex_frag_corr = flex_frag_corrs$correlation
    qgr$mean_frag_corr = mean_frag_corrs$correlation
    qgr$flex_frag_len = flex_frag_corrs$shift

    qdt = as.data.table(qgr)
    to_score = c("signalValue", "qValue", "mean_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len")
    stopifnot(all(to_score %in% colnames(qdt)))

    ngroup = 8
    g = factor(paste0("g", seq_len(ngroup)),
               levels = paste0("g", seq_len(ngroup)))

    for(ts in to_score){
        qdt[, paste0(ts, "_group") := g[ceiling(rank(-get(ts), ties.method = "random") / .N * ngroup)]]
    }

    todo_groups = colnames(qdt)[grepl("group", colnames(qdt))]

    todo_dt = rbindlist(lapply(todo_groups, function(grp){
        data.table(sel = levels(qdt[[grp]]), grp = grp)
    }))

    all_motif = parallel::mclapply(seq_len(nrow(todo_dt)), mc.preschedule = FALSE, function(i){
        Sys.sleep(i / 2) #sql gets overwhelmed without staggering out a bit
        grp = todo_dt[i, ]$grp
        sel = todo_dt[i, ]$sel
        motif_group_key = paste(bam_md5, digest::digest(base_gr),
                                motif_md5,
                                digest::digest(qdt[[grp]]), grp,
                                sel, cach_version, nbases, sep = "_")
        myMotifRes = bfcif(bfc_motif, motif_group_key, function(){
            message("calculating motifs for ",
                    grp, " ", sel)
            k = which(qdt[[grp]] == sel)
            groupReport(subset_MotifEnrichmentResults(motif_res, k, pwm))
        })
        myMotifRes
    })
}


# motif_key = paste(digest::digest(base_gr), digest::digest(pwm), digest::digest(seq), nbases, cach_version, sep = "_")
#
# motif_res = bfcif(bfc_motif, motif_key, function(){
#     message("cached results not found, gathering motif info.")
#     pre_motif(qgr, pwm = pwm, seq_reference = Hsapiens, nbases = nbases)
# })
#
# motif_md5 = digest::digest(motif_res)

# stopifnot(all(read_corrs$id == qgr$name))
# stopifnot(all(flex_frag_corrs$id == qgr$name))
# stopifnot(all(mean_frag_corrs$id == qgr$name))
# qgr$read_corr = read_corrs$correlation
# qgr$flex_frag_corr = flex_frag_corrs$correlation
# qgr$mean_frag_corr = mean_frag_corrs$correlation
# qgr$flex_frag_len = flex_frag_corrs$shift



# qdt = as.data.table(qgr)
# to_score = c("signalValue", "qValue", "mean_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len")
# stopifnot(all(to_score %in% colnames(qdt)))

# ngroup = 8
# g = factor(paste0("g", seq_len(ngroup)),
#            levels = paste0("g", seq_len(ngroup)))
#
# for(ts in to_score){
#     qdt[, paste0(ts, "_group") := g[ceiling(rank(-get(ts), ties.method = "random") / .N * ngroup)]]
# }

# qdt[, rank1 := rank(-qValue, ties.method = "first")]
# qdt[, rank2 := rank(-flex_frag_corr, ties.method = "first")]
# qdt[, rankMax := max(rank1, rank2), by = .(name)]
# qdt[, paste0("combined", "_group") := g[ceiling(rank(rankMax, ties.method = "random") / .N * ngroup)]]

# todo_groups = colnames(qdt)[grepl("group", colnames(qdt))]
#
# todo_dt = rbindlist(lapply(todo_groups, function(grp){
#     data.table(sel = levels(qdt[[grp]]), grp = grp)
# }))
#
# all_motif = parallel::mclapply(seq_len(nrow(todo_dt)), mc.preschedule = FALSE, function(i){
#     Sys.sleep(i / 2) #sql gets overwhelmed without staggering out a bit
#     grp = todo_dt[i, ]$grp
#     sel = todo_dt[i, ]$sel
#     motif_group_key = paste(bam_md5, digest::digest(base_gr),
#                             motif_md5,
#                             digest::digest(qdt[[grp]]), grp,
#                             sel, cach_version, nbases, sep = "_")
#     myMotifRes = bfcif(bfc_motif, motif_group_key, function(){
#         message("calculating motifs for ",
#                 grp, " ", sel)
#         k = which(qdt[[grp]] == sel)
#         groupReport(subset_MotifEnrichmentResults(motif_res, k, pwm))
#     })
#     myMotifRes
# })


names(all_motif) = paste(todo_dt$sel, todo_dt$grp)

dt_motif = rbindlist(lapply(all_motif, function(y){
    dt = as.data.table(as.data.frame(y))
    colnames(dt) = gsub("\\.", "_", colnames(dt))
    dt
}), use.names = TRUE, idcol = "name")

dt_motif[, c("group", "metric") := tstrsplit(name, " ")]

dt_motif$group = factor(dt_motif$group, levels = rev(levels(g)))
dt_motif = dt_motif[order(group)][order(metric)][order(id)]

library("GGally")
ggpairs(qdt[sample(.N, 500), .(signalValue, qValue, mean_frag_corr, flex_frag_corr, read_corr, flex_frag_len)], lower = list(continuous = "density"))


metric_dt = melt(qdt[, c("name", to_score), with = F], id.vars = "name", value.name = "score", variable.name = "metric")
metric_dt[, score_rank := rank(-score, ties.method = "first"), by = .(metric)]
metric_dt[, score_norm := (score - min(score)) / (max(score) - min(score)), by = .(metric) ]
metric_dt[, score_rank_norm := (score_rank - min(score_rank)) / (max(score_rank) - min(score_rank)), by = .(metric) ]

ggplot(metric_dt, aes(x = score, y = score_rank_norm)) + geom_point() + facet_wrap("metric", scales = "free_x")

id_oi = names(motif_res@res$pwms)
id_oi = id_oi[grepl("RUNX", id_oi)]

id_oi = dt_motif[, .(score_range = quantile(top_motif_prop, .75)), by = .(id)][order(score_range, decreasing = TRUE)]$id[1:3]
ggplot(dt_motif[id %in% id_oi],
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

ggplot(dt_motif[id %in% id_oi],
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
# qdt$passLow = FALSE
qdt$combined_passing = FALSE
qdt[signalValue > 5 & flex_frag_corr > .75 & read_corr < .6 &
        flex_frag_len > 50 & flex_frag_len <= 150 & diff_corr > .1, combined_passing := TRUE ]


# qdt[signalValue > 10 & qValue > 10 & flex_frag_corr > .65 & read_corr < .6 &
# flex_frag_len >= 100 & flex_frag_len <= 300, passHigh := TRUE ]

npass = sum(qdt$combined_passing)
npass
nrow(qdt)

todo_groups = colnames(qdt)[grepl("group", colnames(qdt))]
todo_groups_direction = rep(1, length(todo_groups))
names(todo_groups_direction) = todo_groups
todo_groups_direction["read_corr_group"] = -1

for(group in todo_groups){
    qdt[[sub("_group", "_passing", group)]] =
        qdt[, rank(-1*todo_groups_direction[group]*get(sub("_group", "", group)), ties.method = "first") <= npass]
}


todo_passing = colnames(qdt)[grepl("passing", colnames(qdt))]

todo_dt = rbindlist(lapply(todo_passing, function(grp){
    data.table(sel = c(FALSE, TRUE), grp = grp)
}))

pass_motif = parallel::mclapply(seq_len(nrow(todo_dt)), mc.preschedule = FALSE, function(i){
    Sys.sleep(i / 2) #sql gets overwhelmed without staggering out a bit
    grp = todo_dt[i, ]$grp
    sel = todo_dt[i, ]$sel
    motif_group_key = paste(bam_md5, digest::digest(base_gr),
                            motif_md5,
                            digest::digest(qdt[[grp]]), grp,
                            sel, cach_version, nbases, sep = "_")
    myMotifRes = bfcif(bfc_motif, motif_group_key, function(){
        message("calculating motifs for ",
                grp, " ", sel)
        k = which(qdt[[grp]] == sel)
        groupReport(subset_MotifEnrichmentResults(motif_res, k, pwm))
    }, force_overwrite = FALSE)
    myMotifRes
})

names(pass_motif) = paste(todo_dt$sel, todo_dt$grp)

dt_motif_pass = rbindlist(lapply(pass_motif, function(y){
    dt = as.data.table(as.data.frame(y))
    colnames(dt) = gsub("\\.", "_", colnames(dt))
    dt
}), use.names = TRUE, idcol = "name")

dt_motif_pass[, c("group", "metric") := tstrsplit(name, " ")]

dt_motif_pass$group = factor(dt_motif_pass$group, levels = c(TRUE, FALSE))
dt_motif_pass = dt_motif_pass[order(group)][order(metric)][order(id)]


id_oi = dt_motif_pass[, .(score_range = quantile(top_motif_prop, .75)), by = .(id)][order(score_range, decreasing = TRUE)]$id[1]
ggplot(dt_motif_pass[id %in% id_oi],
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

# ggplot(dt_motif_pass[id %in% id_oi],
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

ggplot(dt_motif_pass[id %in% id_oi],
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

ggplot(dt_motif_pass[id %in% id_oi],
       aes(x = metric,
           y = raw_score,
           fill = group,
           group =  paste(group, id, metric))) +
    geom_bar(stat = "identity") +
    # geom_point() +
    facet_wrap("group", ncol = 3) +
    theme(panel.grid.major.y = element_line()) +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 0))

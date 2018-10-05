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

bam_md5 = tools::md5sum(path.expand(bam_file))

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

ggplot(corr_res$full_correlation_results[id %in% unique(id)[1:500]],
       aes(x = shift, y = correlation, group = id)) +
    geom_path(alpha = .05, size = 3)

ggplot(corr_res$full_correlation_results[id %in% unique(id)[1:500]],
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
                               pwm, seq, bam_md5 = bam_md5,
                               cache_path = "~/.cache_motif/",
                               cach_version = cach_version)

library("GGally")
ggpairs(qdt[sample(.N, 500), .(signalValue, qValue, mean_frag_corr,
                               flex_frag_corr, read_corr, flex_frag_len)],
        lower = list(continuous = "density"))

to_score = c("signalValue", "qValue", "mean_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len")
metric_dt = melt(qdt[, c("name", to_score), with = F], id.vars = "name",
                 value.name = "score", variable.name = "metric")
metric_dt[, score_rank := rank(-score, ties.method = "first"), by = .(metric)]
metric_dt[, score_norm := (score - min(score)) / (max(score) - min(score)), by = .(metric) ]
metric_dt[, score_rank_norm := (score_rank - min(score_rank)) / (max(score_rank) - min(score_rank)), by = .(metric) ]

ggplot(metric_dt, aes(x = score, y = score_rank_norm)) +
    geom_point() +
    facet_wrap("metric", scales = "free_x")

id_oi = dt_motif[, .(score_range = quantile(top_motif_prop, .75)),
                 by = .(id)][order(score_range, decreasing = TRUE)]$id[1:3]
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
        flex_frag_len > 50 & flex_frag_len <= 150 & diff_corr > .1,
    combined_passing := TRUE ]


# qdt[signalValue > 10 & qValue > 10 & flex_frag_corr > .65 & read_corr < .6 &
# flex_frag_len >= 100 & flex_frag_len <= 300, passHigh := TRUE ]

npass = sum(qdt$combined_passing)
npass
nrow(qdt)

# todo_groups = colnames(qdt)[grepl("group", colnames(qdt))]
todo_groups_direction = rep(1, length(todo_groups))
names(todo_groups_direction) = todo_groups
todo_groups_direction["read_corr_group"] = -1

for(group in todo_groups){
    qdt[[sub("_group", "_passing", group)]] =
        factor(qdt[, rank(-1*todo_groups_direction[group]*get(sub("_group", "", group)), ties.method = "first") <= npass])
}


todo_passing = colnames(qdt)[grepl("passing", colnames(qdt))]

dt_motif_pass = calcMotifEnrichment(corr_res, base_gr, qdt, todo_passing,
                                    pwm, seq, bam_md5 = bam_md5,
                                    cache_path = "~/.cache_motif/",
                                    cach_version = cach_version)

# todo_dt = rbindlist(lapply(todo_passing, function(grp){
#     data.table(sel = c(FALSE, TRUE), grp = grp)
# }))
#
# pass_motif = parallel::mclapply(seq_len(nrow(todo_dt)), mc.preschedule = FALSE, function(i){
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
#     }, force_overwrite = FALSE)
#     myMotifRes
# })
#
# names(pass_motif) = paste(todo_dt$sel, todo_dt$grp)
#
# dt_motif_pass = rbindlist(lapply(pass_motif, function(y){
#     dt = as.data.table(as.data.frame(y))
#     colnames(dt) = gsub("\\.", "_", colnames(dt))
#     dt
# }), use.names = TRUE, idcol = "name")
#
# dt_motif_pass[, c("group", "metric") := tstrsplit(name, " ")]
#
# dt_motif_pass$group = factor(dt_motif_pass$group, levels = c(TRUE, FALSE))
# dt_motif_pass = dt_motif_pass[order(group)][order(metric)][order(id)]


id_oi_pass = dt_motif_pass[, .(score_range = quantile(top_motif_prop, .75)), by = .(id)][order(score_range, decreasing = TRUE)]$id[1]
ggplot(dt_motif_pass[id %in% id_oi_pass],
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

ggplot(dt_motif_pass[id %in% id_oi_pass],
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

ggplot(dt_motif_pass[id %in% id_oi_pass],
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

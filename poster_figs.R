library(seqsetvis)
library(data.table)
library(magrittr)
library(GenomicRanges)

setwd("~/R/peakrefine/")
setwd(file.path("/slipstream/galaxy/uploads/working/qc_framework"))
peaks = c(
    "output_JR_bookmarking_blocked_RUNX1_U13369masked/MCF10A-blocked_Runx1-4336BF_pooled/MCF10A-blocked_Runx1-4336BF_pooled_peaks.narrowPeak",
    "output_JR_bookmarking_full_RUNX1_U3369_masked/MCF10A-dmso_Runx1_pooled/MCF10A-dmso_Runx1_pooled_peaks.narrowPeak",
    "output_JR_bookmarking_full_RUNX1_U3369_masked/MCF10A-released_Runx1_pooled/MCF10A-released_Runx1_pooled_peaks.narrowPeak"
)
peaks = normalizePath(peaks)
bams = sub("_peaks.narrowPeak$", ".bam", peaks)
inputs = bams %>% gsub("Runx1-4336BF", "input", .) %>% gsub("Runx1", "input", .)

bam_files = c(bams[1], inputs[1])
names(bam_files) = c("blocked_Runx1", "blocked_input")
bam_files = bam_files[1]
setwd("~/R/peakrefine/")

load("results2/motif_res_MCF10A-blocked_Runx1-4336BF_pooled_data.save")

# qdt = qdt[diff_corr > .1 & stable_frag_corr > .9]

qdt = qdt[pValue < 10]
qdt = qdt[stable_frag_corr > .9]
# qgr = GRanges(qdt[order(pValue)])
qgr = GRanges(qdt[stable_frag_corr > read_corr][order(stable_frag_corr)])
qgr = resize(qgr, 1200, fix = "center")
names(qgr) = qgr$name


# pgr = qgr[c(1:6, length(qgr)-5:0)]
# pgr = qgr[c(length(qgr)-11:0)]
pgr = sample(qgr, 6)#[c(length(qgr)-11:0)]


f1_dt = ssvFetchBam(bam_files, pgr, target_strand = "both", return_data.table = TRUE, win_size = 200, fragLens = NA, win_method = "summary")
# f1_dt[, grep = ""]
# f1_dt = applySpline(f1_dt, n = 10, by_ = c("sample", "id", "strand"))
qdt[, id := name]


hiq_ids = c(23, 12, 1128, 5, 17, 102, 120, 1122, 162, 101)
toss_ids = c(1371, 975, 1148)
rescue_ids = c(643, 1179, 1039)
102 #corr decreased with shift
120 #high quality
1122
162
101

1371 #bad by corr
975
1148

643# rescue
1179
1039

f1_dt$id = factor(f1_dt$id, levels = as.character(pgr$name))

max_dt = f1_dt[, .(max_y = max(y)), by = .(id)]
ann_dt = qdt[name %in% f1_dt$id]
ann_dt = merge(ann_dt, max_dt)



ggplot(f1_dt) +
    geom_path(aes(x = x, y = y, color = strand, group = paste(sample, strand)), size = 1.5) +
    geom_label(data = ann_dt, aes(x = -.3, y = max_y*.8, label = round(stable_frag_corr, 3))) +
    geom_label(data = ann_dt, aes(x = 0, y = max_y*.8, label = round(pValue, 3))) +
    geom_label(data = ann_dt, aes(x = .3, y = max_y*.8, label = round(read_corr, 3))) +
    # geom_label(data = qdt[name %in% f1_dt$id], aes(x = -500, y = 20, label = round(stable_frag_corr, 3))) +
    # geom_label(data = qdt[name %in% f1_dt$id], aes(x = 500, y = 20, label = round(pValue, 3))) +
    scale_color_manual(values = c("+" = "darkred", "-" = "darkblue")) +
    facet_wrap("id", scales = "free_y", ncol = 2) + theme_classic()

ggplot(f1_dt[id == "peak_120" & sample == "blocked_Runx1"]) +
    geom_path(aes(x = x, y = y, color = strand, group = paste(sample, strand)), size = 2) +
    # geom_label(data = qdt[name %in% f1_dt$id], aes(x = 0, y = 30, label = stable_frag_corr)) +
    scale_color_manual(values = c("+" = "darkred", "-" = "darkblue")) +
    facet_grid("id~sample") + theme_classic()


ggplot(f1_dt[id == "peak_102" & sample == "blocked_Runx1"]) +
    geom_path(aes(x = x, y = y, color = strand, group = paste(sample, strand)), size = 2) +
    # geom_label(data = qdt[name %in% f1_dt$id], aes(x = 0, y = 30, label = stable_frag_corr)) +
    scale_color_manual(values = c("+" = "darkred", "-" = "darkblue")) +
    facet_grid("id~sample") + theme_classic()

load("results2/motif_res_MCF10A-blocked_Runx1-4336BF_pooled_data.save")
qdt

myplot = function(qid, bam_file = bams[1], show_stats = TRUE){#, corr_res = corr_res, qdt = qdt){
    qid = paste0("peak_", qid)
    rl = corr_res$read_length
    fl = corr_res$flex_fragment_correlation[treatment == "pulldown" & id == qid]$shift


    cr = corr_res$full_correlation_results[id == qid]
    cr = cr[order(shift)][order(treatment)]

    yrng = range(cr[treatment == "pulldown"]$correlation)
    purp = "#ad33b3"
    blu = "#009eeb"

    ann_dt = qdt[name %in% qid]
    fc = round(ann_dt[, stable_frag_corr], 3)
    rc = round(ann_dt[, read_corr], 3)
    pv = round(ann_dt[, pValue], 3)
    sv = round(ann_dt[, signalValue], 3)

    p_corr = ggplot(cr[treatment == "pulldown"], aes(x = shift, y = correlation)) +
        annotate("line", x = rep(rl, 2), y = yrng, color = purp, size = 1.5, lty = 2)
    if(fl > 100){
        p_corr = p_corr + annotate("line", x = rep(fl, 2), y = yrng, color = blu, size = 1.5, lty = 2)
    }
    p_corr = p_corr + geom_path(size = 1.5) + labs(y = "SCC", x= "shift") + theme_classic() +
        annotate("text", x = rl, y = mean(yrng), color = purp, label = round(rl, 3), hjust  = 1.1, size = 6)
    if(fl > 100){
        p_corr = p_corr + annotate("text", x = fl, y = mean(yrng), color = blu, label = round(fl, 3), hjust  = -.1, size = 6)
    }
    p_corr = p_corr +
        theme(text = element_text(size = 14),
              plot.subtitle = element_text(size = 8),
              axis.line = element_line(color = 'darkgray', size = 1),
              axis.ticks = element_line(color = 'darkgray', size = 1),
              axis.text = element_text(color = "darkgray"))

    qgr = GRanges(qdt[name == qid])
    names(qgr) = qgr$name
    bam_dt = ssvFetchBam(bam_file, qgr, target_strand = "both", return_data.table = TRUE, win_size = 40, fragLens = NA, win_method = "summary")
    bam_dt[, x := (start + end) / 2]
    bam_dt[, x := x - mean(x)]
    bam_dt = applySpline(bam_dt, 10, by_ = c("strand"))

    max_dt = bam_dt[, .(max_y = max(y)), by = .(id)]

    ann_dt[, id := name]
    ann_dt = merge(ann_dt, max_dt)

    p_pile = ggplot(bam_dt) +
        geom_path(aes(x = x, y = y, color = strand, group = paste(sample, strand)), size = 1.5) +
        # labs(title = qid) +
        # labs(subtitle = paste0("fragment r2 = ", fc, "\nread r2 = ", rc, "\n-log10 pval = ", pv, "\nsigalValue = ", sv), x = "bp", y = "read pileup") +
        scale_color_manual(values = c("+" = "#ff0000", "-" = "#0072eb")) +
        theme_classic() + guides(color = "none") +
        theme(text = element_text(size = 14),
              plot.subtitle = element_text(size = 8),
              axis.line = element_line(color = 'darkgray', size = 1),
              axis.ticks = element_line(color = 'darkgray', size = 1),
              axis.text = element_text(color = "darkgray")) +
        labs(x = "bp", y = "read pileup")
    # facet_wrap("id", scales = "free_y", ncol = 2) +

    # p_pile
    if(show_stats){
        p_pile = p_pile +
            labs(title = qid) +
            labs(subtitle = paste0("fragment r2 = ", fc,
                                   "\nread r2 = ", rc,
                                   "\n-log10 pval = ", pv,
                                   "\nsigalValue = ", sv))
        p_corr = p_corr +
            labs(title = qid) +
            labs(subtitle = paste0("fragment r2 = ", fc,
                                   "\nread r2 = ", rc,
                                   "\n-log10 pval = ", pv,
                                   "\nsigalValue = ", sv))
    }

    cowplot::plot_grid(p_pile, p_corr)
}
myplot(612)

# pdf("fig2_hiq.pdf", width = 4.5, height  = 2)
# for(id in c(toss_ids, rescue_ids, hiq_ids)){
#     p = myplot(id)
#     print(p)
# }
# dev.off()

pdf("fig2_hiq_stat.pdf", width = 4.5, height  = 7)
plist = lapply(c(23, 5, 643, 1179), function(x)myplot(x, show_stats = T))
cowplot::plot_grid(plotlist = plist, ncol = 1)
dev.off()

pdf("fig2_loq_stat.pdf", width = 4.5, height  = 7)
plist = lapply(c(1128, 1148, 17, 102), function(x)myplot(x, show_stats = T))
cowplot::plot_grid(plotlist = plist, ncol = 1)
dev.off()

ggplot(corr_res$full_correlation_results[id == "peak_102"], aes(x = shift, y = correlation, color = treatment)) + geom_path()
qdt[id == "peak_102"]

cr = corr_res$full_correlation_results[, .(y = mean(correlation, na.rm = TRUE)), by = .(treatment, shift)]
ggplot(cr[order(shift)], aes(x = shift, y = y, color = treatment)) + geom_path()

print(load("results3/motif_res_MCF10A-blocked_Runx1-4336BF_pooled_data.save"))
hist(qdt$read_corr, breaks = 100)
hist(qdt$read_corr, breaks = 100)
qdt$combined_passing = NULL
qdt$combined_passing = TRUE
qdt[read_corr > .8, combined_passing := FALSE]
table(qdt$combined_passing)
qdt[stable_frag_corr - read_corr < .1, combined_passing := FALSE]
table(qdt$combined_passing)
hist(qdt[combined_passing == TRUE]$flex_frag_len)
hist(qdt[combined_passing == FALSE]$flex_frag_len)
qdt[flex_frag_len < 100, combined_passing := FALSE]
qdt[flex_frag_len > 300, combined_passing := FALSE]
table(qdt$combined_passing)
qdf = qdt[pValue < 40, .(pValue, signalValue, stable_frag_corr, read_corr, flex_frag_len)]

GGally::ggpairs(qdf, mapping = aes(alpha = .002, shape = "21"))
ggsave("metric_matrix.pdf", width = 9, height = 9)

dt_motif_pass[id == "Hsapiens-jolma2013-RUNX3-2"][group == "FALSE"]

ggplot(qdt, aes(x = combined_passing, y = signalValue)) + geom_boxplot() + coord_cartesian(ylim = c(0,15))

# rtracklayer::import.bed("/slipstream/home/joeboyd/R/jrb_csaw/hg38.blacklist.bed")
black_gr = rtracklayer::import.bed("/slipstream/galaxy/uploads/working/IDEAS/black_lists/hg38.sort_blacklist.bed")
black_kept = subsetByOverlaps(black_gr, GRanges(qdt[combined_passing == TRUE]))
subsetByOverlaps(black_gr, GRanges(qdt[combined_passing == FALSE]))



dropSeqlevels(black_kept, setdiff(seqlevels(black_kept), as.character(seqnames(black_kept))))
black_kept = keepSeqlevels(black_kept, as.character(seqnames(black_kept)))
seqlevels(black_kept)
subsetByOverlaps(GRanges(qdt[combined_passing == TRUE]), black_kept)
seqsetvis::ssvSignalLineplot(seqsetvis::ssvFetchBam(bam_files, subsetByOverlaps(GRanges(qdt[combined_passing == TRUE]), black_kept), target_strand = "both", win_size = 10, fragLens = 150), color_ = "strand", group_ = "strand")


mets = c("combined_passing", "flex_frag_len_passing", "qValue_passing", "read_corr_passing", "signalValue_passing", "stable_frag_corr_passing")
drop_rate= t(sapply(mets, function(m){
    black_kept = subsetByOverlaps(black_gr, GRanges(qdt[get(m) == TRUE]))
    black_drop = subsetByOverlaps(black_gr, GRanges(qdt[get(m) == FALSE]))
    len = length(black_drop) / (length(black_drop) + length(black_kept))
    c(len)
}))

drop_dt = data.table(metric = colnames(drop_rate), drop_rate = drop_rate[1,])

motif_dt = dt_motif_pass[id == "Hsapiens-jolma2013-RUNX3-2"][order(raw_score)]
motif_dt$group = factor(motif_dt$group, levels = c("TRUE", "FALSE"))
motif_dt = motif_dt[metric %in% mets]
motif_dt$metric = factor(motif_dt$metric, levels = mets)
p_passing = ggplot(motif_dt, aes(x = metric, fill = group, y = top_motif_prop)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x = "", title = "Runx as top motif for best 700 peaks by metric") +
    scale_y_continuous(breaks = c(0, .4, .8)) +
    scale_fill_manual(values = c("TRUE" = "darkgreen", "FALSE" = "darkgray")) +
    theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
p_passing
ggsave("fig_passing.pdf", p_passing, width = 5, height = 3)
drop_dt$metric = factor(drop_dt$metric, levels = mets)
p_drop = ggplot(drop_dt, aes(x = metric, y = drop_rate)) +
    labs(title = "drop rate of 41 peaks in blacklisted regions", y = "drop rate", x = "") +
    geom_bar(stat = "identity") + scale_y_continuous(breaks = 0:5/5, limits = c(0,1)) +
    theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
p_drop
ggsave("fig_drop.pdf", p_drop, width = 2.5, height = 3)
# library(GenomicRanges)
# library(BiocFileCache)
# base_gr = GRanges(qdt)
# library(PWMEnrich.Hsapiens.background)
# data(PWMLogn.hg19.MotifDb.Hsap)
# library(BSgenome.Hsapiens.UCSC.hg38)
# pwm = PWMLogn.hg19.MotifDb.Hsap
# seq = Hsapiens
# cach_version = "x"
# ncores = 16
# todo_groups = "combined_passing"
# qdt$combined_passing = factor(qdt$combined_passing)
# dt_motif = calcMotifEnrichment(corr_res, base_gr, qdt, todo_groups,
#                                pwm, seq, bam_md5 = "manual1",
#                                cache_path = "~/.cache_motif/",
#                                cach_version = cach_version, ncores = ncores, force_overwrite_pre = F, force_overwrite = F)
#
# dt_motif[order(top_motif_prop)][id == "Hsapiens-jolma2013-RUNX3-2"]

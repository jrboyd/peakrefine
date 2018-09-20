# source("https://bioconductor.org/biocLite.R")
# biocLite("PWMEnrich")
# biocLite("PWMEnrich.Hsapiens.background")
# biocLite("BSgenome.Hsapiens.UCSC.hg38")

library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)

library(seqsetvis)
library(data.table)
library(peakrefine)
data("PWMLogn.hg19.MotifDb.Hsap")
setwd("~/R/peakrefine/")
# source("functions_motif.R")

ncores = 32

registerCoresPWMEnrich(ncores)
useBigMemoryPWMEnrich(TRUE)


# setwd(rdir)
# qgr = easyLoad_narrowPeak("MCF10A_CTCF_pooled/MCF10A_CTCF_pooled_peaks.narrowPeak")[[1]]
# bam_file = file.path(rdir, "MCF10A_CTCF_pooled/MCF10A_CTCF_pooled.bam")
# bam_input = file.path(rdir, "MCF10A_input_pooled/MCF10A_input_pooled.bam")
# todo_fl = c(50, 65, 100, 150, 195, 200, 205, 215, 230, 260, 300)
# id_oi = unique(names(PWMLogn.hg19.MotifDb.Hsap$pwms))
# id_oi = id_oi[grepl("CTCF", id_oi)]
# gray_fl = 100
# red_fl = 205
# out_pref = "AF_MCF10A_CTCF"

grepPWM = function(regex, pwm){
    names(pwm$pwms)[grepl("CTCF", names(pwm$pwms))]
}

myPWM = PWMLogn.hg19.MotifDb.Hsap

rdir = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/"
setwd(rdir)
pr_10a_ctcf = PeakRefiner(peak_set = rtracklayer::import("MCF10A_CTCF_pooled/MCF10A_CTCF_pooled_peaks.narrowPeak", format = "narrowPeak"),
                          bam_treat_file = file.path(rdir, "MCF10A_CTCF_pooled/MCF10A_CTCF_pooled.bam"),
                          bam_input_file = file.path(rdir, "MCF10A_input_pooled/MCF10A_input_pooled.bam"),
                          pwm = PWMLogn.hg19.MotifDb.Hsap,
                          target_pwm_names = grepPWM("CTCF", myPWM),
                          fragment_lengths = c(50, 65, 100, 150, 195, 200, 205, 215, 230, 260, 300),
                          color_overrides = c("100" = "gray", "205" = "red"),
                          output_prefix = "AF_MCF10A_CTCF")
setwd("~/R/peakrefine/")
load("cache/MCF10A_CTCF_pooled/cache_motif_res_200nbases_51679seq.save")
ScoreMotif = function(pr){
score_motif(pr@bam_treat_file, pr@bam_input_file, qgr = pr@peak_set, fl = pr@fragment_lengths[1], pwm = pr@pwm, out_dir = file.path("cache", pr@output_prefix), motif_res = motif_res)
}

# rdir = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/"
# setwd(rdir)
# qgr = easyLoad_narrowPeak("MCF10CA1_CTCF_pooled/MCF10CA1_CTCF_pooled_peaks.narrowPeak")[[1]]
# bam_file = file.path(rdir, "MCF10CA1_CTCF_pooled/MCF10CA1_CTCF_pooled.bam")
# bam_input = file.path(rdir, "MCF10CA1_input_pooled/MCF10CA1_input_pooled.bam")
# todo_fl = c(50, 100, 150, 195, 200, 205, 215, 240, 260, 300)
# id_oi = unique(names(PWMLogn.hg19.MotifDb.Hsap$pwms))
# id_oi = id_oi[grepl("CTCF", id_oi)]
# gray_fl = 100
# red_fl = 205
# out_pref = "AF_MCF10CA1_CTCF"

# rdir = "/slipstream/galaxy/uploads/working/qc_framework/output_MCF7_ESR1_enhancers"
# setwd(rdir)
# qgr = easyLoad_narrowPeak("MCF7-E2_ESR1_pooled/MCF7-E2_ESR1_pooled_peaks.narrowPeak")[[1]]
# bam_file = file.path(rdir, "MCF7-E2_ESR1_pooled/MCF7-E2_ESR1_pooled.bam")
# bam_input = file.path(rdir, "MCF7-E2_input_pooled/MCF7-E2_input_pooled.bam")
# todo_fl = c(37,100, 110, 120, 200)
# id_oi = unique(names(PWMLogn.hg19.MotifDb.Hsap$pwms))
# id_oi = id_oi[grepl("ESR[1A]", id_oi)]
# gray_fl = 37
# red_fl = 100
# out_pref = "ESR1_carroll"

# rdir = "/slipstream/galaxy/uploads/working/qc_framework/output_JR_bookmarking_blocked_RUNX1_U13369masked/"
# setwd(rdir)
# qgr = easyLoad_narrowPeak("MCF10A-blocked_Runx1-4336BF_pooled/MCF10A-blocked_Runx1-4336BF_pooled_peaks.narrowPeak")[[1]]
# bam_file = file.path(rdir, "MCF10A-blocked_Runx1-4336BF_pooled/MCF10A-blocked_Runx1-4336BF_pooled.bam")
# bam_input = file.path(rdir, "MCF10A-blocked_input_pooled/MCF10A-blocked_input_pooled.bam")
# todo_fl = c(100, 150, 200, 280)
# id_oi = c("Hsapiens-jolma2013-RUNX2-3")
# gray_fl = 100
# red_fl = 150
# out_pref = "JR_Runx1-mitotic"

##AF RUNX1
# rdir = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_RUNX1_ChIP/"
# setwd(rdir)
# qgr = easyLoad_narrowPeak("AF-MCF10A_RUNX1_pooled/AF-MCF10A_RUNX1_pooled_peaks.narrowPeak")[[1]]
# bam_file = file.path(rdir, "AF-MCF10A_RUNX1_pooled/AF-MCF10A_RUNX1_pooled.bam")
# bam_input = file.path(rdir, "AF-MCF10A_input_pooled/AF-MCF10A_input_pooled.bam")
# todo_fl = c(50, 65, 100, 175, 200, 280)
# id_oi = c("Hsapiens-jolma2013-RUNX2-3")
# gray_fl = 65
# red_fl = 175
# out_pref = "AF_MCF10A-Runx1"

# rdir = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_RUNX1_ChIP/"
# setwd(rdir)
# qgr = easyLoad_narrowPeak("AF-MCF10CA1_RUNX1_pooled/AF-MCF10CA1_RUNX1_pooled_peaks.narrowPeak")[[1]]
# bam_file = file.path(rdir, "AF-MCF10CA1_RUNX1_pooled/AF-MCF10CA1_RUNX1_pooled.bam")
# bam_input = file.path(rdir, "AF-MCF10CA1_input_pooled/AF-MCF10CA1_input_pooled.bam")
# todo_fl = c(50, 65, 100, 150, 175, 200, 280)
# gray_fl = 65
# red_fl = 150
# id_oi = c("Hsapiens-jolma2013-RUNX2-3")
# out_pref = "AF_MCF10CA1a-Runx1"

# bam_file = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled.bam"
# bam_input = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_input_pooled/MDA231_input_pooled.bam"
# qgr = easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled_peaks.narrowPeak")[[1]]
# todo_fl = c(50, 150, 180, 200, 220, 280)
# gray_fl = 150
# red_fl = 200
# id_oi = c("Hsapiens-jolma2013-RUNX2-3")
# out_pref = "MK_MDA231-Runx2"


qgr = subset(qgr, seqnames(qgr) != "chrU13369.1")
qgr$id = paste0("region_", seq_along(qgr))
names(qgr) = qgr$id
bgr = qgr

sc = ssvStrandCorr(bam_file, qgr, frag_max = 350)
sc$read_length
sc$frag_length
sc$sample_plot + labs(title = basename(bam_file),
                      subtitle = paste("read length:", sc$read_length,
                                       "\nfragment length:", sc$frag_length))

nbases = 200
if(!exists("motif_res")){
    cfile = file.path("~/R/peakrefine/cache", sub(".bam", "", basename(bam_file)), paste0("cache_motif_res_", nbases, "nbases_", length(qgr), "seq.save"))
    if(file.exists(cfile)){
        load(cfile)
    }else{
        motif_res = pre_motif(qgr)
    }
}
# todo_fl = c(50, 100, 175, 200, 280)
all_res = list()
for(fl in todo_fl){
    dt_motif = score_motif(bam_file, ncores = 16,
                           bam_input,
                           qgr,
                           fl,
                           motif_res = motif_res,
                           include_fl_independent = fl == todo_fl[1])
    all_res[[paste0("fl_", fl)]] = dt_motif
}



dt = rbindlist(all_res, use.names = TRUE, idcol = "fragLen")
dt$fragLen = factor(dt$fragLen, levels = paste0("fl_", sort(todo_fl)))
dt[, numFragLen := as.numeric(sub("fl_", "", fragLen))]
dt$metric = sub("_group", "", dt$metric)
dt$metric = factor(dt$metric, levels = c("pValue", "signalValue",
                                         "frag", "diff", "count", "input_count",
                                         "input_frag", "input_read", "read"))

mycol = c("gray", "black", "red")[rep(2, length(todo_fl))]
mycol[todo_fl == gray_fl] = "gray"
mycol[todo_fl == red_fl] = "red"
names(mycol) = paste0("fl_", sort(todo_fl))


# id_oi = id_oi[grepl("CTCF", id_oi)]
pdf(file.path("results", pdf_name), width = 18, height = length(id_oi)*3.5)
ggplot(dt[id %in% id_oi], aes(x = group, y = top_motif_prop, color = fragLen, group =  paste(id, metric, fragLen))) +
    geom_path() +
    geom_point() +
    ggrepel::geom_text_repel(data = dt[as.numeric(group) == max(as.numeric(group)) & id %in% id_oi],
                             mapping = aes(label = numFragLen)) +
    facet_grid("id~metric", scales = "free_y") +
    theme(panel.grid.major.y = element_line()) +
    # scale_x_discrete(breaks = c(1, 8), labels = c("worst", "best")) +
    scale_color_manual(values = mycol) + theme_classic()


ggplot(dt[id %in% id_oi], aes(x = group, y = raw_score, color = fragLen, group =  paste(id, metric, fragLen))) +
    geom_path() +
    geom_point() +
    ggrepel::geom_text_repel(data = dt[as.numeric(group) == max(as.numeric(group)) & id %in% id_oi],
                             mapping = aes(label = numFragLen)) +
    facet_grid("id~metric", scales = "free_y") +
    theme(panel.grid.major.y = element_line()) +
    scale_color_manual(values = mycol) + theme_classic()
dev.off()


strand_file  = dir(file.path("~/R/peakrefine/cache", sub(".bam", "", basename(bam_file))), pattern = paste0("strandRes_", as.character(red_fl)), full.names = TRUE)
load(strand_file)
sr = as.data.table(strandRes)
png(file.path("results", sub("metrics.pdf", "matrix.png", pdf_name)), width = 8, height = 8, units = "in", res= 300)
plot(sr[, .(read_corr, frag_corr, diff_corr = frag_corr - read_corr, signalValue, pValue)],
     pch = 16, col = rgb(0,0,0,.02))
dev.off()


pdt = data.table(x = seq_len(nrow(sr)), y = sort(sr$frag_corr))
plot(pdt)
pdt[, x := (x - min(x)) / (max(x) - min(x))]
pdt[, y := (y - min(y)) / (max(y) - min(y))]
pdt[, y := seqsetvis:::movingAverage(y,n = 100)]
ls = loess(pdt$y ~ pdt$x)
pr.loess <- predict(ls)
lines(pr.loess~pdt$x, col = "red")
fit3 <- lm(pdt$y~poly(pdt$x,3,raw=TRUE))
pr.fit3 = predict(fit3)
lines(pr.fit3~pdt$x, col = "blue")


stagger = cbind(pdt[-nrow(pdt)], pdt[-1,])
colnames(stagger) = paste0(colnames(stagger), rep(1:2, each = 2))
stagger[, m := (y2 - y1) / (x2 - x1)]
stagger[, xm := (x1 + x2) / 2]
plot(stagger[, .(xm, m)], ylim = c(0,50))
plot(seqsetvis:::movingAverage(stagger$m, n = 100), ylim = c(0,10), pch = 16, col = rgb(0,0,0,.02))
# min_p = min(dt[p_value > 0]$p_value)
# dt[p_value == 0, p_value := min_p]
#
# id_oi = "ABCF2"
# ggplot(dt[(id %in% id_oi)], aes(x = group, y = -log10(p_value), color = fragLen, group =  paste(id, metric, fragLen))) +
#     geom_path() +
#     geom_point() +
#     ggrepel::geom_text_repel(data = dt[as.numeric(group) == max(as.numeric(group)) & id %in% id_oi],
#                              mapping = aes(label = numFragLen)) +
#     facet_grid("id~metric") +
#     theme(panel.grid.major.y = element_line()) +
#     scale_color_manual(values = mycol) + theme_classic()
#
# id_oi2 = dt[, .N, by = .(id)][N == 168]$id
#
# resClust = ssvSignalClustering(dt[metric == "frag_group" & id %in% id_oi2],
#                     row_ = "id",
#                     column_ = "group",
#                     facet_ = "fragLen",
#                     fill_ = "top_motif_prop", max_rows = Inf)
# resMap = seqsetvis::ssvSignalHeatmap(resClust,
#                                   row_ = "id",
#                                   column_ = "group",
#                                   facet_ = "fragLen",
#                                   fill_ = "top_motif_prop")
# resMap
# resClust[cluster_id == "1"][order(top_motif_prop, decreasing = TRUE)][, unique(id)]
#
# dt[, lg_p_value := -log10(p_value)]
# resClust = ssvSignalClustering(dt[metric == "frag_group" & id %in% id_oi2], nclust = 12,
#                                row_ = "id",
#                                column_ = "group",
#                                facet_ = "fragLen",
#                                fill_ = "lg_p_value", max_rows = Inf)
# resMap = seqsetvis::ssvSignalHeatmap(resClust,
#                                      row_ = "id",
#                                      column_ = "group",
#                                      facet_ = "fragLen",
#                                      fill_ = "lg_p_value")

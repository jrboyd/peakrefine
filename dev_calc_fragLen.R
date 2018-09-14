library(seqsetvis)
library(peakrefine)
library(GenomicRanges)
library(ggplot2)
library(data.table)
load("../KZ_Runx1_BRG1_ESR1_overlap/tmp.runx.save")

# qgr = CTCF_in_10a_overlaps_gr[1:5]


theme_set(theme_classic() + theme(strip.background = element_blank()))



# qgr = GenomicRanges::shift(qgr, 5000)
# bam_file = system.file("extdata/test.bam", package = "seqsetvis")

todo = "MCF10A_CTCF"
todo = "T47D_BRG1"
todo = "MDA231_RUNX2"
switch(todo,
       MCF10A_CTCF = {
           qgr = easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled_peaks_passIDR.05.narrowPeak")[[1]]
           bam_file = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled.bam"
           bam_input = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_input_pooled/MCF10A_input_pooled.bam"
           roi = paste0("region_", c(15131, 15345, 16463, 17975, 19150, 33105))
       },
       T47D_BRG1 = {
           qgr = easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_GSE112491_T47D_BRG1_rDNAmasked/A1A3_BRG1_R1/A1A3_BRG1_R1_peaks.narrowPeak")[[1]]
           bam_file = "/slipstream/galaxy/uploads/working/qc_framework/output_GSE112491_T47D_BRG1_rDNAmasked/A1A3_BRG1_R1/A1A3_BRG1_R1.bam"
           bam_input = "/slipstream/galaxy/uploads/working/qc_framework/output_GSE112491_T47D_BRG1_rDNAmasked/A1A3_input_R1/A1A3_input_R1.bam"
           roi = paste0("region_", c(3281, 3291, 3292, 3260, 5486, 3263))
       },
       MDA231_RUNX2 = {
           qgr = easyLoad_narrowPeak(runx2_files$peaks_pooled)[[1]]
           bam_file = runx2_files$bam_pooled[1]
           bam_input = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_input_pooled/MDA231_input_pooled.bam"
       }

)




qgr = centerFixedSizeGRanges(qgr, 800)
qgr$id = paste0("region_", seq_along(qgr))
names(qgr) = qgr$id

n = 6
revbest = FALSE
# sdt = ssvRecipes::myFetchStrandedBam(bam_file, qgr[order(qgr$pValue, decreasing = revbest)][1:n])
# p1 = ggplot(sdt[id %in% unique(sdt$id)[1:n]], aes(x = x, y = y, color = strand)) +
#     geom_path() +
#     facet_wrap("id", scales = "free_y") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + labs(title = "no ext")
# #
# ext = 200
# sdt = ssvRecipes::myFetchStrandedBam(bam_file, qgr[order(qgr$pValue, decreasing = revbest)][1:n], fragLens = ext)
# p2 = ggplot(sdt[id %in% unique(sdt$id)[1:n]], aes(x = x, y = y, color = strand)) +
#     geom_path() +
#     facet_wrap("id", scales = "free_y") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + labs(title = paste(ext, "ext"))
# cowplot::plot_grid(p1, p2)


badgr = GenomicRanges::shift(qgr, 5000)
names(badgr) = paste0(names(qgr), "bad")

biggr = GenomicRanges::resize(qgr, 5000, fix = "center")

bgr = qgr[order(qgr$pValue, decreasing = TRUE)][1:12]



cc = ssvCrossCorr(bam_file, bgr, shift_min = 0, shift_max = 300, step = 10, nbest = 12, revbest = revbest)


corr_res = ssvStrandCorr(bam_file, bgr, frag_min = 1, frag_max = 300)

p1 = ggplot(cc, aes(x = shiftLen, y = corr)) +
    facet_wrap("id") +
    geom_path() +
    labs(title = "cross corr")
p2 = corr_res$sample_plot +
    labs(title = "frag corr")


fl = corr_res$frag_length
fl = 180
rl = corr_res$read_length
cowplot::plot_grid(p1 +
                       annotate("line", x = rep(fl, 2), y = c(-.5, 1), color = "green") +
                       annotate("line", x = rep(rl, 2), y = c(-.5, 1), color = "red"),
                   p2 +
                       annotate("line", x = rep(fl, 2), y = c(-.5, 1), color = "green") +
                       annotate("line", x = rep(rl, 2), y = c(-.5, 1), color = "red")
)

bdt = ssvRecipes::myFetchStrandedBam(bam_file, qgr[roi], fragLens = NA, max_dupes = 1)
p1 = ggplot(bdt, aes(x = x, y = y, color = strand)) + geom_path() + labs(title = roi) +
    facet_wrap("id")

bdt = ssvRecipes::myFetchStrandedBam(bam_file, qgr[roi], fragLens = fl, max_dupes = 1)
p2 = ggplot(bdt, aes(x = x, y = y, color = strand)) + geom_path() + labs(title = roi) +
    facet_wrap("id")

cowplot::plot_grid(p1 + labs(title = paste0("reads (", rl, ")")),
                   p2 + labs(title = paste0("fragments (", fl, ")"))
)

strandRes = ssvStrandCorrFull(bam_file, qgr, fragLen = fl, ncores = 16, output_withGRanges = TRUE)
strandResInput = ssvStrandCorrFull(bam_input, qgr, fragLen = fl, ncores = 16, output_withGRanges = TRUE)

ggplot(as.data.table(strandRes), aes(x = frag_corr, y = pValue)) + geom_point(alpha = .2, shape = 16) + geom_density2d()
ggplot(as.data.table(strandRes), aes(x = frag_corr - read_corr, y = pValue)) + geom_point(alpha = .2, shape = 16) + geom_density2d()

ggplot(as.data.table(strandResInput), aes(x = frag_corr, y = pValue)) + geom_point(alpha = .2, shape = 16) + geom_density2d()
ggplot(as.data.table(strandResInput), aes(x = frag_corr - read_corr, y = pValue)) + geom_point(alpha = .2, shape = 16) + geom_density2d()

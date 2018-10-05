source("dev_pipeline_fun.R")
library(magrittr)

cach_version = "v6"

# bam_file = "~/ENCODE_EGR1/K562_EGR1_myers_rep1_ENCFF645IYJ.bam"
# qgr = rtracklayer::import("~/ENCODE_EGR1/K562_EGR1_myers_rep1_ENCFF178PTI.bed", format = "narrowPeak")
# pipeline(bam_file, qgr, cach_version)
#
#
# bam_file = "~/ENCODE_EGR1/K562_EGR1_myers_rep2_ENCFF021CMK.bam"
# qgr = rtracklayer::import("~/ENCODE_EGR1/K562_EGR1_myers_rep2_ENCFF913QDG.bed", format = "narrowPeak")
# pipeline(bam_file, qgr, cach_version)

#RUNX META
peaks = dir("/slipstream/home/jonathan/RUNX_META-ANALYSIS/MACS/", full.names = TRUE) %>%
    dir(full.names = TRUE, pattern = "Peak$")# %>% dir(full.names = TRUE, pattern = "Peak$")

bams = paste0(file.path("/slipstream/home/jonathan/RUNX_META-ANALYSIS/BAMS", basename(dirname(peaks))), ".bam")
file.exists(bams)
todo_df = data.frame(bam = bams,
                     peak = peaks,
                     sample = sapply(basename(bams) %>% strsplit(., "_"), function(x)paste(x[1:4], collapse = "_")),
                     stringsAsFactors = FALSE)

all_res = lapply(seq_len(nrow(todo_df)), function(i){
    bam_file = todo_df$bam[i]
    print(bam_file)
    frag_max = 250


    if(file.exists(paste0(bam_file, ".bai"))){
        qgr = rtracklayer::import(todo_df$peak[i], format = "narrowPeak")
        # pipeline(bam_file, qgr, cach_version)
        get_corr_res(bam_file, qgr, cach_version)
    }
})
names(all_res) = todo_df$sample
format(object.size(all_res), units = "GB")
k = sapply(all_res, function(x)x$read_length) == 35
todo_df[k,]
table(sapply(all_res, function(x)x$read_length))
all_corr = lapply(all_res, function(x)x$full_correlation_results)
all_corr = lapply(all_corr, function(x){
    uid = unique(x$id)
    x[id %in% sample(uid, min(500, length(uid)))]
})
all_corr = rbindlist(all_corr, use.names = TRUE, idcol = "sample")
# ggplot(all_corr[sample == unique(all_corr$sample)[12]],
#        aes(x = shift, y = correlation, group = shift)) +
#     geom_boxplot() +
#     facet_wrap("sample")

all_corr = all_corr[!is.nan(correlation)]
range(all_corr$correlation)


library(ggplot2)
theme_set(theme_classic())
pdf("tmp.pdf")
for(i in seq_len(nrow(todo_df))){
    samp = todo_df$sample[i]
    # samp = unique(all_corr$sample)[i]
    rl = all_res[[i]]$read_length
    fl = all_res[[i]]$fragment_length
    p = seqsetvis::ssvSignalBandedQuantiles(all_corr[sample == samp],
                                        x_ = "shift", y_ = "correlation", by_ = "sample",
                                        hsv_symmetric = TRUE, hsv_grayscale = TRUE, hsv_reverse = TRUE) +
        labs(title = samp) + guides(fill = "none") +
        annotate("line", x = rep(rl, 2), y = range(all_corr[sample == samp]$correlation), color = "blue", size = 2) +
        annotate("label", x = rl, y = mean(range(all_corr[sample == samp]$correlation)), color = 'blue', label = "read") +
        annotate("line", x = rep(fl, 2), y = range(all_corr[sample == samp]$correlation), color = "red", size = 2) +
        annotate("label", x = fl, y = mean(range(all_corr[sample == samp]$correlation)), color = 'red', label = "fragment")
    print(p)
}
dev.off()

# for(i in seq_len(nrow(todo_df))){
#     bam_file = todo_df$bam[i]
#     print(bam_file)
#     if(file.exists(paste0(bam_file, ".bai"))){
#         qgr = rtracklayer::import(todo_df$peak[i], format = "narrowPeak")
#         pipeline(bam_file, qgr, cach_version)
#     }
# }

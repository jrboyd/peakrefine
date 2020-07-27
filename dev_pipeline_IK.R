setwd("~/R/peakrefine/")
source("dev_pipeline_fun.R")
library(magrittr)

theme_set(theme_classic())
cach_version = "v1"
# data_source = "bookmarking_vs_blocked"
data_source = "SF_IK"
ncores = 32
# bam_file = "~/ENCODE_EGR1/K562_EGR1_myers_rep1_ENCFF645IYJ.bam"
# qgr = rtracklayer::import("~/ENCODE_EGR1/K562_EGR1_myers_rep1_ENCFF178PTI.bed", format = "narrowPeak")
# pipeline(bam_file, qgr, cach_version)
#
#
# bam_file = "~/ENCODE_EGR1/K562_EGR1_myers_rep2_ENCFF021CMK.bam"
# qgr = rtracklayer::import("~/ENCODE_EGR1/K562_EGR1_myers_rep2_ENCFF913QDG.bed", format = "narrowPeak")
# pipeline(bam_file, qgr, cach_version)

#RUNX META

suppressWarnings({
    remove("gen")
    remove("peaks")
    remove("bams")
    remove("inputs")
})
gen = "mm10"
skip_motif = FALSE
k = TRUE
frag_max = 400
if(data_source == "SF_IK"){
    peaks = ik_peak_files = dir("~/R/SF_AutoImmune_ssv/data/chromatin/IKAROS/peaks/macs2/", full.names = TRUE)
    peaks = peaks[grepl("IK[0-9]_", peaks)]

    bams = dir("~/R/SF_AutoImmune_ssv/data/chromatin/IKAROS/bam/", pattern = ".bam$", full.names = TRUE)
    inputs = bams[grepl("IgG\\.", bams)]

    bams = bams[grepl("IK[0-9]\\.", bams)]

    inputs = rep(inputs, length(bams))

    peaks = normalizePath(peaks)
    bams = normalizePath(bams)
    inputs = normalizePath(inputs)

    odir = paste0("result_", data_source, "_v1")

}
dir.create(odir, showWarnings = FALSE)
setwd("~/R/SF_AutoImmune_ssv/")
stopifnot(all(file.exists(peaks)))
stopifnot(all(file.exists(bams)))
force_overwrite_motif = TRUE
if(exists("inputs")){
    stopifnot(all(file.exists(inputs)))
    todo_df = data.frame(bam = bams,
                         peak = peaks,
                         inputs = inputs,
                         sample = sapply(basename(bams) %>% strsplit(., "_"), function(x)paste(x[1:4], collapse = "_")),
                         stringsAsFactors = FALSE)
    k = which(grepl("Hollen", todo_df$sample) | grepl("Chig", todo_df$sample))
    length(k)
    all_res = lapply(seq_len(nrow(todo_df)), function(i){
        bam_file = todo_df$bam[i]
        input_file = todo_df$inputs[i]
        message(i, " ", bam_file)


        if(file.exists(paste0(bam_file, ".bai"))){
            if(exists("over_gr")){
                qgr = over_gr
            }else{
                qgr = rtracklayer::import(todo_df$peak[i], format = "narrowPeak")
            }

            if(i %in% k){
                pipeline(bam_file, qgr, cach_version, frag_max = frag_max,
                         inputs_file = input_file, gen = gen, ncores = ncores,
                         output_dir = odir,
                         to_score = c("signalValue", "qValue",
                                      "stable_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len",
                                      "stable_frag_corr_input", "flex_frag_corr_input", "read_corr_input", "flex_frag_len_input"),
                         skip_motif = skip_motif, force_overwrite_motif = force_overwrite_motif)
                # get_corr_res(bam_file, qgr, cach_version, frag_max = 500)
            }else{
                pipeline(bam_file, qgr, cach_version,
                         inputs_file = input_file,
                         output_dir = odir,
                         gen = gen, ncores = ncores, frag_max = frag_max,
                         to_score = c("signalValue", "qValue",
                                      "stable_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len",
                                      "stable_frag_corr_input", "flex_frag_corr_input", "read_corr_input", "flex_frag_len_input"),
                         skip_motif = skip_motif, force_overwrite_motif = force_overwrite_motif)
                # get_corr_res(bam_file, qgr, cach_version)
            }

            # get_corr_res(bam_file, qgr, cach_version, force_overwrite = TRUE)
        }
    })
}else{
    todo_df = data.frame(bam = bams,
                         peak = peaks,
                         sample = sapply(basename(bams) %>% strsplit(., "_"), function(x)paste(x[1:4], collapse = "_")),
                         stringsAsFactors = FALSE)


    k = which(grepl("Hollen", todo_df$sample) | grepl("Chig", todo_df$sample))
    length(k)
    all_res = lapply(seq_len(nrow(todo_df)), function(i){
        bam_file = todo_df$bam[i]
        message(i, " ", bam_file)


        if(file.exists(paste0(bam_file, ".bai"))){
            if(exists("over_gr")){
                qgr = over_gr
            }else{
                qgr = rtracklayer::import(todo_df$peak[i], format = "narrowPeak")
            }

            if(i %in% k){
                pipeline(bam_file, qgr, cach_version,
                         output_dir = odir,
                         frag_max = frag_max, gen = gen,
                         ncores = ncores,
                         skip_motif = skip_motif, force_overwrite_motif = force_overwrite_motif)
                # get_corr_res(bam_file, qgr, cach_version, frag_max = 500)
            }else{
                pipeline(bam_file, qgr, cach_version,
                         output_dir = odir,
                         gen = gen, ncores = ncores,
                         frag_max = frag_max,
                         skip_motif = skip_motif, force_overwrite_motif = force_overwrite_motif)
                # get_corr_res(bam_file, qgr, cach_version)
            }

            # get_corr_res(bam_file, qgr, cach_version, force_overwrite = TRUE)
        }
    })
}

names(all_res) = todo_df$sample
format(object.size(all_res), units = "GB")
# k = sapply(all_res, function(x)x$read_length) == 35
# todo_df[k,]
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

pdf(paste0(odir, "/crosscorr_", data_source, ".pdf"))
for(i in seq_len(nrow(todo_df))){
    samp = todo_df$sample[i]
    # samp = unique(all_corr$sample)[i]
    rl = all_res[[i]]$read_length
    fl = all_res[[i]]$fragment_length
    ylim = c(-.5, 1)
    p = seqsetvis::ssvSignalBandedQuantiles(all_corr[sample == samp],
                                            x_ = "shift", y_ = "correlation", by_ = "sample",
                                            hsv_symmetric = TRUE, hsv_grayscale = TRUE, hsv_reverse = TRUE) +
        labs(title = samp) + guides(fill = "none") +
        coord_cartesian(ylim = ylim) +
        annotate("line", x = rep(rl, 2), y =ylim, color = "blue", size = 2) +
        annotate("label", x = rl, y = mean(ylim), color = 'blue', label = "read") +
        annotate("line", x = rep(fl, 2), y = ylim, color = "red", size = 2) +
        annotate("label", x = fl, y = mean(ylim), color = 'red', label = "fragment")
    print(p)

    p = ggplot(all_corr[sample == samp][id %in% sample(unique(id), 500)], aes(x = shift, y = correlation, group = id)) + geom_path(alpha = .05, size = 3, color = "black") +
        labs(title = samp) + guides(fill = "none")
    print(p)

    p = seqsetvis::ssvSignalHeatmap(all_corr[sample == samp], max_cols = Inf, nclust = 2, max_rows = Inf,
                                    column_ = "shift", row_ = "id",
                                    facet_ = "sample", fill_ = "correlation")
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


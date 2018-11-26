setwd("~/R/peakrefine/")
source("dev_pipeline_fun.R")
library(magrittr)

theme_set(theme_classic())
cach_version = "v8"
data_source = "bookmarking"

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
gen = "hg38"
k = TRUE
if(data_source == "Jonathan"){
    peaks = dir("/slipstream/home/jonathan/RUNX_META-ANALYSIS/MACS/", full.names = TRUE) %>%
        dir(full.names = TRUE, pattern = "Peak$")# %>% dir(full.names = TRUE, pattern = "Peak$")

    bams = paste0(file.path("/slipstream/home/jonathan/RUNX_META-ANALYSIS/BAMS", basename(dirname(peaks))), ".bam")
}else if(data_source == "bookmarking"){
    cach_version = "v8"
    setwd(file.path("/slipstream/galaxy/uploads/working/qc_framework"))
    peaks = c(
        "output_JR_bookmarking_blocked_RUNX1_U13369masked/MCF10A-blocked_Runx1-4336BF_pooled/MCF10A-blocked_Runx1-4336BF_pooled_peaks.narrowPeak",
        "output_JR_bookmarking_full_RUNX1_U3369_masked/MCF10A-dmso_Runx1_pooled/MCF10A-dmso_Runx1_pooled_peaks.narrowPeak",
        "output_JR_bookmarking_full_RUNX1_U3369_masked/MCF10A-released_Runx1_pooled/MCF10A-released_Runx1_pooled_peaks.narrowPeak"
    )
    peaks = normalizePath(peaks)
    bams = sub("_peaks.narrowPeak$", ".bam", peaks)
    inputs = bams %>% gsub("Runx1-4336BF", "input", .) %>% gsub("Runx1", "input", .)


}else if(data_source == "k27ac"){
    setwd(file.path("/slipstream/home/joeboyd/jonathan_MSC_k27ac_timecourse/"))
    peaks = c("narrowcall/h3k27ac_bmsc_d0_pooled_peaks.narrowPeak",
              "narrowcall/h3k27ac_bmsc_d21_pooled_peaks.narrowPeak")
    peaks = normalizePath(peaks)
    bams = c("MS_MSC_HISTONE_CODE_DATA_HWJG/H3K27AC/h3k27ac_bmsc_d0_COMB_GAII_sequence_sort_bowtie2_mm10.bam",
             "MS_MSC_HISTONE_CODE_DATA_HWJG/H3K27AC/h3k27ac_bmsc_d21_COMB_GAII_sequence_sort_bowtie2_mm10.bam")
    bams = normalizePath(bams)
    inputs = c("MS_MSC_HISTONE_CODE_DATA_HWJG/INPUT/input_bmsc_SR100_d00_COMB_SE100_sequence_sort_bowtie2_mm10.bam",
               "MS_MSC_HISTONE_CODE_DATA_HWJG/INPUT/input_bmsc_SR100_d21_COMB_SE100_sequence_sort_bowtie2_mm10.bam")
    inputs = normalizePath(inputs)
    gen = "mm10"
    # library(GenomicRanges)
    # library(data.table)
    # load("~/jonathan_MSC_k27ac_timecourse/DB_csaw_results.JG_k27ac_v5.save")
    # sel_gr = GRanges(comb_res$h3k27ac$`from d0 to d21`)
    # library(seqsetvis)
    # over_gr = ssvOverlapIntervalSets(easyLoad_narrowPeak(peaks), ext = 500)
    # over_gr = subsetByOverlaps(over_gr, sel_gr)
    # hist(width(over_gr), xlim = c(0,20000), breaks = 50)
}
setwd("~/R/peakrefine/")
stopifnot(all(file.exists(bams)))
if(exists("inputs")){
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
        frag_max = 250


        if(file.exists(paste0(bam_file, ".bai"))){
            if(exists("over_gr")){
                qgr = over_gr
            }else{
                qgr = rtracklayer::import(todo_df$peak[i], format = "narrowPeak")
            }

            if(i %in% k){
                pipeline(bam_file, qgr, cach_version, frag_max = 500, inputs_file = input_file, gen = gen,
                         to_score = c("signalValue", "qValue",
                                      "stable_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len",
                                      "stable_frag_corr_input", "flex_frag_corr_input", "read_corr_input", "flex_frag_len_input"))
                # get_corr_res(bam_file, qgr, cach_version, frag_max = 500)
            }else{
                pipeline(bam_file, qgr, cach_version, inputs_file = input_file, gen = gen,
                         to_score = c("signalValue", "qValue",
                                      "stable_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len",
                                      "stable_frag_corr_input", "flex_frag_corr_input", "read_corr_input", "flex_frag_len_input"))
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
        frag_max = 250


        if(file.exists(paste0(bam_file, ".bai"))){
            if(exists("over_gr")){
                qgr = over_gr
            }else{
                qgr = rtracklayer::import(todo_df$peak[i], format = "narrowPeak")
            }

            if(i %in% k){
                pipeline(bam_file, qgr, cach_version, frag_max = 500, gen = gen)
                # get_corr_res(bam_file, qgr, cach_version, frag_max = 500)
            }else{
                pipeline(bam_file, qgr, cach_version, gen = gen)
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

pdf(paste0("results2/crosscorr_", data_source, ".pdf"))
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


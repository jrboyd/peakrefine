setwd("~/R/peakrefine/")
source("dev_pipeline_fun.R")
library(magrittr)

theme_set(theme_classic())
cach_version = "v8"
data_source = "sim"

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
skip_motif = FALSE
k = TRUE
if(data_source == "sim"){
    gen = "simGenome10M_v4"
    skip_motif = TRUE
    peaks = dir("simulation/genomes/simGenome10M_v5/output/", full.names = TRUE, pattern = "pooled$") %>%
        dir(full.names = TRUE, pattern = "mediumLoose.+Peak$")

    bams = dir("simulation/genomes/simGenome10M_v5/output/", full.names = TRUE, pattern = "pooled$") %>%
        dir(full.names = TRUE, pattern = "_1fe_pooled.bam$")
    inputs = gsub("_[0-9]+fe_", "_input_", bams)

    peaks = normalizePath(peaks)
    bams = normalizePath(bams)
    inputs = normalizePath(inputs)

}
setwd("~/R/peakrefine/")
stopifnot(all(file.exists(peaks)))
stopifnot(all(file.exists(bams)))
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
                                      "stable_frag_corr_input", "flex_frag_corr_input", "read_corr_input", "flex_frag_len_input"),
                         skip_motif = skip_motif)
                # get_corr_res(bam_file, qgr, cach_version, frag_max = 500)
            }else{
                pipeline(bam_file, qgr, cach_version, inputs_file = input_file, gen = gen,
                         to_score = c("signalValue", "qValue",
                                      "stable_frag_corr", "flex_frag_corr", "read_corr", "flex_frag_len",
                                      "stable_frag_corr_input", "flex_frag_corr_input", "read_corr_input", "flex_frag_len_input"),
                         skip_motif = skip_motif, force_overwrite_corr = TRUE)
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
                pipeline(bam_file, qgr, cach_version, frag_max = 500, gen = gen,
                         skip_motif = skip_motif)
                # get_corr_res(bam_file, qgr, cach_version, frag_max = 500)
            }else{
                pipeline(bam_file, qgr, cach_version, gen = gen,
                         skip_motif = skip_motif)
                # get_corr_res(bam_file, qgr, cach_version)
            }

            # get_corr_res(bam_file, qgr, cach_version, force_overwrite = TRUE)
        }
    })
}

names(all_res) = todo_df$sample
format(object.size(all_res), units = "GB")

all_qdt = lapply(all_res, function(x)x$qdt)
all_corr = lapply(all_res, function(x)x$corr_res)



for(i in seq_along(all_res)){
    ver = sub("_.+", "", names(all_res)[i]) %>% sub("V", "", .)
    ground_file = paste0("simulation/genomes/simGenome10M_v4/peaks/simPeaks_v", ver, ".bed")
    ground_gr = rtracklayer::import.bed(ground_file)
    corr_gr = GRanges(all_res[[i]]$qdt)
    olaps = findOverlaps(corr_gr, ground_gr)
    corr_gr$is_true = FALSE
    corr_gr$is_true[queryHits(olaps)] = TRUE
    corr_gr$true_weight = 0
    corr_gr$true_weight[queryHits(olaps)] = ground_gr$score[subjectHits(olaps)]

    ground_gr$is_called = FALSE
    ground_gr$is_called[subjectHits(olaps)] = TRUE

    df = as.data.frame(ground_gr)
    ggplot(df, aes(x = is_called, y = score)) + geom_point()
}


# all_corr = rbindlist(all_res, use.names = TRUE, idcol = "sample")
# ggplot(all_corr[sample == unique(all_corr$sample)[12]],
#        aes(x = shift, y = correlation, group = shift)) +
#     geom_boxplot() +
#     facet_wrap("sample")

# all_corr = all_corr[!is.nan(correlation)]
# range(all_corr$correlation)


full_corr = lapply(all_corr, function(x)x$full_correlation_results)
names(full_corr) = todo_df$sample
full_corr = rbindlist(full_corr, use.names = TRUE, idcol = "sample")

library(ggplot2)
dir.create("results_sim", showWarnings = FALSE)
pdf(paste0("results_sim/crosscorrSim_", data_source, ".pdf"))
full_corr[is.nan(correlation), correlation := 0]
for(i in seq_len(nrow(todo_df))){
    samp = todo_df$sample[i]
    # samp = unique(all_corr$sample)[i]
    rl = all_corr[[i]]$read_length
    fl = all_corr[[i]]$fragment_length
    ylim = c(-.5, 1)
    # full_corr = all_corr[[i]]$full_correlation_results
    p = seqsetvis::ssvSignalBandedQuantiles(full_corr[sample == samp],
                                            x_ = "shift", y_ = "correlation", by_ = "sample",
                                            hsv_symmetric = TRUE, hsv_grayscale = TRUE, hsv_reverse = TRUE) +
        labs(title = samp) + guides(fill = "none") +
        coord_cartesian(ylim = ylim) +
        annotate("line", x = rep(rl, 2), y =ylim, color = "blue", size = 2) +
        annotate("label", x = rl, y = mean(ylim), color = 'blue', label = "read") +
        annotate("line", x = rep(fl, 2), y = ylim, color = "red", size = 2) +
        annotate("label", x = fl, y = mean(ylim), color = 'red', label = "fragment")
    print(p)

    p = ggplot(full_corr[sample == samp][id %in% sample(unique(id), min(500, length(unique(id))))], aes(x = shift, y = correlation, group = id)) + geom_path(alpha = .05, size = 3, color = "black") +
        labs(title = samp) + guides(fill = "none")
    print(p)

    p = seqsetvis::ssvSignalHeatmap(full_corr[sample == samp], max_cols = Inf, nclust = 2, max_rows = Inf,
                                    column_ = "shift", row_ = "id",
                                    facet_ = "treatment", fill_ = "correlation")
    print(p)
}
dev.off()

all_qdt = lapply(all_qdt, function(x)x[, !grepl("group", colnames(x)), with = FALSE])
names(all_qdt) = todo_df$sample
all_qdt = rbindlist(all_qdt, use.names = TRUE, idcol = "sample")
all_qdt[, ver := sub("_.+", "", sample)]
all_qdt[, name := paste(ver, name)]

all_qdt[, seqnames := ver]

true_gr =  seqsetvis::easyLoad_bed(paste0("simulation/genomes/simGenome10M_v4/peaks/simPeaks_v", 1:10, ".bed"))
names(true_gr) = paste0("V", 1:10)
true_dt = lapply(true_gr, as.data.table)
true_dt = rbindlist(true_dt, use.names = TRUE, idcol = "ver")
true_dt[, seqnames := ver]
true_gr = GRanges(true_dt)

olaps = findOverlaps(GRanges(all_qdt), true_gr)
all_qdt$is_true = FALSE
all_qdt$is_true[queryHits(olaps)] = TRUE

all_qdt$is_refined = FALSE
all_qdt[stable_frag_corr > .8 &
            # flex_frag_corr_input < .8 &
            flex_frag_len  > 150 & flex_frag_len < 225,
        is_refined := TRUE]


ggplot(all_qdt[grepl("_10fe", sample)], aes(x = stable_frag_corr ,
                    y = log10(qValue),
                    color = is_true)) +
    geom_point(alpha = .2) + labs(title = "called peaks") + facet_wrap("is_true")

ggplot(all_qdt, aes(x = flex_frag_len,
                    y = log10(qValue),
                    color = is_true)) +
    geom_point() + labs(title = "called peaks") + facet_grid("is_true~is_refined")
# for(i in seq_len(nrow(todo_df))){
#     bam_file = todo_df$bam[i]
#     print(bam_file)
#     if(file.exists(paste0(bam_file, ".bai"))){
#         qgr = rtracklayer::import(todo_df$peak[i], format = "narrowPeak")
#         pipeline(bam_file, qgr, cach_version)
#     }
# }


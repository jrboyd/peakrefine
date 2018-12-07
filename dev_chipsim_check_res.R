#peak call
# for f in simulation/output_simV3/*ed/mac*error; do cmd=$(head $f -n 2 | tail -n 1); cmd=${cmd/"# Command line:"/macs2}; cmd=$cmd" --nomodel --extsize 200"; echo $cmd; $cmd; done

library(seqsetvis)
library(magrittr)
library(GenomicRanges)
library(data.table)
res_dir = "~/R/peakrefine/simulation/genomes/simGenome10M_v5"
ver = 1

all_qgr = lapply(1:10, function(ver){
    peaks = dir(file.path(res_dir, "peaks"), pattern = paste0("simPeaks_v", ver, ".bed"), full.names = TRUE)
    bg = dir(file.path(res_dir, "peaks"), pattern = paste0("simBg_v", ver, ".bed"), full.names = TRUE)
    qgr = rtracklayer::import.bed(peaks)
    qgr
})

sum(lengths(all_qgr))

dt_stats = lapply(1:10, function(ver){
    peaks = dir(file.path(res_dir, "peaks"), pattern = paste0("simPeaks_v", ver, ".bed"), full.names = TRUE)
    # bg = dir(file.path(res_dir, "peaks"), pattern = paste0("simBg_v", ver, ".bed"), full.names = TRUE)
    qgr = rtracklayer::import.bed(peaks)
    qgr = resize(qgr, 1200, fix = "center")
    # bgr = rtracklayer::import.bed(bg)


    # vals = c(qgr$score, bgr$score/20)
    # grps = c(rep("peak", length(qgr)), rep("bg", length(bgr)))
    # boxplot(vals ~ grps)


    seqlevels(qgr) = "chrSim1"
    # bams = dir(file.path(res_dir, "output"), pattern = paste0("V", ver, "_.+R1$"), full.names = TRUE) %>% dir(pattern = ".bam$", full.names = TRUE)
    #
    # nps = dir(file.path(res_dir, "output"), pattern = paste0("V", ver, "_.+R1$"), full.names = TRUE) %>% dir(pattern = "loose.+Peak$", full.names = TRUE)
    # names(nps) = basename(nps) %>% sub("_R1_loose_peaks.narrowPeak", "", .) %>% sub("000reads", "k", .)

    bams = dir(file.path(res_dir, "output"), pattern = paste0("V", ver, "_.+pooled$"), full.names = TRUE) %>%
        dir(pattern = ".bam$", full.names = TRUE)

    nps = dir(file.path(res_dir, "output"), pattern = paste0("V", ver, "_.+pooled$"), full.names = TRUE) %>%
        dir(pattern = "medium.+Peak$", full.names = TRUE)
    names(nps) = basename(nps) %>% sub("peaks.narrowPeak", "", .) %>% sub("000reads", "k", .)

    dt = data.table(name = names(nps))
    dt[, c("ver", "nreads", "fe") := tstrsplit(name, "_", keep = 1:3)]
    dt[, nreads := as.numeric(sub("k", "", nreads))]
    dt[, fe := as.numeric(sub("fe", "", fe))]
    dt = dt[order(fe)][order(nreads)]

    np_grs = easyLoad_narrowPeak(nps[dt$name])
    # np_grs = lapply(np_grs, function(x){
    #     if(length(x) == 0) return(x)
    #     subset(x, qValue > 15)
    # })
    np_grs = append(list(ground = qgr), np_grs)

    # np_grs$true = qgr#subset(qgr, score > 350)
    olaps = ssvOverlapIntervalSets(np_grs)
    print(ssvFeatureBinaryHeatmap(olaps))
    dt_olaps = as.data.table(mcols(olaps))
    dt_olaps = data.table::melt(dt_olaps, id.vars = "ground")
    dt_stats = dt_olaps[, .(total = sum(value),
                            true_pos = sum(ground == TRUE & value == TRUE),
                            false_pos = sum(ground == FALSE & value == TRUE),
                            true_neg = sum(ground == FALSE & value == FALSE),
                            false_neg = sum(ground == TRUE & value == FALSE)), by = .(variable)]

    print(dt_stats$true_pos / length(qgr))

    dt_stats = melt(dt_stats, id.vars = "variable", variable.name = "stat")

    dt_stats[, c("version", "nreads", "fe") := tstrsplit(variable, "_", keep = 1:3)]

    numsortFactor = function(fact){
        fact = factor(fact)
        lev = levels(fact)
        lev = lev[order(as.numeric(gsub("[a-zA-Z]", "", lev)))]
        fact = factor(fact, levels = lev)
        fact
    }



    dt_stats$version = numsortFactor(dt_stats$version)
    dt_stats$nreads = numsortFactor(dt_stats$nreads)
    dt_stats$fe = numsortFactor(dt_stats$fe)



    dt_stats
})

dt_stats = rbindlist(dt_stats)
dt_stats = dt_stats[, .(value = sum(value)), by = .(stat, nreads, fe)]

ggplot(dt_stats[stat != "true_neg" & stat != "false_neg" & stat != "total"],
       aes(x = fe, y = value, color = stat, size = nreads)) +
    geom_point(shape = 21) + theme_classic() + facet_wrap("stat", scales = "free_y")

dt_stats[, key := paste(nreads, fe)]

mdt = merge(dt_stats[stat == "true_pos", 4:5],
      dt_stats[stat == "false_pos", 4:5], by = "key")
colnames(mdt)[2:3] = c("true_pos", "false_pos")

mdt[, c("nreads", "fe") := tstrsplit(key, " ")]

mdt[order(false_pos)][]


# ggplot(dt_stats, aes(x = fe, y = true_pos, color = version, size = nreads)) + geom_point()
# ggplot(dt_stats, aes(x = fe, y = total, color = version, size = nreads)) + geom_point()
# ggplot(dt_stats, aes(x = fe, y = total, color = version, size = nreads)) + geom_point()

# ssvFeatureBinaryHeatmap(olaps)
#
# bam_dt = ssvFetchBam(bams, qgr, target_strand = "both", fragLens = 200, return_data.table = TRUE, win_size = 10)
# bam_dt[, c("version", "depth_str", "fe_str") := tstrsplit(sample, "_", keep = 1:3)]
# bam_dt[, depth := as.numeric(sub("reads", "", depth_str))]
#
# bam_dt$fe_str = factor(bam_dt$fe_str)#, levels = c("input", "50fe", "500fe", "5000fe"))
# lev = levels(bam_dt$fe_str)
# ni = lev != "input"
# o = order(as.numeric(sub("[a-zA-Z]", "", lev[ni])))
# lev = c("input", lev[ni][o])
# bam_dt$fe_str = factor(bam_dt$fe_str, levels = lev)
# # bam_dt$fe_str = factor(bam_dt$fe_str, levels = c("input", "10fe", "50fe", "100fe", "500fe"))
#
# ggplot(bam_dt[id == "region_1"], aes(x = x, y = y, color = strand, group = paste(sample, strand))) +
#     geom_path() + facet_wrap("fe_str~depth", scales = "free_y", ncol = 3) + labs(x = "bp", y = "reads")
#
# # ggplot(bam_dt[id == "region_2"], aes(x = x, y = y, color = strand, group = paste(sample, strand))) +
# #     geom_path() + facet_grid("depth~fe_str") + labs(x = "bp", y = "reads")
#


# pdf(paste0("tmp_", ver, ".pdf"), height = 16, width = 6)
# ssvFeatureBinaryHeatmap(olaps)
# for(i in 1:20){
#     p = ggplot(bam_dt[id == paste0("region_", i)], aes(x = x, y = y, color = strand, group = paste(sample, strand))) +
#         geom_path() + facet_wrap("fe_str~depth", scales = "free_y", ncol = 3) + labs(x = "bp", y = "reads")
#     print(p)
# }
# dev.off()



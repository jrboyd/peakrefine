library(seqsetvis)
library(GenomicRanges)
library(parallel)
library(data.table)
bam_file = "simulation/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled.bam"
peak_file = "simulation/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled_peaks_passIDR.05.narrowPeak"
np_gr = easyLoad_narrowPeak(peak_file)[[1]]
np_gr = resize(np_gr, 600, fix = "center")
names(np_gr) = np_gr$name

st_normal = system.time({
    bam_dt = ssvFetchBam(bam_file, np_gr, target_strand = "both", fragLens = NA, return_data.table = TRUE)
})

ncores = 20
mc_grps = ceiling(seq_along(np_gr) / length(np_gr) * ncores)

mc_gr = split(np_gr, mc_grps)
options("mc.cores" = ncores)

st_mc = system.time({
    bam_list = mclapply(seq_len(ncores), function(i){
        qgr = mc_gr[[i]]
        names(qgr) = qgr$name
        ssvFetchBam(bam_file, qgr, target_strand = "both", fragLens = NA, return_data.table = TRUE)
    })
    bam_dt.mc = rbindlist(bam_list)
})
message("normal :")
st_normal
message("mc :")
st_mc

is_match = all(bam_dt[order(x)][order(id)] ==
        bam_dt.mc[order(x)][order(id)])
stopifnot(is_match)

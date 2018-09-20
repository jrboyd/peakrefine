library(peakrefine)
library(testthat)
library(data.table)
library(GenomicRanges)

bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
bam_input = system.file("extdata", "MCF10A_input.random5.bam", package = "peakrefine")
np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
qgr = rtracklayer::import(np, format = "narrowPeak")
strand(qgr) = c("+", "-", "-", "+", "-")

bam_dt = .fetch_bam_stranded(bam_file, qgr)
GRanges(bam_dt[strand == "+"])
S4Vectors::split(GRanges(bam_dt[strand == "+"]), bam_dt$which_label)
bam_gr = GRanges(coverage(S4Vectors::split(GRanges(bam_dt[strand == "+"]), bam_dt$which_label)))
#test varying size regions
vgr = qgr
end(vgr) = end(vgr) + 0:4*100


test_that("viewGrangeWinSample_dt ids match input", {
    names(qgr) = paste0("peak_", seq_along(qgr))
    sample_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(sample_dt$id)))
})

test_that("viewGrangeWinSample_dt unnamed qgr still creates id", {
    names(qgr) = NULL
    sample_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 100)
    expect_true(!is.null(sample_dt$id))
})

test_that("viewGrangeWinSample_dt ids match input", {
    names(qgr) = paste0("peak_", seq_along(qgr))
    summary_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(summary_dt$id)))

    names(qgr) = seq_along(qgr)
    summary_dt = viewGRangesWinSample_dt(bam_gr, qgr, window_size = 100)
    expect_equal(sort(unique(names(qgr))), sort(unique(summary_dt$id)))
})

test_that("viewGRangesWinSample_dt sizes vary, viewGRangesWinSummary_dt don't", {
    sample_dt = viewGRangesWinSample_dt(bam_gr, vgr, window_size = 100, anchor = "left")
    expect_gt(length(unique(sample_dt[, .N, by = id]$N)), 1)
})

test_that(".remove_duplicates removes duplicates", {
    #create fake read spans where:
    #position 1 has 1 read
    #position 2 has 2 reads
    #position 3 has 3 reads
    #etc.
    dt = data.table(which_label = 1:10, seqnames = "chr1", strand = c("+", "-"), start = 1:10, end = 1:10+10)
    dt[, width := end - start + 1]
    make_dupes = seq_len(nrow(dt))
    make_dupes = rep(make_dupes, make_dupes)
    dt = dt[make_dupes]
    dtl = lapply(1:10, function(x)seqsetvis:::.remove_duplicates(dt, max_dupes = x))
    # max_dupes 1 should yield 10 unique entries
    expect_true(all(dtl[[1]]$which_label == 1:10))
    # max_dupes 10 should perform no dupe removal in this case
    expect_true(all(dtl[[10]]$which_label == make_dupes))
    lens = sapply(dtl, nrow)
    # as max_dupes increases, returned entries increases
    expect_true(all(lens[-length(lens)] - lens[-1] < 0))
})

test_that(".fetch_bam_stranded basic", {
    bam_dt = .fetch_bam_stranded(bam_file, qgr)
    expect_equal(length(unique(bam_dt$which_label)), length(qgr))

})

test_that("strandsCoverage basic", {
    bam_dt = .fetch_bam_stranded(bam_file, qgr)
    strandsCoverage(bam_dt, qgr)
})

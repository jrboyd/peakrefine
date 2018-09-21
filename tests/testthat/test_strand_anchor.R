# flipping viewGranges
library(peakrefine)
library(GenomicRanges)
library(testthat)
library(data.table)

bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
bam_input = system.file("extdata", "MCF10A_input.random5.bam", package = "peakrefine")
np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
qgr = rtracklayer::import(np, format = "narrowPeak")
strand(qgr) = c("+", "-", "-", "+", "-")

reads_dt = peakrefine:::.fetch_bam_stranded(bam_file, qgr)
# bam_dt = .calc_stranded_coverage(reads_dt, qgr)
# bam_gr = GRanges(bam_dt)
# bam_gr = bam_gr[strand(bam_gr) == "+"]
# bam_gr$score = bam_gr$y

ext_cov = GenomicRanges::coverage(GenomicRanges::split(
    GenomicRanges::GRanges(reads_dt[strand == "+"]), reads_dt$which_label))
score_gr = GenomicRanges::GRanges(ext_cov)

view_gr_by_window_sample = peakrefine:::.view_gr_by_window_sample

bam_dt = view_gr_by_window_sample(score_gr, qgr, 1,
                                 anchor = "center_unstranded")
bam_gr = GRanges(bam_dt)
bam_gr$score = bam_gr$y
shift_anchor = peakrefine:::.shift_anchor

test_that(".shift_anchor center", {
    s_dt = shift_anchor(bam_dt, 1, "center")
    expect_lt(min(s_dt$x), 0)
    expect_gt(max(s_dt$x), 0)
    expect_equal(min(s_dt$x), -max(s_dt$x))
})

test_that(".shift_anchor center_unstranded", {
    s_dt = shift_anchor(bam_dt, 1, "center_unstranded")
    expect_lt(min(s_dt$x), 0)
    expect_gt(max(s_dt$x), 0)
    expect_equal(min(s_dt$x), -max(s_dt$x))
})

test_that(".shift_anchor center vs center_unstranded", {
    c_dt = copy(shift_anchor(bam_dt, 1, "center"))
    cu_dt = shift_anchor(bam_dt, 1, "center_unstranded")

    expect_true(all(c_dt[strand == "+"] == cu_dt[strand == "+"]))
    expect_true(!all(c_dt[strand == "-"] == cu_dt[strand == "-"]))
})


test_that(".shift_anchor left", {
    s_dt = shift_anchor(bam_dt, 1, "left")
    expect_equal(min(s_dt$x), 0)
    expect_gt(max(s_dt$x), 0)
})

# test_that(".shift_anchor left_unstranded", {
#     s_dt = shift_anchor(bam_dt, 1, "left_unstranded")
# })
#
# test_that(".shift_anchor left vs left_unstranded", {
#     c_dt = copy(shift_anchor(bam_dt, 1, "left"))
#     cu_dt = shift_anchor(bam_dt, 1, "left_unstranded")
#
#     expect_true(all(c_dt[strand == "+"] == cu_dt[strand == "+"]))
#     expect_true(!all(c_dt[strand == "-"] == cu_dt[strand == "-"]))
# })

# #sampling
# test_that(".view_gr_by_window_sample center", {
#     dt = peakrefine:::.view_gr_by_window_sample(bam_gr, qgr, window_size = 50, anchor = "center")
#     # dt$group = "stranded"
#     dt_us = peakrefine:::.view_gr_by_window_sample(bam_gr, qgr, window_size = 50, anchor = "center_unstranded")
#     # dt_us$group = "unstranded"
#     # ssvSignalLineplot(rbind(dt, dt_us), sample_ = "id", color_ = "strand", group_ = "group")
#     expect_true(all(dt[strand == "+"][order(x)]$y == dt_us[strand == "+"][order(x)]$y))
#     expect_false(all(dt[strand == "-"][order(x)]$y == dt_us[strand == "-"][order(x)]$y))
#     expect_false(all(dt$x > 0))
#     expect_false(all(dt$x < 0))
#     expect_false(all(dt_us$x > 0))
#     expect_false(all(dt_us$x < 0))
# })
#
# test_that(".view_gr_by_window_sample left", {
#     dt = peakrefine:::.view_gr_by_window_sample(bam_gr, qgr, window_size = 50, anchor = "left")
#     # dt$group = "stranded"
#     dt_us = peakrefine:::.view_gr_by_window_sample(bam_gr, qgr, window_size = 50, anchor = "left_unstranded")
#     # dt_us$group = "unstranded"
#     # ssvSignalLineplot(rbind(dt, dt_us), sample_ = "id", color_ = "strand", group_ = "group")
#     expect_equal(dt[strand == "+"][order(x)]$y, dt_us[strand == "+"][order(x)]$y)
#     expect_failure(expect_equal(dt[strand == "-"][order(x)]$y, dt_us[strand == "-"][order(x)]$y))
#     expect_true(all(dt$x > 0))
#     expect_true(all(dt_us$x > 0))
# })


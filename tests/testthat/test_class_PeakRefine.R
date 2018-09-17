# flipping viewGranges
library(seqsetvis)
library(GenomicRanges)
library(testthat)
library(data.table)

bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
bam_input = system.file("extdata", "MCF10A_input.random5.bam", package = "peakrefine")
np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
np = rtracklayer::import(np, format = "narrowPeak")

#sampling
test_that("PeakRefiner init", {
    pr = new("PeakRefiner", np, bam_file, bam_input, fragment_lengths = c(50, 200), pwm = PWMLogn.hg19.MotifDb.Hsap)
    expect_s4_class(pr, "PeakRefiner")
    expect_equal(pr@bam_treat_file, bam_file)
    expect_equal(pr@bam_input_file, bam_input)
    expect_equal(pr@fragment_lengths, c(50, 200))
    expect_equal(pr@fl2color, rep("black", 2))
    expect_equal(pr@output_prefix, sub("\\.bam", "", basename(bam_file)))
})

test_that("PeakRefiner constructor", {
    pr = PeakRefiner(np, bam_file, bam_input, fragment_lengths = c(50, 200), pwm = PWMLogn.hg19.MotifDb.Hsap)
    expect_s4_class(pr, "PeakRefiner")
    expect_equal(pr@bam_treat_file, bam_file)
    expect_equal(pr@bam_input_file, bam_input)
    expect_equal(pr@fragment_lengths, c(50, 200))
    expect_equal(pr@fl2color, rep("black", 2))
    expect_equal(pr@output_prefix, sub("\\.bam", "", basename(bam_file)))
})

test_that("PeakRefiner auto frag length", {
    pr = PeakRefiner(np, bam_file, bam_input, pwm = PWMLogn.hg19.MotifDb.Hsap,
                     auto_frag_len_FUN = function(bf, np){
                         ssvStrandCorr(bf, np, step = 50, small_step = 25)
    })
    expect_equal(pr@fragment_lengths, c(101, 195))
    expect_s3_class(pr@plots[[1]], "ggplot")
})


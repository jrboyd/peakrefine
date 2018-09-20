library(peakrefine)
library(GenomicRanges)
library(testthat)
library(data.table)

bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
bam_input = system.file("extdata", "MCF10A_input.random5.bam", package = "peakrefine")
np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
np = rtracklayer::import(np, format = "narrowPeak")

# test_that("crossCorrByExtension", {
#     crossCorrByExtension(bam_file, np, step = 50, small_step = 25)
# })

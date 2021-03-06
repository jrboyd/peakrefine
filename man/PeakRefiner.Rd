% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_peakrefine.R
\name{PeakRefiner}
\alias{PeakRefiner}
\title{constructor for PeakRefiner}
\usage{
PeakRefiner(peak_set, bam_treat_file, bam_input_file, pwm,
  target_pwm_names, fragment_lengths = NULL, color_overrides = NULL,
  color_default = "black", auto_frag_len_FUN = NULL,
  output_prefix = NULL)
}
\arguments{
\item{peak_set}{GRanges formatted peak set to calculate motifs on}

\item{bam_treat_file}{.bam file of aligned reads from ChIP-seq pulldown.
Index must be at .bam.bai}

\item{bam_input_file}{.bam file of aligned reads from input control. Index
must be at .bam.bai}

\item{pwm}{Position Weight Matrix from PWMEnrich}

\item{target_pwm_names}{names of PWMs to include in figures by default}

\item{fragment_lengths}{fragment lengths to consider}

\item{color_overrides}{character. colors to use for each item in
fragment_lengths.  should contain valid hex ("#000000") or R colors
("black") and be named with items in fragment_lengths.  non-overriden items
will be black.}

\item{color_default}{single character. black.}

\item{auto_frag_len_FUN}{function to use to auto calculate fragment length.
must accept two argument, bam_file and peak_set.}

\item{output_prefix}{prefix to use for output files.}
}
\value{
a new valid object of class PeakRefiner
}
\description{
constructor for PeakRefiner
}
\examples{
bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
bam_input = system.file("extdata", "MCF10A_input.random5.bam", package = "peakrefine")
np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
qgr = rtracklayer::import(np, format = "narrowPeak")

library(PWMEnrich.Hsapiens.background)
data("PWMLogn.hg19.MotifDb.Hsap")
pwm = PWMLogn.hg19.MotifDb.Hsap[1:2]
PeakRefiner(qgr, bam_file, bam_input, pwm, names(pwm$pwms)[1])
}

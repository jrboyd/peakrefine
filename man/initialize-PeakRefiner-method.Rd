% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_peakrefine.R
\docType{methods}
\name{initialize,PeakRefiner-method}
\alias{initialize,PeakRefiner-method}
\title{initialize a new PeakRefiner}
\usage{
\S4method{initialize}{PeakRefiner}(.Object, peak_set, bam_treat_file,
  bam_input_file, pwm, target_pwm_names, fragment_lengths = NULL,
  color_overrides = NULL, color_default = "black",
  auto_frag_len_FUN = NULL, output_prefix = NULL)
}
\arguments{
\item{.Object}{empty PeakRefiner}

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
valid PeakRefiner
}
\description{
initialize a new PeakRefiner
}

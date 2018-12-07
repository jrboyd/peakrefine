# peakrefine
refining peaks using metrics derived from aligned read strand information

Installation:
devtools::install_github("jrboyd/peakrefine")

Running:
You just need an indexed bam file and a set of peaks in a GRanges.
peakrefine::calcCorrMetrics(bam_file, peaks_gr, frag_min = 50, frag_max = 300)

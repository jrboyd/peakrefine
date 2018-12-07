peakrefine
================

# peakrefine

Refining peaks using metrics derived from Strand Cross-Correlation (SCC)

devtools::install\_github(“jrboyd/peakrefine”)

## Installation:

``` r
if(!require(devtools))
    install.packages("devtools")
devtools::install_github("jrboyd/peakrefine")
```

## Usage

``` r
bam_file = "path/to/your/bam" #must be indexed
peaks_gr = rtracklayer::import("path/to/your/peaks")
peakrefine::calcCorrMetrics(bam_file, peaks_gr, frag_min = 50, frag_max = 300)
```

## GitHub Documents

This is an R Markdown format used for publishing markdown documents to
GitHub. When you click the **Knit** button all R code chunks are run and
a markdown file (.md) suitable for publishing to GitHub is generated.

## Including Code

You can include R code in the document as follows:

``` r
summary(cars)
```

## Including Plots

You can also embed plots, for example:

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.

peakrefine
================

Refining peaks using metrics derived from Strand Cross-Correlation (SCC)

## Installation:

``` r
if(!require(devtools))
    install.packages("devtools")
devtools::install_github("jrboyd/peakrefine")
```

## Basic Usage

``` r
bam_file = "path/to/your/bam" #must be indexed
peaks_gr = rtracklayer::import("path/to/your/peaks")
corr_res = peakrefine::calcSCCMetrics(bam_file, peaks_gr, frag_sizes = 50:250)
```

Run this way, 100 peaks will take a bit under a minute.

Results will be automatically cached using BiocFileCache such that
repeated calls with identical parameters load previous results.

## Running in parallel

peakrefine uses the mc.cores option to automatically split up peaks for
calculations. If you run into memory or similar issues you should set
n\_splits to a multiple of mc.cores (this applies to running in serial
as well).  
The default is 1 job per core which may prove too large.

``` r
ncores = max(1, parallel::detectCores() - 1)
options(mc.cores = ncores)
peakrefine::calcSCCMetrics(bam_file, peaks_gr, frag_sizes = 50:250, n_splits = ncores * 3)
```

## Using the output

Output is a named list:

``` r
cbind(names(corr_res))
```

    ##      [,1]                         
    ## [1,] "read_length"                
    ## [2,] "fragment_length"            
    ## [3,] "read_correlation"           
    ## [4,] "flex_fragment_correlation"  
    ## [5,] "stable_fragment_correlation"
    ## [6,] "full_correlation_results"

Getting read and fragment lengths:

``` r
corr_res$read_length
```

    ## [1] 65

``` r
corr_res$fragment_length
```

    ## [1] 162

Other list items are data.tables:

``` r
sapply(corr_res, class)
```

    ## $read_length
    ## [1] "numeric"
    ## 
    ## $fragment_length
    ## [1] "numeric"
    ## 
    ## $read_correlation
    ## [1] "data.table" "data.frame"
    ## 
    ## $flex_fragment_correlation
    ## [1] "data.table" "data.frame"
    ## 
    ## $stable_fragment_correlation
    ## [1] "data.table" "data.frame"
    ## 
    ## $full_correlation_results
    ## [1] "data.table" "data.frame"

We’ll need a couple more libraries here.

``` r
library(data.table)
library(ggplot2)
```

``` r
metrics_dt = rbindlist(corr_res[3:5], use.names = TRUE, idcol = "metric")
theme_set(theme_classic())
ggplot(metrics_dt, aes(x = metric, y = correlation, color = metric)) + geom_boxplot()
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Stable fragment correlation uses the calculated fragment size while flex
always uses the fragment size with the maximum SCC. Flex will report
high correlation from artifact peaks but is useful for assessing
fragment size
distribution.

``` r
ggplot(metrics_dt, aes(x = metric, y = shift, color = metric)) + geom_boxplot()
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

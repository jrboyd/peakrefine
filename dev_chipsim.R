### R code from vignette source 'ChIPsimIntro.Rnw'

library(ChIPsim)
library(BiocFileCache)
library(digest)
library(data.table)
source("dev_chipsim_functions.R")

bfcif = peakrefine:::bfcif
bfc_sim = BiocFileCache("~/.cache_sim")
set.seed(1)

#size of genome
gen_size = 100e6
read_depths = c(5e4, 10e4, 50e4, 100e4)
chip_enrichments = c(2, 5, 10, 50)

n_reads = read_depths[1]
f_enrich = chip_enrichments[1]



make_genome = function(){
    Biostrings::DNAStringSet(c(CHR=paste(sample(Biostrings::DNA_BASES, gen_size, replace = TRUE), collapse = "")))
}

genome <- bfcif(bfc_sim, rname = digest(list("genome", gen_size, "v1")), make_genome)
names(genome) = "chrSim1"
Biostrings::writeXStringSet(genome, file = "simulation/simGenome.fa")


ver = 2
fo = FALSE
options(scipen=999)
all_sim = list()
for(n_reads in read_depths){
    for(f_enrich in chip_enrichments){

        this_make_sim = function(){
            make_sim(genome, n_reads, f_enrich, seed = ver)
        }
        message(n_reads, " ", f_enrich)
        sim = bfcif(bfc = bfc_sim, digest(list("sim", genome, n_reads, f_enrich, paste0("v", ver))), this_make_sim, force_overwrite = fo)
        sim_fastq = file.path("simulation/fastqs", paste0("simV", ver, "_", n_reads, "_", f_enrich, ".fastq"))
        # if(!file.exists(sim_fastq)){
            message("writing fastq")
            my_writeFASTQ.dt(file = sim_fastq,
                          read = sim$readSequence$sequence, quality = sim$readSequence$quality, name = sim$readSequence$name)
        # }

        all_sim[[paste(n_reads, f_enrich)]] = sim
    }

    this_make_input = function(){
        make_sim(genome, n_reads, f_enrich, seed = ver, bind_p = 0)
    }
    message(n_reads, " ", "input")
    sim = bfcif(bfc = bfc_sim, digest(list("sim", genome, n_reads, "input", paste0("v", ver))), this_make_input, force_overwrite = TRUE)
    sim_fastq = file.path("simulation/fastqs", paste0("simV", ver, "_", n_reads, "_", "input", ".fastq"))
    # if(!file.exists(sim_fastq)){
        message("writing fastq")
        my_writeFASTQ.dt(file = sim_fastq,
                      read = sim$readSequence$sequence, quality = sim$readSequence$quality, name = sim$readSequence$name)
    # }
}

all_simPeaks = lapply(all_sim, function(sim){
    x = sim$features
    k = sapply(x[[1]], function(y)class(y)[1]) == "Binding"
    binders = x[[1]][which(k)]
    starts = sapply(binders, function(b){
        b$start
    })
    widths = sapply(binders, function(b){
        b$length
    })
    weights = sapply(binders, function(b){
        b$weight
    })
    library(GenomicRanges)
    simPeaks = GRanges("chrSim", IRanges(starts, starts + widths), weight = weights)
    simPeaks
})

#check if all starts equal
all_starts = lapply(all_simPeaks, function(x)start(x))
stopifnot(all(unique(unlist(all_starts)) == all_starts$`50000 2`))
#check if weights correlate perfectly
mat_weight = sapply(all_simPeaks, function(x)x$weight)
stopifnot(all(round(cor(mat_weight), digits = 5) == 1))

sim_peak = as.data.table(all_simPeaks$`50000 2`)
sim_peak$weight = apply(mat_weight, 1, max)
sim_peak[, weight := weight / min(weight)]
sim_peak$width = NULL
sim_peak[, name := paste0("true_", seq(.N))]
sim_peak[, weight := round(weight * 100)]
sim_peak = sim_peak[, c(1:3, 6, 5, 4)]
write.table(sim_peak, file.path("simulation/peaks", paste0("simPeaks_v", ver, ".bed")), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


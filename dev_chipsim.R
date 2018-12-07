### R code from vignette source 'ChIPsimIntro.Rnw'
setwd("~/R/peakrefine")
library(ChIPsim)
library(BiocFileCache)
library(digest)
library(data.table)
source("dev_chipsim_functions.R")
options(scipen=999)
#ver
# 8 use system time to seed reads

ver = 3
ver_g = 5
bfcif = peakrefine:::bfcif
bfc_sim = BiocFileCache("~/.cache_sim")
set.seed(ver)

make_genome = function(){
    message("making genome ", gen_size, " bp")
    set.seed(ver_g)
    Biostrings::DNAStringSet(c(CHR=paste(sample(Biostrings::DNA_BASES, gen_size, replace = TRUE), collapse = "")))
}

gen_size = 10e6
genome_name = paste0("10M_v", ver_g)
genome <- bfcif(bfc_sim, rname = digest(list("genome", gen_size, genome_name, ver_g)), make_genome)
names(genome) = paste0("chrSim", 1)

gen_dir = paste0("simulation/genomes/simGenome", genome_name)
dir.create(gen_dir, showWarnings = FALSE, recursive = TRUE)
gen_dir = normalizePath(gen_dir)
dir.create(file.path(gen_dir, "fastqs"), showWarnings = FALSE)
dir.create(file.path(gen_dir, "peaks"), showWarnings = FALSE)
genome_file = file.path(gen_dir, paste0(basename(gen_dir), ".fa"))
Biostrings::writeXStringSet(genome, file = genome_file)
genome_file = normalizePath(genome_file)

for(ver in 1:10){
    message("making sim version ", ver)


    #size of genome


    file.copy("dev_chipsim.R", paste0("dev_chipsim_V", ver, "_", genome_name, ".R"))

    read_depths = c(50e4)
    chip_enrichments = c(1500)
    fraction_enrich = c(0, .005, .01, .1, 1)



    if(!dir.exists(file.path("/slipstream/galaxy/data", basename(gen_dir))))
        system(paste("ln -s", normalizePath(gen_dir), "/slipstream/galaxy/data"))

    idx_dir =  file.path(gen_dir, "star_index")
    if(!dir.exists(idx_dir)){
        dir.create(idx_dir, recursive = TRUE, showWarnings = FALSE)
        setwd(idx_dir)
        system(paste("bash ~/scripts/make_STAR_indexes_small.sh", genome_file))
        setwd("~/R/peakrefine/")
    }


    bfcif.path = function(bfc, rname){
        # is rname in cache?
        if(nrow(bfcquery(bfc, query = rname, field = "rname")) == 0){
            cache_path = bfcnew(bfc, rname = rname)

        }else{
            cache_path = bfcrpath(bfc, rname)
        }
        cache_path
    }


    fo = FALSE
    message("allocate keys...")
    todo_df = as.data.frame(expand.grid(read_depths, chip_enrichments, fraction_enrich))
    todo_df$path = sapply(seq_len(nrow(todo_df)), function(i){
        n_reads = round(todo_df[i, 1] * todo_df[i, 3])
        f_enrich = todo_df[i, 2]
        key = digest(list("sim", genome, n_reads, f_enrich, paste0("v", ver)))
        bfcif.path(bfc = bfc_sim, digest(list("sim", genome, n_reads, f_enrich, paste0("v", ver))))
    })

    todo_df$path_bg = sapply(seq_len(nrow(todo_df)), function(i){
        n_reads = todo_df[i, 1] - round(todo_df[i, 1] * todo_df[i, 3])
        f_enrich = todo_df[i, 2]
        key = digest(list("simBg", genome, n_reads, f_enrich, paste0("v", ver)))
        bfcif.path(bfc = bfc_sim, digest(list("simBg", genome, n_reads, f_enrich, paste0("v", ver))))
    })

    # options(mc.cores = min(nrow(todo_df), 8))
    options(mc.cores = nrow(todo_df))

    all_sim = parallel::mclapply(seq_len(nrow(todo_df)), function(i){
    #all_sim = lapply(seq_len(nrow(todo_df)), function(i){
        n_reads = round(todo_df[i, 1] * todo_df[i, 3])
        f_enrich = todo_df[i, 2]
        fraction_enrich = todo_df[i, 3]
        path = todo_df$path[i]
        if(file.exists(path)){
            load(path)
        }else{
            this_make_sim = function(){
                if(n_reads == 0){
                    list(readSequence =
                             list(
                                 sequence = character(),
                                 quality = character(),
                                 name = character()
                             )

                    )
                }else{
                    make_sim(genome, n_reads, f_enrich, seed = ver)
                }

            }
            message(n_reads, " ", f_enrich)
            sim = this_make_sim()
            save(sim, file = path)
        }

        path_bg = todo_df$path_bg[i]
        if(file.exists(path_bg)){
            load(path_bg)
        }else{
            n_reads_bg = todo_df[i, 1] - n_reads
            this_make_sim_bg = function(){
                if(n_reads_bg == 0){
                    list(
                        read = character(),
                        quality = character(),
                        name = character()
                    )
                }else{
                    make_unif_sim(genome, n_reads_bg)
                }
            }
            message(n_reads_bg, " ", "bg")
            sim_bg = this_make_sim_bg()
            save(sim_bg, file = path_bg)
        }
        # browser()
        if(fraction_enrich == 0){
            sim_fastq = file.path(gen_dir, "fastqs", paste0("simV", ver, "_", todo_df[i, 1] , "_input.fastq"))
        }else{
            sim_fastq = file.path(gen_dir, "fastqs", paste0("simV", ver, "_", todo_df[i, 1] , "_", 100*fraction_enrich, ".fastq"))
        }
        message("writing fastq ", sim_fastq)
        my_writeFASTQ.dt(file = sim_fastq,
                         read = c(sim$readSequence$sequence, sim_bg$read),
                         quality = c(sim$readSequence$quality, sim_bg$quality),
                         name = c(sub("read_", "signal_", sim$readSequence$name), sub("read_", "bg_", sim_bg$name)))
        sim
    })
    names(all_sim) = paste(todo_df$Var1, todo_df$Var2)

    all_simPeaks = lapply(all_sim[-1], function(sim){
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
        simPeaks = GRanges("chrSim1", IRanges(starts, starts + widths), weight = weights)
        simPeaks
    })

    #check if all starts equal
    all_starts = lapply(all_simPeaks, function(x)start(x))
    stopifnot(all(unique(unlist(all_starts)) == all_starts[[1]]))
    #check if weights correlate perfectly
    mat_weight = sapply(all_simPeaks, function(x)x$weight)
    stopifnot(all(round(cor(mat_weight), digits = 5) == 1))

    sim_peak = as.data.table(all_simPeaks[[1]])
    sim_peak$weight = apply(mat_weight, 1, max)
    sim_peak[, weight := weight / min(weight)]
    sim_peak$width = NULL
    sim_peak[, name := paste0("true_", seq(.N))]
    sim_peak[, weight := round(weight * 100)]
    sim_peak = sim_peak[, c(1:3, 6, 5, 4)]
    write.table(sim_peak, file.path(gen_dir, "peaks", paste0("simPeaks_v", ver, ".bed")), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

    # all_simPeaks = lapply(all_sim[-1], function(sim){
    #     x = sim$features
    #     k = sapply(x[[1]], function(y)class(y)[1])== "Binding"
    #     binders = x[[1]][which(k)]
    #     starts = sapply(binders, function(b){
    #         b$start
    #     })
    #     widths = sapply(binders, function(b){
    #         b$length
    #     })
    #     weights = sapply(binders, function(b){
    #         b$weight
    #     })
    #     browser()
    #     library(GenomicRanges)
    #     simPeaks = GRanges("chrSim1", IRanges(starts, starts + widths), weight = weights)
    #     simPeaks
    # })
    #
    # #check if all starts equal
    # all_starts = lapply(all_simPeaks, function(x)start(x))
    # stopifnot(all(unique(unlist(all_starts)) == all_starts$`50000 2`))
    # #check if weights correlate perfectly
    # mat_weight = sapply(all_simPeaks, function(x)x$weight)
    # stopifnot(all(round(cor(mat_weight), digits = 5) == 1))
    #
    # sim_peak = as.data.table(all_simPeaks[[1]])
    # sim_peak$weight = apply(mat_weight, 1, max)
    # sim_peak[, weight := weight / min(weight)]
    # sim_peak$width = NULL
    # sim_peak[, name := paste0("true_", seq(.N))]
    # sim_peak[, weight := round(weight * 100)]
    # sim_peak = sim_peak[, c(1:3, 6, 5, 4)]
    # write.table(sim_peak, file.path(gen_dir, "peaks", paste0("simBg_v", ver, ".bed")), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
}

fastq_dir = file.path(gen_dir, "fastqs")
qc_file = file.path(gen_dir, "fastqs", "qc_config.csv")
in_dir = paste0("in=", normalizePath(fastq_dir))
out_dir = paste0("out=", file.path(normalizePath(gen_dir), paste0("output")))
if(dir.exists(out_dir)) unlink(out_dir)
write.table(rbind(in_dir, out_dir), file = qc_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
qc_file = normalizePath(qc_file)


qc_dt = data.table(fastq = dir(fastq_dir, pattern = ".+fastq$"))#pattern = paste0("V", ver, ".+fastq$")))
qc_dt[, c("cell", "nreads", "fe") := tstrsplit(fastq, "[_\\.]", keep = 1:3)]
qc_dt$cell = sub("sim", "", qc_dt$cell)
qc_dt[, nreads := paste0(nreads, "reads")]
qc_dt[fe != "input", fe := paste0(fe, "fe")]
qc_dt$rep = "R1"
fwrite(qc_dt, file = qc_file, append = TRUE)

dir.create(file.path(file.path(normalizePath(gen_dir), paste0("output"))), showWarnings = FALSE)

setwd("~/R/peakrefine/simulation/qc_framework")
cmd = paste("bash qc_framework.sh", qc_file, basename(gen_dir), " >", file.path(file.path(normalizePath(gen_dir), paste0("output")), "submit.log"))
print(cmd)
# system()
setwd("~/R/peakrefine")

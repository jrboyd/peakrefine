library(magrittr)
pipe_output = "/slipstream/home/joeboyd/R/peakrefine/simulation/genomes/simGenome10M_v5/output"
all_bams = dir(pipe_output, full.names = TRUE, pattern = "pooled$") %>% dir(pattern = "bam$", full.names = TRUE)

base_cmd = "macs2 callpeak -t TREAT -c CTRL -g 10000000 --outdir OUTDIR -n NAME -p PVAL --bdg --nomodel --extsize 200"

for(treat_bam in all_bams[!grepl("input", all_bams)]){
    message(treat_bam)
    pval = ".01"
    fe = strsplit(basename(treat_bam), "_")[[1]][3]
    ctrl_bam = gsub(fe, "input", treat_bam)
    stopifnot(file.exists(treat_bam))
    stopifnot(file.exists(ctrl_bam))
    cmd = base_cmd
    cmd = sub("TREAT", treat_bam, cmd)
    cmd = sub("CTRL", ctrl_bam, cmd)
    cmd = sub("OUTDIR", dirname(treat_bam), cmd)
    cmd = sub("NAME", paste0(sub(".bam", "", basename(treat_bam)), "_mediumLoose"), cmd)
    cmd = sub("PVAL", pval, cmd)
    hide = system(cmd, intern = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE)
}

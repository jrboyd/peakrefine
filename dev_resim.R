# chipsim has strand offset background which isn't ok
# combine chip-sim signal with randomized background (random position and random read orientation)
source("dev_chipsim_functions.R")
library(GenomicRanges)
library(Biostrings)
library(pbapply)
# gr = GRanges(gsub(",", "", "chr1:119,000,000-123,000,000"))
# gr =GRanges(gsub(",", "", "chr1:142,500,000-142,511,000"))
# library(BSgenome.Hsapiens.UCSC.hg38)
# mySeq = getSeq(Hsapiens, gr)

#' makes a simulated set of reads for input fasta
#' every position is equally likely to be output
#'
#' @param gen_fasta path to fasta file
#' @param nreads number of reads to return
#' @param readSize size of reads to return
#' @param plusStrandRatio ratio of plus to minus reads. 1 means 1:1 - an even split. 2 means 2:1 - 2/3 of reads will be plus.
#'
#' @return list of read info
#' read - sequence
#' name - name of reads
#' qual - quality - default is I.
#' @export
#'
#' @examples
make_unif_sim = function(
    gen_fasta,
    nreads,
    readSize = 50,
    qual_char = "I",
    plusStrandRatio = 1){

    strand_cut = (plusStrandRatio) / (1 + plusStrandRatio)
    mySeq = readDNAStringSet("simulation/genomes/simGenome10M_v3/simGenome10M_v3.fa")

    myComp = complement(mySeq)

    mySeq = sub("^N+", "", mySeq)
    myComp = sub("^N+", "", myComp)

    mySplit = strsplit(mySeq, "")[[1]]
    mySplitComp = strsplit(myComp, "")[[1]]

    todo_df = data.frame("index" = sample(MAX, size = nreads, replace = TRUE),
                         "strand" = sample(c("+","-"), nreads, prob = c(plusStrandRatio, 1), replace = TRUE))

    MAX = nchar(mySeq) - readSize + 1

    tmp = matrix(unlist(pblapply(seq_len(nrow(todo_df)), function(i){#~5s per 1k
        itr = todo_df[i, "index"]
        rn = paste0("read_", i)
        if(todo_df[i, "strand"] == "+"){
            seq = paste(x = mySplit[itr:(itr + readSize  - 1)], collapse = "")
        }else{
            seq = paste(x = mySplitComp[itr:(itr + readSize  - 1)], collapse = "")
        }
        c(rn, seq)
    })), ncol = 2, byrow = TRUE)

    qual = paste(rep(qual_char, nchar(tmp[1,2])), collapse = "")

    list(read = tmp[,2], name = tmp[,1], quality = rep(qual, nrow(tmp)))
}
# my_writeFASTQ.dt(read = tmp[,2], name = tmp[,1], quality = rep(qual, nrow(tmp)), file = paste0("random_", nreads, ".fastq"))

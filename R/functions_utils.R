#' Ensure that regions in GRanges are compatible wtih seqlengths of bam_file
#'
#' @param gr GRanges that possibly extend to invalid start or end values.
#' @param bam_file bam_file to retrieve seqlengths info from
#'
#' @return GRanges shifted such that all start at 1 or greater and end prior to chromosome end
#' @export
#'
#' @examples
harmonize_seqlengths = function(gr, bam_file){
    chr_lengths = scanBamHeader(bam_file)[[1]]$targets
    seqlengths(gr) = chr_lengths[names(seqlengths(gr))]
    too_long = end(gr) > seqlengths(gr)[as.character(seqnames(gr))]
    if(any(too_long)){
        message(sum(too_long), " region shifted for extending beyond seqlengths")
        fix_gr = gr[too_long]
        shift_by = -(end(fix_gr) - seqlengths(fix_gr)[as.character(seqnames(fix_gr))])
        gr[too_long] = GenomicRanges::shift(fix_gr, shift_by)
    }
    too_short = start(gr) < 1
    if(any(too_short)){
        message(sum(too_short), " region shifted for starting before seqlengths")
        fix_gr = gr[too_short]
        shift_by = 1 - start(fix_gr)
        gr[too_short] = GenomicRanges::shift(fix_gr, shift_by)
    }
    gr
}

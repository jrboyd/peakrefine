#S4 class to load and refine peak sets evaluated by motif enrichment
# validate config
# use cached data as possible

.datatable.aware = TRUE

#' class to handle peak refinement
#'
#' @slot peak_set GRanges.
#' @slot bam_treat_file character.
#' @slot bam_input_file character.
#' @slot fragment_lengths integer.
#' @slot fl2color character.
#' @slot pwm PWMLognBackground.
#' @slot target_pwm_names character.
#' @slot output_prefix character.
#' @slot plots list.
#'
#' @export
#' @importClassesFrom  PWMEnrich PWMLognBackground
#' @importClassesFrom  GenomicRanges GRanges
#' @import PWMEnrich
setClass(Class = "PeakRefiner",

         slots = c(
             peak_set = "GRanges",
             bam_treat_file = "character",
             bam_input_file = "character",
             fragment_lengths = "integer",
             fl2color = "character",
             pwm = "PWMLognBackground",
             target_pwm_names = "character",
             output_prefix = "character",
             plots = "list"
         ),

         validity = function(object){
             errors <- character()
             # mat_cnames = c("i", "j", "val")
             # if (length(intersect(colnames(object@hic_2d), mat_cnames)) != length(mat_cnames)){
             #     msg <- "colnames of hic_2d must be c(i, j, val)"
             #     errors <- c(errors, msg)
             # }
             # reg_cnames = c("seqnames", "start", "end", "index")
             # if (length(intersect(colnames(object@hic_1d), reg_cnames)) != length(reg_cnames)){
             #     msg <- "colnames of hic_1d must be c(seqnames, start, end, index)"
             #     errors <- c(errors, msg)
             # }
             if (length(errors) == 0) TRUE else errors
         }
)

#' initialize a new PeakRefiner
#'
#' @param .Object empty PeakRefiner
#' @param peak_set GRanges formatted peak set to calculate motifs on
#' @param bam_treat_file .bam file of aligned reads from ChIP-seq pulldown.
#'   Index must be at .bam.bai
#' @param bam_input_file .bam file of aligned reads from input control. Index
#'   must be at .bam.bai
#' @param pwm Position Weight Matrix from PWMEnrich
#' @param target_pwm_names names of PWMs to include in figures by default
#' @param fragment_lengths fragment lengths to consider
#' @param color_overrides character. colors to use for each item in
#'   fragment_lengths.  should contain valid hex ("#000000") or R colors
#'   ("black") and be named with items in fragment_lengths.  non-overriden items
#'   will be black.
#' @param color_default single character. black.
#' @param auto_frag_len_FUN function to use to auto calculate fragment length.
#'   must accept two argument, bam_file and peak_set.
#' @param output_prefix prefix to use for output files.
#'
#' @return valid PeakRefiner
#'
#' @importFrom methods validObject
setMethod("initialize", "PeakRefiner", function(.Object,
                                                peak_set,
                                                bam_treat_file,
                                                bam_input_file,
                                                pwm,
                                                target_pwm_names,
                                                fragment_lengths = NULL,
                                                color_overrides = NULL,
                                                color_default = "black",
                                                auto_frag_len_FUN = NULL,
                                                output_prefix = NULL) {
    # if(missing(matrix_file) & missing(regions_file) & missing(parameters)){
    #     return(.Object)
    # }
    if(is.null(peak_set$id)){
        if(is.null(names(peak_set))){
            peak_set$id = paste0("peak_", seq_along(peak_set))
        }else{
            peak_set$id = names(peak_set)
        }
    }

    .Object@peak_set = peak_set
    .Object@bam_treat_file = bam_treat_file
    .Object@bam_input_file = bam_input_file
    .Object@pwm = pwm
    .Object@target_pwm_names = target_pwm_names
    plots = list()

    if(is.null(fragment_lengths)){
        warning("fragment_lengths not set,  attempting to determine read and fragment sizes...")
        # browser()
        if(is.null(auto_frag_len_FUN)){
            # auto_frag_len_FUN = crossCorrByExtension()
            warning("TODO replace auto FUN")
            auto_frag_len_FUN = function(bf, ps){
                return(list(
                    read_length = 101,
                    frag_length = 180,
                    sample_plot = ggplot()
                ))
            }
        }
        sc = auto_frag_len_FUN(bam_treat_file, peak_set)
        fragment_lengths = c(sc$read_length, sc$frag_length)
        fl2color = c("gray", "red")
        names(fl2color) = fragment_lengths
        plots$fragment_length_calculation =
            sc$sample_plot +
            labs(title = basename(bam_treat_file),
                 subtitle = paste("read length:", sc$read_length,
                                  "\nfragment length:", sc$frag_length))
    }else{
        fl2color = rep(color_default, length(fragment_lengths))
        names(fl2color) = fragment_lengths
        if(!is.null(color_overrides)){
            fl2color[names(color_overrides)] = color_overrides
        }
    }


    .Object@fragment_lengths = as.integer(fragment_lengths)
    .Object@fl2color = fl2color
    .Object@plots = plots

    if(is.null(output_prefix)){
        output_prefix = sub("\\.bam$", "", basename(bam_treat_file))
    }
    .Object@output_prefix = output_prefix
    validObject(.Object)
    .Object

})



#' constructor for PeakRefiner
#'
#' @param peak_set GRanges formatted peak set to calculate motifs on
#' @param bam_treat_file .bam file of aligned reads from ChIP-seq pulldown.
#'   Index must be at .bam.bai
#' @param bam_input_file .bam file of aligned reads from input control. Index
#'   must be at .bam.bai
#' @param pwm Position Weight Matrix from PWMEnrich
#' @param target_pwm_names names of PWMs to include in figures by default
#' @param fragment_lengths fragment lengths to consider
#' @param color_overrides character. colors to use for each item in
#'   fragment_lengths.  should contain valid hex ("#000000") or R colors
#'   ("black") and be named with items in fragment_lengths.  non-overriden items
#'   will be black.
#' @param color_default single character. black.
#' @param auto_frag_len_FUN function to use to auto calculate fragment length.
#'   must accept two argument, bam_file and peak_set.
#' @param output_prefix prefix to use for output files.
#'
#' @return a new valid object of class PeakRefiner
#' @export
#' @importFrom methods new
#' @examples
#' bam_file = system.file("extdata", "MCF10A_CTCF.random5.bam", package = "peakrefine")
#' bam_input = system.file("extdata", "MCF10A_input.random5.bam", package = "peakrefine")
#' np = system.file("extdata", "MCF10A_CTCF.random5.narrowPeak", package = "peakrefine")
#' qgr = rtracklayer::import(np, format = "narrowPeak")
#'
#' library(PWMEnrich.Hsapiens.background)
#' data("PWMLogn.hg19.MotifDb.Hsap")
#' pwm = PWMLogn.hg19.MotifDb.Hsap[1:2]
#' PeakRefiner(qgr, bam_file, bam_input, pwm, names(pwm$pwms)[1])
PeakRefiner = function(peak_set,
                       bam_treat_file,
                       bam_input_file,
                       pwm,
                       target_pwm_names,
                       fragment_lengths = NULL,
                       color_overrides = NULL,
                       color_default = "black",
                       auto_frag_len_FUN = NULL,
                       output_prefix = NULL){
    new("PeakRefiner",
        peak_set = peak_set,
        bam_treat_file = bam_treat_file,
        bam_input_file = bam_input_file,
        pwm = pwm,
        target_pwm_names = target_pwm_names,
        fragment_lengths = fragment_lengths,
        color_overrides = color_overrides,
        color_default = color_default,
        auto_frag_len_FUN = auto_frag_len_FUN,
        output_prefix = output_prefix)
}

setMethod("show", "PeakRefiner",
          function(object) {
              message(paste(length(object@pwm$pwms), "motif PWMs"))
              message(paste(length(object@peak_set), "peaks/regions"))
              message("ChIP-seq bam file: ", object@bam_treat_file)
              message("input bam file: ", object@bam_input_file)
              message("fragment lengths are: ", paste(object@fragment_lengths, collapse = ", "))

              # nr_mat = nrow(object@hic_2d)
              # nr_reg = nrow(object@hic_1d)
              # covered = nr_mat / ((nr_reg^2 - nr_reg) / 2)
              # print(paste("size is", format(object.size(object), units = "GB")))
              # print(paste0(round(covered*100, 2), "% of bins have signal"))
          }
)

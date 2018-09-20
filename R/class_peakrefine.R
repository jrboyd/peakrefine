#S4 class to load and refine peak sets evaluated by motif enrichment
# validate config
# use cached data as possible

.datatable.aware = TRUE

#' Title
#'
#' @slot peak_set GRanges.
#' @slot bam_treat_file character.
#' @slot bam_input_file character.
#' @slot fragment_lengths integer.
#' @slot fl2color character.
#' @slot pwm PWMLognBackground.
#' @slot output_prefix character.
#' @slot plots list.
#'
#' @return
#' @export
#' @importClassesFrom  PWMEnrich PWMLognBackground
#' @importClassesFrom  GenomicRanges GRanges
#' @import PWMEnrich
#'
#' @examples
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

#' Title
#'
#' @param PeakRefiner
#'
#' @return
#' @export
#'
#' @importFrom methods validObject
#'
#' @examples
setMethod("initialize", "PeakRefiner", function(.Object,
                                                peak_set,
                                                bam_treat_file,
                                                bam_input_file,
                                                pwm,
                                                target_pwm_names,
                                                fragment_lengths = NULL,
                                                color_overrides = NULL,
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
            auto_frag_len_FUN = crossCorrByExtension()
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
        fl2color = rep("black", length(fragment_lengths))
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

#' HiC_matrix constructor
#'
#' @param matrix_file HiC-Pro matrix file
#' @param regions_file HiC-Pro region bed file matching matrix file
#' @param parameters HiC_parameters object
#'
#' @return
#' @export
#'
#' @importFrom methods new
#'
#' @examples
PeakRefiner = function(peak_set,
                       bam_treat_file,
                       bam_input_file,
                       pwm,
                       target_pwm_names,
                       fragment_lengths = NULL,
                       color_overrides = NULL,
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

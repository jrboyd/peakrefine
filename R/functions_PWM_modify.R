#' Calculate a large motif enrichment once for later rounds of
#' subsetting and analysis
#'
#' @param query_gr GRanges.  Regions to calculate cross correlation for.
#' @param pwm a position weight matrix object to pass to
#'   \code{\link{[PWMEnrich]{motifEnrichment}}}.
#' @param nbases number of bases around center of query_gr regions to consider
#'   for motif enrichment.
#' @param ncores number of computer cores to use.
#'
#' @return MotifEnrichmentResults for each region in \code{query_gr} and every
#'   motif in \code{pwm}
#' @export
#' @import PWMEnrich Biostrings
#' @examples
pre_motif = function(query_gr,
                     pwm,
                     seq_reference,
                     nbases = 200,
                     ncores = 32,
                     verbose = TRUE){
    # @importClassesFrom PWMEnrich MotifEnrichmentResults
    #TODO biocfilecache and digest args to auto cache, helper get cache file function
    options(mc.cores = ncores)
    PWMEnrich::registerCoresPWMEnrich(ncores)
    PWMEnrich::useBigMemoryPWMEnrich(TRUE)
    if(verbose) message("Fetching sequences...")
    sequences = Biostrings::getSeq(seq_reference,
                       resize(query_gr, nbases, fix = "center"))
    if(verbose) message("Calculating enrichments...")
    motif_res = PWMEnrich::motifEnrichment(sequences, pwm, verbose = verbose)
    if(verbose) message("Done! Please cache results.")
    return(motif_res)
}

#' a non-exported method from PWMEnrich, here for stability
PWMEnrich.affinitySequenceSet = function (scores, seq.len, pwm.len)
{
    if (is.vector(scores))
        scores = matrix(scores, nrow = 1, dimnames = list(NULL,
                                                          names(scores)))
    final = structure(rep(0, ncol(scores)), names = colnames(scores))
    for (i in 1:ncol(scores)) {
        seq.len.pwm = seq.len - pwm.len[i] + 1
        final[i] = sum(scores[, i] * seq.len.pwm)/sum(seq.len.pwm)
    }
    return(final)
}

#' a non-exported method from PWMEnrich, here for stability
PWMEnrich.cloverScore = function (scores, lr3 = FALSE, verbose = FALSE)
{
    if (is.vector(scores))
        return(scores)
    clover = matrix(0, nrow = nrow(scores), ncol = ncol(scores))
    colnames(clover) = colnames(scores)
    rownames(clover) = rownames(scores)
    for (m.inx in 1:ncol(scores)) {
        if (verbose)
            message(paste("Calculating Clover score for motif",
                          m.inx, "/", ncol(scores)))
        s = scores[, m.inx]
        N = length(s)
        A = matrix(NA, ncol = N + 1, nrow = N + 1)
        A[1, ] = 1
        for (i in 2:ncol(A)) {
            A[i, i - 1] = 0
        }
        for (i in 2:(N + 1)) {
            for (j in 2:(N + 1)) {
                if (j >= i) {
                    di = i - 1
                    dj = j - 1
                    A[i, j] = (di * s[dj] * A[i - 1, j - 1] + (dj -
                                                                   di) * A[i, j - 1])/dj
                }
            }
        }
        clover[, m.inx] = A[2:(N + 1), N + 1]
    }
    if (lr3)
        clover
    else colMeans(clover)
}


.recalc_MotifEnrichmentResults = function(res, pwmobj){
    sequences = res$sequences
    score = res$score
    pwms = res$pwms
    verbose = FALSE
    cutoff = NULL
    bg = res$bg
    group.only = FALSE

    ### from motifEnrichment
    if (score == "affinity") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        # res$sequence.nobg = motifScores(sequences, pwms, verbose = verbose)
        res$group.nobg = PWMEnrich.affinitySequenceSet(res$sequence.nobg,
                                                       seq.len, pwm.len)
    } else if (score == "clover") {
        # res$sequence.nobg = motifScores(sequences, pwms, verbose = verbose)
        res$group.nobg = PWMEnrich.cloverScore(res$sequence.nobg, verbose = verbose)
    } else if (score == "cutoff") {
        res$params = list(cutoff = cutoff)
        # res$sequence.nobg = motifScores(sequences, pwms, cutoff = cutoff,
        # verbose = verbose)
        res$group.nobg = colSums(res$sequence.nobg)
    } else {
        stop(paste("Unknown scoring algorithm: '", score, "'. Please select one of: 'affinity', 'cutoff', 'clover'",
                   sep = ""))
    }

    if (bg == "none") {
        # res$sequence.bg = NULL
        res$group.bg = NULL
    } else if (bg == "logn") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        # res$sequence.bg = PWMEnrich:::logNormPval(res$sequence.nobg, seq.len,
        #                               pwm.len, pwmobj@bg.mean, pwmobj@bg.sd, pwmobj@bg.len)
        # colnames(res$sequence.bg) = names(pwms)
        # res$sequence.norm = apply(res$sequence.bg, 1:2, qlnorm,
        #                           lower.tail = FALSE)
        if (score == "affinity") {
            if (is.matrix(pwmobj@bg.mean)) {
                res$group.bg = apply(res$sequence.bg, 2, function(x) {
                    pchisq(-2 * sum(log(x)), 2 * length(x), lower.tail = FALSE)
                })
                res$group.norm = sapply(res$group.bg, qlnorm,
                                        lower.tail = FALSE)
            } else {
                res$group.bg = PWMEnrich:::logNormPvalSequenceSet(res$sequence.nobg,
                                                      seq.len, pwm.len, pwmobj@bg.mean, pwmobj@bg.sd,
                                                      pwmobj@bg.len)
                res$group.norm = sapply(res$group.bg, qlnorm,
                                        lower.tail = FALSE)
            }
        } else if (score == "clover") {
            res$group.bg = PWMEnrich:::cloverScore(res$sequence.norm, verbose = verbose)
        }
    } else if (bg == "z") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        # res$sequence.bg = cutoffZscore(res$sequence.nobg, seq.len,
        #                                pwm.len, pwmobj@bg.P)
        res$group.bg = PWMEnrich:::cutoffZscoreSequenceSet(res$sequence.nobg,
                                               seq.len, pwm.len, pwmobj@bg.P)
    } else if (bg == "pval") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        usecutoff = NULL
        if (score == "cutoff")
            usecutoff = cutoff
        if (!group.only) {
            # res$sequence.bg = t(sapply(1:length(seq.len), function(i) empiricalPvalue(res$sequence.nobg[i,
            #                                                                                             ], seq.len[i], pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev,
            #                                                                           cutoff = usecutoff, B = B, verbose = verbose)))
        }
        if (score == "clover")
            res$group.bg = PWMEnrich:::cloverPvalue1seq(res$sequence.nobg,
                                            seq.len, pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev,
                                            B = B, verbose = verbose, clover = res$group.nobg)
        else res$group.bg = PWMEnrich:::empiricalPvalueSequenceSet(res$sequence.nobg,
                                                       seq.len, pwm.len, pwmobj@bg.fwd, pwmobj@bg.rev, cutoff = usecutoff,
                                                       B = B, verbose = verbose)
    } else if (bg == "ms") {
        usecutoff = NULL
        if (score == "cutoff")
            usecutoff = cutoff
        # res$sequence.bg = matrixShuffleZscorePerSequence(res$sequence.nobg,
        #                                                  sequences, pwms, cutoff = usecutoff, B = motif.shuffles)
    } else if (bg == "gev") {
        seq.len = sapply(sequences, length)
        pwm.len = sapply(pwms, length)
        # res$sequence.bg = gevPerSequence(res$sequence.nobg, seq.len,
        #                                  pwm.len, pwmobj@bg.loc, pwmobj@bg.scale, pwmobj@bg.shape)
        res$group.bg = NULL
    } else {
        stop(paste("Uknown background correction algorithm: '",
                   bg, "', Please select one of: 'none', 'logn', 'z', 'pval', 'ms'",
                   sep = ""))
    }

    if (group.only) {
        seq.n = grep("^sequence[.]", names(res))
        if (length(seq.n) > 0) {
            res = res[-seq.n]
        }
    }
    if ("sequence.nobg" %in% names(res)) {
        rownames(res$sequence.nobg) = names(sequences)
    }
    if ("sequence.bg" %in% names(res)) {
        rownames(res$sequence.bg) = names(sequences)
    }
    return(res)
}

#' subset MotifEnrichmentResults by region and update statistics
#'
#' @param me_res MotifEnrichmentResults
#' @param k logical or integer (as from which(logical)) to subset results by region.
#'
#' @return MotifEnrichmentResults subsetted by k
#' @export
#'
#' @importFrom stats pchisq qlnorm
#'
#' @examples
subset_MotifEnrichmentResults = function(me_res, k, pwmobj){
    ### top enriched can be cached but score and pval are group calcs
    res = me_res@res
    res$sequences = res$sequences[k]
    res$sequence.nobg = res$sequence.nobg[k, ]
    res$sequence.bg = res$sequence.bg[k, ]
    res$sequence.norm = res$sequence.norm[k, ]

    me_res@res = .recalc_MotifEnrichmentResults(res, pwmobj)
    return(me_res)
}

#' combine two MotifEnrichmentResults dervied from unique region sets and update statistics
#'
#' @param me_res1 MotifEnrichmentResults
#' @param me_res2 MotifEnrichmentResults
#'
#' @return MotifEnrichmentResults of combined inputs
#' @export
#'
#' @importFrom stats pchisq qlnorm
#'
#' @examples
add_MotifEnrichmentResults = function(me_res1, me_res2, pwmobj){
    ### top enriched can be cached but score and pval are group calcs
    res = me_res1@res
    res2 = me_res2@res
    res$sequences = c(res$sequences, res2$sequences)
    res$sequence.nobg = rbind(res$sequence.nobg, res2$sequence.nobg)
    res$sequence.bg = rbind(res$sequence.bg, res2$sequence.bg)
    res$sequence.norm = rbind(res$sequence.norm, res2$sequence.norm)

    me_res1@res = .recalc_MotifEnrichmentResults(res, pwmobj)
    return(me_res1)
}

qgr = seqsetvis::easyLoad_narrowPeak("/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled_peaks_passIDR.05.narrowPeak")[[1]]
qgr = dropSeqlevels(qgr, "chrU13369.1")
qgr = sample(qgr, 50)
bam_file = "/slipstream/galaxy/uploads/working/qc_framework/output_AF_MCF10_CTCF/MCF10A_CTCF_pooled/MCF10A_CTCF_pooled.bam"

###ChIPQC param setup
library(Rsamtools)
bamFile = bam_file
Window=400
FragmentLength=50
shiftWindowStart=50
shiftWindowEnd=300
verboseT=TRUE
runCrossCor = TRUE

### ChIPQC corr method
ChrLengths <- scanBamHeader(bamFile)[[1]]$targets
#
# if(length(ChrLengths[ChrLengths < shiftWindowEnd - shiftWindowStart]) > 0){
#     message("Removing ",length(ChrLengths[ChrLengths < shiftWindowEnd - shiftWindowStart]),
#             " chromosomes with length less than cross-coverage shift")
#     ChrLengths <- ChrLengths[!ChrLengths < shiftWindowEnd - shiftWindowStart]
# }
# if(verboseT == T){
#     message("Bam file has ",length(names(ChrLengths))," contigs")
# }
# if(!all(ChrOfInterest %in% names(ChrLengths))){
#     stop("Contigs of interest are not all in Bam file!!")
# }
# if(!is.null(ChrOfInterest)){
#     ChrLengths <- ChrLengths[names(ChrLengths) %in% ChrOfInterest]
# }

ShiftMatCor <- NULL
# for(k in 1:length(ChrLengths)){
    Param <- ScanBamParam(which=GRanges(seqnames=names(ChrLengths)[k],IRanges(start=1,end=unname(ChrLengths[names(ChrLengths) == names(ChrLengths)[k]])-shiftWindowEnd)),
                          what=c("flag","mapq"))
    temp <- GenomicAlignments::readGAlignments(bamFile,param=Param)
    if(k == 1){
        tocheckforreads <- 1000
        readlength=round(mean(width(temp[1:tocheckforreads])))
    }
    if(runCrossCor==TRUE){
        Sample_GIT <- GNCList(GRanges(seqnames=seqnames(temp),ranges=ranges(temp),strand=strand(temp),elementMetadata(temp)))
        PosCoverage <- coverage(IRanges(start(Sample_GIT[strand(Sample_GIT)=="+"]),start(Sample_GIT[strand(Sample_GIT)=="+"])),width=ChrLengths[k])
        NegCoverage <- coverage(IRanges(end(Sample_GIT[strand(Sample_GIT)=="-"]),end(Sample_GIT[strand(Sample_GIT)=="-"])),width=ChrLengths[k])

        ShiftsCorTemp <- shiftApply(seq(shiftWindowStart,shiftWindowEnd),
                                    PosCoverage,NegCoverage,cor, verbose = verboseT)
        ShiftMatCor <- cbind(ShiftMatCor,ShiftsCorTemp)
        colnames(ShiftMatCor)[ncol(ShiftMatCor)] <- names(ChrLengths)[k]
    }
# }
Weights <- ChrLengths
ShiftsCorAv <- apply(ShiftMatCor,1,function(x)weighted.mean(x,Weights[colnames(ShiftMatCor)],na.rm=TRUE))
plot(ShiftsCorAv, type = "l")

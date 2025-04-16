#' @title Estimate the allelic fraction of the pruned SNVs for a specific
#' profile
#'
#' @description The function estimates the allelic fraction of the
#' SNVs for a specific profile and add the information to the associated
#' Profile GDS file. The allelic fraction estimation method is adapted to
#' the type of study (DNA or RNA).
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Reference GDS file.
#'
#' @param gdsProfile an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Profile GDS file.
#'
#' @param currentProfile a \code{character} string corresponding to
#' the sample identifier as used in \code{\link{pruningSample}} function.
#'
#' @param studyID a \code{character} string corresponding to the name of
#' the study as
#' used in \code{\link{pruningSample}} function.
#'
#' @param chrInfo a \code{vector} of \code{integer} values representing
#' the length of the chromosomes. See 'details' section.
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "DNA" and "RNA". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"DNA"}.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#'
#' @param minProb a single \code{numeric} between 0 and 1 representing the
#' probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#'
#' @param eProb a single \code{numeric} between 0 and 1 representing the
#' probability of sequencing error. Default: \code{0.001}.
#'
#' @param cutOffLOH a single \code{numeric} representing the cutoff, in log,
#' for the homozygote score to assign a region as LOH.
#' Default: \code{-5}.
#'
#' @param cutOffHomoScore a single \code{numeric} representing the cutoff, in
#' log, that the SNVs in a block are called homozygote by error.
#' Default: \code{-3}.
#'
#' @param wAR a single positive \code{integer} representing the size-1 of
#' the window used to compute an empty box. Default: \code{9}.
#'
#' @param cutOffAR a single \code{numeric} representing the cutoff, in
#' log score, that the SNVs in a gene are allelic fraction different 0.5
#' Default: \code{3}.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Reference SNV Annotation GDS file.
#' This parameter is RNA specific.
#' Default: \code{NULL}.
#'
#' @param blockID a \code{character} string corresponding to the block
#' identifier in \code{gdsRefAnnot}. **This parameter is RNA specific.**
#' Default: \code{NULL}
#'
#' @param verbose a \code{logicial} indicating if the function should print
#' message when running. Default: \code{FALSE}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @details
#'
#' The `chrInfo` parameter contains the length of the chromosomes. The
#' length of the chromosomes can be obtain through the
#' \code{\link[GenomeInfoDb]{seqlengths}}
#' library.
#'
#' As example, for hg38 genome:
#'
#' ```
#'
#' if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'      chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#' }
#'
#' ```
#'
#' @examples
#'
#' ## Required library for GDS
#' library(gdsfmt)
#'
#' ## Path to the demo 1KG GDS file located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG.gds")
#'
#' ## Profile GDS file for one profile
#' fileProfile <- file.path(tempdir(), "ex1.gds")
#'
#' ## Copy the Profile GDS file demo that has been pruned and annotated
#' ## into current directory
#' file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                             fileProfile)
#'
#' ## Open the reference GDS file (demo version)
#' gds1KG <- snpgdsOpen(fileGDS)
#'
#' ## Profile GDS file for one profile
#' profileGDS <- openfn.gds(fileProfile, readonly=FALSE)
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#'
#'     ## Estimate the allelic fraction of the pruned SNVs
#'     estimateAllelicFraction(gdsReference=gds1KG, gdsProfile=profileGDS,
#'         currentProfile="ex1", studyID="MYDATA", chrInfo=chrInfo,
#'         studyType="DNA", minCov=10L, minProb=0.999, eProb=0.001,
#'         cutOffLOH=-5, cutOffHomoScore=-3, wAR=9, cutOffAR=3,
#'         gdsRefAnnot=NULL, blockID=NULL)
#'
#'     ## The allelic fraction is saved in the 'lap' node of Profile GDS file
#'     ## The 'lap' entry should be present
#'     profileGDS
#'
#'     ## Close both GDS files (important)
#'     closefn.gds(profileGDS)
#'     closefn.gds(gds1KG)
#'
#'     ## Remove Profile GDS file (created for demo purpose)
#'     unlink(fileProfile, force=TRUE)
#'
#' }
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
estimateAllelicFraction <- function(gdsReference, gdsProfile,
    currentProfile, studyID, chrInfo, studyType=c("DNA", "RNA"), minCov=10L,
    minProb=0.999, eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=9,
    cutOffAR=3, gdsRefAnnot=NULL, blockID=NULL, verbose=FALSE) {

    validateEstimateAllelicFraction(gdsReference=gdsReference,
        gdsProfile=gdsProfile, currentProfile=currentProfile, studyID=studyID,
        chrInfo=chrInfo, studyType=studyType, minCov=minCov, minProb=minProb,
        eProb=eProb, cutOffLOH=cutOffLOH, cutOffHomoScore=cutOffHomoScore,
        wAR=wAR, cutOffAR=cutOffAR, gdsRefAnnot=gdsRefAnnot, blockID=blockID,
        verbose=verbose)

    ## Set study type
    studyType <- arg_match(studyType)

    snpPos <- NULL

    ## The type of study affects the allelic fraction estimation
    if(studyType == "DNA") {
        snpPos <- computeAllelicFractionDNA(gdsReference=gdsReference,
            gdsSample=gdsProfile, currentProfile=currentProfile,
            studyID=studyID, chrInfo=chrInfo, minCov=minCov, minProb=minProb,
            eProb=eProb, cutOffLOH=cutOffLOH, cutOffHomoScore=cutOffHomoScore,
            wAR=wAR, verbose=verbose)
    } else if(studyType == "RNA") {
        snpPos <- computeAllelicFractionRNA(gdsReference=gdsReference,
            gdsSample=gdsProfile, gdsRefAnnot=gdsRefAnnot,
            currentProfile=currentProfile, studyID=studyID, blockID=blockID,
            chrInfo=chrInfo, minCov=minCov, minProb=minProb, eProb=eProb,
            cutOffLOH=cutOffLOH, cutOffAR=cutOffAR, verbose=verbose)
    }



    ## Calculate the cumulative sum for each chromosome
    cumSumResult <- lapply(unique(snpPos$snp.chr), function(i) {
        snpChr <- snpPos[snpPos$snp.chr == i, ]
        tmp <- c(0, abs(snpChr[2:nrow(snpChr), "lap"] -
                    snpChr[seq_len(nrow(snpChr)- 1),  "lap"]) > 1e-3)
        return(cumsum(tmp))
    })

    # Find segment with same lap
    snpPos$seg <- rep(0, nrow(snpPos))
    k <- 1
    for(i in seq_len(length(unique(snpPos$snp.chr)))) {
        ##snpChr <- snpPos[snpPos$snp.chr == chr, ]
        ##tmp <- c(0, abs(snpChr[2:nrow(snpChr), "lap"] -
        ##                snpChr[seq_len(nrow(snpChr)- 1),  "lap"]) > 1e-3)
        chr <- unique(snpPos$snp.chr)[i]
        snpPos$seg[snpPos$snp.chr == chr] <- cumSumResult[[i]] + k
        k <- max(snpPos$seg[snpPos$snp.chr == chr]) + 1
    }

    ## Save information into the "lap" node in the Profile GDS file
    ## Save information into the "segment" node in the Profile GDS file
    ## Suppose we keep only the pruned SNVs
    addUpdateLap(gdsProfile, snpPos$lap[which(snpPos$pruned == TRUE)])
    addUpdateSegment(gdsProfile, snpPos$seg[which(snpPos$pruned == TRUE)])

    return(0L)
}


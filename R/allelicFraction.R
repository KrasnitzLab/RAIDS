#' @title TODO
#'
#' @description TODO
#'
#' @param snp.pos a \code{data.frame} containing the genotype information for
#' a SNV dataset. TODO
#'
#' @param w a single positive \code{numeric} representing the size of the
#' window to compute the allelic fraction.
#' Default: \code{10}.
#'
#' @param cutOff a single \code{numeric} representing TODO. Default: \code{-3}.
#'
#' @return a \code{matrix} of \code{numeric} with 3 columns where each
#' row represent a segment
#' of imbalanced SNVs. The first column represents the position, in
#' \code{snp.pos}, of the first
#' SNV in the segment. The second column represents the position, in the
#' \code{snp.pos}, of the last SNV in the segment. The third column represents
#' the lower allelic frequency of the segment and is \code{NA} when the value
#' cannot be calculated. The value \code{NULL} is
#' returned when none of the SNVs
#' tested positive for the imbalance.
#'
#' @examples
#'
#' ## Data frame with SNV information for the specified chromosome (chr 1)
#' snpInfo <- data.frame(cnt.tot=c(41, 17, 27, 15, 11, 37, 16, 32),
#'     cnt.ref=c(40, 17, 27, 15, 4, 14, 16, 32),
#'     cnt.alt=c(0, 0, 0, 0, 7, 23, 0, 0),
#'     snp.pos=c(3722256, 3722328, 3767522, 3868160, 3869467, 4712655,
#'         6085318, 6213145),
#'     snp.chr=c(rep(1, 8)),
#'     normal.geno=c(rep(1, 8)),
#'     pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, rep(TRUE, 3)),
#'     snp.index=c(160, 162, 204, 256, 259, 288, 366, 465),
#'     keep=rep(TRUE, 8), hetero=c(rep(FALSE, 4), TRUE, TRUE, rep(FALSE, 2)),
#'     homo=c(rep(TRUE, 4), FALSE, FALSE, TRUE, TRUE),
#'     lap=rep(-1, 8), LOH=rep(0, 8), imbAR=rep(-1, 8),
#'     stringAsFactor=FALSE)
#'
#' ## The function returns NULL when there is not imbalanced SNVs
#' RAIDS:::computeAlleleFraction(snp.pos=snpInfo, w=10, cutOff=-3)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom stats median
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeAlleleFraction <- function(snp.pos, w=10, cutOff=-3) {

    listBlockAR <- list()
    j <- 1
    tmp <- as.integer(snp.pos$imbAR == 1)
    z <- cbind(c(tmp[1], tmp[-1] - tmp[seq_len(length(tmp) -1)]),
            c(tmp[-1] - tmp[seq_len(length(tmp) -1)], tmp[length(tmp)] * -1))

    ## Split SNVs by segment of continuous imbalanced SNVs
    ## There must be at least one segment with imbalanced SNVs to go one
    if(length(which(z[,1] == 1)) > 0) {
        ## Find segmentsof imbalanced SNVs
        segImb <- data.frame(start=seq_len(nrow(snp.pos))[which(z[,1] > 0)],
                                end=seq_len(nrow(snp.pos))[which(z[,2] < 0)])

        for(i in seq_len(nrow(segImb))) {
            # index of the segment
            listSeg <- (segImb$start[i]):(segImb$end[i])
            # index hetero segment
            listHetero  <- listSeg[snp.pos[listSeg,"hetero"] == TRUE]
            # SNP hetero for the segment
            snp.hetero <- snp.pos[listHetero,]

            if(nrow(snp.hetero) >= 2 * w) {
                lapCur <- median(apply(snp.hetero[seq_len(w),
                            c("cnt.ref", "cnt.alt")], 1, min) /
                    (rowSums(snp.hetero[seq_len(w),c("cnt.ref", "cnt.alt")])))

                start <- 1
                k <- w + 1
                while(k < nrow(snp.hetero)) {
                    # We have (k+w-1) <= nrow(snp.hetero)
                    # Case 1 true because (nrow(snp.hetero) >= 2 * w
                    # Other case nrow(snp.hetero) >= w+k - 1
                    curWin <- testAlleleFractionChange(snp.hetero[k:(k+w-1),
                                    c("cnt.ref", "cnt.alt")], cutOff, lapCur)

                    if(curWin$pCut1 == 1){ # new Region the allelicFraction
                        # table of the index of the block with lapCur
                        listBlockAR[[j]] <- c(listHetero[start],
                                                listHetero[k], lapCur)

                        lapCur <- median(apply(snp.hetero[k:(k+w-1),
                                    c("cnt.ref", "cnt.alt")], 1, min) /
                                    (rowSums(snp.hetero[k:(k+w-1),
                                                c("cnt.ref", "cnt.alt")])))

                        start <- k

                        if(nrow(snp.hetero) - start < w) { # Close the segment
                            lapCur <-
                                median(apply(snp.hetero[start:nrow(snp.hetero),
                                    c("cnt.ref", "cnt.alt")], 1, min) /
                                (rowSums(snp.hetero[start:nrow(snp.hetero),
                                    c("cnt.ref", "cnt.alt")])))

                            listBlockAR[[j]] <- c(listHetero[start],
                                                    segImb$end[i], lapCur)

                            j <- j+1
                            k <- nrow(snp.hetero)
                        }else{ # nrow(snp.hetero) >= w+k
                            k<- k + 1
                            j <- j + 1

                        }
                    }else{ # keep the same region
                        if((nrow(snp.hetero) - k ) < w){ # close
                            lapCur <-
                                median(apply(snp.hetero[start:nrow(snp.hetero),
                                    c("cnt.ref", "cnt.alt")], 1, min) /
                                (rowSums(snp.hetero[start:nrow(snp.hetero),
                                    c("cnt.ref", "cnt.alt")])))

                            listBlockAR[[j]] <- c(listHetero[start],
                                                    segImb$end[i], lapCur)

                            j <- j + 1

                            k <- nrow(snp.hetero)
                        } else{ # continue nrow(snp.hetero) >= w+k
                            lapCur <- median(apply(snp.hetero[start:k,
                                c("cnt.ref", "cnt.alt")], 1, min) /
                                    (rowSums(snp.hetero[start:k,c("cnt.ref",
                                                        "cnt.alt")])))

                            k <- k + 1
                        }
                    }
                }# End while
            }else {
                lapCur <- median(apply(snp.hetero[, c("cnt.ref", "cnt.alt")],
                                1, min) / (rowSums(snp.hetero[,c("cnt.ref",
                                                        "cnt.alt")])))

                listBlockAR[[j]] <- c(segImb$start[i], segImb$end[i], lapCur)

                j <- j + 1
            }
        }
    }

    # note NULL if length(listBlockAR) == 0
    listBlockAR <- do.call(rbind, listBlockAR)

    return(listBlockAR)
}


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
#' @param minProb a single \code{numeric} between 0 and 1 representing TODO.
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
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
#'
#' ## Profile GDS file for one profile
#' fileProfile <- file.path("ex1.gds")
#'
#' ## Example can only be run if the current directory is in writing mode
#' if (file.access(getwd()) == 0 && !file.exists(fileProfile))  {
#'     ## Copy the Profile GDS file demo that has been pruned and annotated
#'     ## into  current directory
#'     file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                         "ex1.gds")
#'
#'     ## Open the reference GDS file (demo version)
#'     gds1KG <- snpgdsOpen(fileGDS)
#'
#'     ## Profile GDS file for one profile
#'     profileGDS <- openfn.gds(fileProfile, readonly=FALSE)
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'         181538259L, 170805979L, 159345973L, 145138636L, 138394717L,
#'         133797422L, 135086622L, 133275309L, 114364328L, 107043718L,
#'         101991189L, 90338345L, 83257441L,  80373285L,  58617616L,
#'         64444167L,  46709983L, 50818468L, 156040895L, 57227415L,  16569L)
#'
#'     ## A formal way to get the chormosome length information
#'     ## library(GenomeInfoDb)
#'     ## library(BSgenome.Hsapiens.UCSC.hg38)
#'     ## chrInfo <- GenomeInfoDb::seqlengths(Hsapiens)[1:25]
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
#'     ## Unlink Profile GDS file (created for demo purpose)
#'     unlink("ex1.gds", force=TRUE)
#' }
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

    snp.pos <- NULL

    ## The type of study affects the allelic fraction estimation
    if(studyType == "DNA") {
        snp.pos <- computeAllelicFractionDNA(gdsReference=gdsReference,
            gdsSample=gdsProfile, currentProfile=currentProfile,
            studyID=studyID, chrInfo=chrInfo, minCov=minCov, minProb=minProb,
            eProb=eProb, cutOffLOH=cutOffLOH, cutOffHomoScore=cutOffHomoScore,
            wAR=wAR, verbose=verbose)
    } else if(studyType == "RNA") {
        snp.pos <- computeAllelicFractionRNA(gdsReference=gdsReference,
            gdsSample=gdsProfile, gdsRefAnnot=gdsRefAnnot,
            currentProfile=currentProfile, studyID=studyID, blockID=blockID,
            chrInfo=chrInfo, minCov=minCov, minProb=minProb, eProb=eProb,
            cutOffLOH=cutOffLOH, cutOffAR=cutOffAR, verbose=verbose)
    }

    snp.pos$seg <- rep(0, nrow(snp.pos))
    k <- 1
    # Find segment with same lap
    for(chr in seq_len(22)) {
        snpChr <- snp.pos[snp.pos$snp.chr == chr, ]
        tmp <- c(0, abs(snpChr[2:nrow(snpChr), "lap"] -
                        snpChr[seq_len(nrow(snpChr)- 1),  "lap"]) > 1e-3)
        snp.pos$seg[snp.pos$snp.chr == chr] <- cumsum(tmp) + k
        k <- max(snp.pos$seg[snp.pos$snp.chr == chr]) + 1
    }

    ## Save information into the "lap" node in the Profile GDS file
    ## Save information into the "segment" node in the Profile GDS file
    ## Suppose we keep only the pruned SNVs
    addUpdateLap(gdsProfile, snp.pos$lap[which(snp.pos$pruned == TRUE)])
    addUpdateSegment(gdsProfile, snp.pos$seg[which(snp.pos$pruned == TRUE)])

    return(0L)
}


#' @title Extract the genotype information for a SNV dataset using
#' the Profile GDS file and the Reference GDS file
#'
#' @description The function generates a \code{data.frame} containing the
#' genotype information from a initial list of SNVs associated to a specific
#' profile. The function uses the information present in the Reference GDS file
#' and the Profile GDS file.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}} (a
#' GDS file), the opened Reference GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Profile GDS file.
#'
#' @param currentProfile a \code{character} string corresponding to
#' the sample identifier used in \code{\link{pruningSample}} function.
#'
#' @param studyID a \code{character} string corresponding to the study
#' identifier used in \code{\link{pruningSample}} function.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#'
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#'
#' @param eProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error. Default: \code{0.001}.
#'
#' @param verbose a \code{logicial} indicating if messages should be printed
#' when the function is running.
#'
#' @return a \code{data.frame} containing:
#' \itemize{
#' \item{cnt.tot} {a single \code{integer} representing the total coverage for
#' the SNV.}
#' \item{cnt.ref} {a single \code{integer} representing the coverage for
#' the reference allele.}
#' \item{cnt.alt} {a single \code{integer} representing the coverage for
#' the alternative allele.}
#' \item{snpPos} {a single \code{integer} representing the SNV position.}
#' \item{snp.chr} {a single \code{integer} representing the SNV chromosome.}
#' \item{normal.geno} {a single \code{numeric} indicating the genotype of the
#' SNV. The possibles are: \code{0} (wild-type homozygote), \code{1}
#' (heterozygote), \code{2} (altenative homozygote), \code{3} indicating that
#' the normal genotype is unknown.}
#' \item{pruned} { a \code{logical}}
#' \item{snp.index} {a \code{vector} of \code{integer} representing the
#' position of the SNVs in the Reference GDS file.}
#' \item{keep} {a \code{logical} }
#' \item{hetero} {a \code{logical} }
#' \item{homo} {a \code{logical} }
#' }
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
#'
#' ## Temporary Profile GDS file for one profile
#' fileProfile <- file.path(getwd(), "ex1.gds")
#'
#' ## Example can only be run if the current directory is in writing mode
#' if (file.access(getwd()) == 0 && !file.exists(fileProfile)) {
#'
#'     ## Copy the Profile GDS file demo that has been pruned and annotated
#'     file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                  fileProfile)
#'
#'     ## Open the reference GDS file (demo version)
#'     gds1KG <- snpgdsOpen(fileGDS)
#'
#'     ## Open Profile GDS file for one profile
#'     profileGDS <- openfn.gds(fileProfile)
#'
#'     ## The function returns a data frame containing the SNVs information
#'     result <- RAIDS:::getTableSNV(gdsReference=gds1KG, gdsSample=profileGDS,
#'         currentProfile="ex1", studyID="MYDATA", minCov=10L, minProb=0.999,
#'         eProb=0.001, verbose=FALSE)
#'     head(result)
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
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
getTableSNV <- function(gdsReference, gdsSample, currentProfile, studyID,
                minCov=10, minProb=0.999, eProb=0.001, verbose) {

    ## Extract study information (data.frame) from GDS Sample file
    study.annot <- read.gdsn(index.gdsn(node=gdsSample, path="study.annot"))

    ## Retain the specified sample in the specified study
    posCur <- which(study.annot$data.id == currentProfile &
                        study.annot$study.id == studyID)

    ## Extract SNV coverage from Profile GDS file
    cnt.total <- read.gdsn(node=index.gdsn(gdsSample, "Total.count"),
                            start=c(1, posCur), count=c(-1, 1))

    ## Only retained the SNV with the minimum required coverage
    listKeep <- cnt.total@i[which(cnt.total@x >= minCov)] + 1

    ## Create the data.frame with the required information
    snpPos <- data.frame(cnt.tot=cnt.total[listKeep],
                    cnt.ref=read.gdsn(index.gdsn(gdsSample, "Ref.count"),
                        start=c(1, posCur), count=c(-1, 1))[listKeep],
                    cnt.alt=read.gdsn(index.gdsn(gdsSample, "Alt.count"),
                        start=c(1, posCur), count=c(-1, 1))[listKeep],
                    snp.pos=read.gdsn(index.gdsn(node=gdsReference,
                                    "snp.position"))[listKeep],
                    snp.chr=read.gdsn(index.gdsn(node=gdsReference,
                                    "snp.chromosome"))[listKeep],
                    normal.geno=rep(3, length(listKeep)),#Norm genotype unknown
                    pruned=rep(FALSE, length(listKeep)), #bit(length(listKeep))
                    snp.index=listKeep, stringsAsFactors=FALSE)

    snp.pruned <- read.gdsn(index.gdsn(node=gdsSample, "snp.index"))

    listKeepPruned <- which(listKeep %in% snp.pruned)
    snpPos$pruned[listKeepPruned] <- TRUE

    rm(cnt.total, snp.pruned, listKeepPruned)

    # Add snv info for snv in normal
    if ("normal.geno" %in% ls.gdsn(node=gdsSample)) {
        # if normal.geno exist mean there is count not in the ref
        # I have other genotype than 1KG
        if (verbose) {
            message("Genotype")
        }

        cnt.total <- read.gdsn(index.gdsn(gdsSample, "Total.count.o"))
        listKeep.o <- which(cnt.total >= minCov)

        snpPosO <- data.frame(cnt.tot=cnt.total[listKeep.o],
            cnt.ref=read.gdsn(index.gdsn(gdsSample, "Ref.count.o"))[listKeep.o],
            cnt.alt=read.gdsn(index.gdsn(gdsSample, "Alt.count.o"))[listKeep.o],
            snp.pos=read.gdsn(index.gdsn(gdsReference,
                "snp.position.o"))[listKeep.o],
            snp.chr=read.gdsn(index.gdsn(gdsReference,
                "snp.chromosome.o"))[listKeep.o],
            normal.geno=read.gdsn(index.gdsn(gdsSample,
                "normal.geno"))[listKeep.o], pruned=rep(0, length(listKeep)),
            snp.index=rep(0, length(listKeep.o)), stringsAsFactors=FALSE)
        listChr <- unique(snpPosO$snp.chr)
        listUnion <- list()

        # if snpPosO intersect snpPos and normal.geno != 3 (we know
        # the genotype of normal) change in snpPos normal.geno
        z <- cbind(c(snpPosO$snp.chr, snpPos$snp.chr, snpPosO$snp.chr),
                c(snpPosO$snp.pos, snpPos$snp.pos, snpPosO$snp.pos),
                c(seq_len(nrow(snpPosO)), 0, -1*seq_len(nrow(snpPosO))),
                c(rep(0, nrow(snpPosO)), seq_len(nrow(snpPos)),
                        rep(0, nrow(snpPosO))))
        z <- z[order(z[,1], z[,2], z[,3]), ]
        vCum <- cumsum(z[,3])

        snpPos[z[ vCum < 0 & z[,3] == 0, 4], "normal.geno"] <-
            snpPosO[vCum[vCum < 0 & z[, 3] == 0], "normal.geno"]
        rm(z)

        # Keep the snpPosO not in snpPos
        z <- cbind(c(snpPos$snp.chr, snpPosO$snp.chr, snpPos$snp.chr),
                c(snpPos$snp.pos, snpPosO$snp.pos, snpPos$snp.pos),
                c(seq_len(nrow(snpPos)), 0, -1*seq_len(nrow(snpPos))),
                c(rep(0, nrow(snpPos)), seq_len(nrow(snpPosO)),
                                rep(0, nrow(snpPos))))
        z <- z[order(z[,1], z[,2], z[,3]), ]
        # merge snpPos with snpPosO
        snpPos <- rbind(snpPos,
                            snpPosO[z[cumsum(z[,3] == 0 & z[,3] == 0),4],])

    }
    listCnt <- unique(snpPos$cnt.tot)
    listCnt <- listCnt[order(listCnt)]

    cutOffA <- data.frame(count=unlist(vapply(listCnt,
        FUN=function(x, minProb, eProb){
                return(max(2, qbinom(minProb, x,eProb)))},
                        FUN.VALUE = numeric(1), minProb=minProb, eProb=eProb)),
                            allele = unlist(vapply(listCnt,
                            FUN=function(x, minProb, eProb){
                                return(max(2,qbinom(minProb, x,eProb)))},
                        FUN.VALUE = numeric(1), minProb=minProb, eProb=eProb)))
    row.names(cutOffA) <- as.character(listCnt)

    snpPos$keep <- rowSums(snpPos[, c("cnt.ref", "cnt.alt")]) >=
        snpPos$cnt.tot - cutOffA[as.character(snpPos$cnt.tot), "count"]

    snpPos$hetero <- snpPos$keep == TRUE &
        rowSums(snpPos[, c("cnt.ref", "cnt.alt")] >=
            cutOffA[as.character(snpPos$cnt.tot), "allele"]) == 2

    # We set to homo if 2th allele can be explain by error
    # can switch low allelic fraction to LOH which is less a problem
    # then reduce the allelic ratio by seq error
    snpPos$homo <- snpPos$keep == TRUE &
        rowSums(snpPos[, c("cnt.ref", "cnt.alt")] >=
            cutOffA[as.character(snpPos$cnt.tot), "allele"]) == 1

    ## If we know the normal is hetero then we call hetero
    ## if the cnt.alt and cnt.ref > 0
    listHeteroN <- which(snpPos$homo == TRUE &
                rowSums(snpPos[, c("cnt.ref", "cnt.alt")] > 0) == 2 &
                snpPos$normal.geno == 1)

    if (length(listHeteroN) > 0) {
        snpPos$hetero[listHeteroN] <- TRUE
        snpPos$homo <- FALSE
    }

    return(snpPos)
}


#' @title Identify regions of LOH on one chromosome using homozygote SNVs
#'
#' @description The function identifies regions of LOH on a specific
#' chromosome using the homozygote SNVs present on the chromosome.
#'
#' @param gdsReference an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, an
#' opened Reference GDS file.
#'
#' @param chrInfo a \code{vector} of \code{integer} representing the length of
#' the chromosomes. As an example, the information ca be obtained from
#' package 'BSgenome.Hsapiens.UCSC.hg38'.
#'
#' @param snpPos a \code{data.frame} containing the SNV information for the
#' chromosome specified by the \code{chr} argument. The \code{data.frame} must
#' contain:
#' \itemize{
#' \item{cnt.tot} {a single \code{integer} representing the total coverage for
#' the SNV.}
#' \item{cnt.ref} {a single \code{integer} representing the coverage for
#' the reference allele.}
#' \item{cnt.alt} {a single \code{integer} representing the coverage for
#' the alternative allele.}
#' \item{snp.pos} {a single \code{integer} representing the SNV position.}
#' \item{snp.chr} {a single \code{integer} representing the SNV chromosome.}
#' \item{normal.geno} {a single \code{numeric} indicating the genotype of the
#' SNV. The possibles are: \code{0} (wild-type homozygote), \code{1}
#' (heterozygote), \code{2} (altenative homozygote), \code{3} indicating that
#' the normal genotype is unknown.}
#' \item{pruned} {a \code{logical} indicating if the SNV is retained after
#' pruning}
#' \item{snp.index} {a \code{integer} representing the index position of the
#' SNV in the Reference GDS file that contains all SNVs}
#' \item{keep} {a \code{logical} indicating if the genotype exists for the SNV}
#' \item{hetero} {a \code{logical} indicating if the SNV is heterozygote}
#' \item{homo} {a \code{logical} indicating if the SNV is homozygote}
#' }
#'
#' @param chr a single positive \code{integer} for the current chromosome. The
#' \code{chrInfo} parameter must contain the value for the specified
#' chromosome.
#'
#' @param  genoN a single \code{numeric} between 0 and 1 representing the
#' probability of sequencing error. Default: \code{0.0001}.
#'
#' @return a \code{data.frame} with the informations about LOH on a specific
#' chromosome. The \code{data.frame} contains those columns:
#' \itemize{
#' \item{chr} {a \code{integer} representing the current chromosome}
#' \item{start} {a \code{integer} representing the starting position on the
#' box containing only homozygote SNVs (or not SNV). The first box starts at
#' position 1.}
#' \item{end} {a \code{integer} representing the end position on the
#' box containing only homozygote SNVs (or not SNV). The last box ends at the
#' length of the chromosome.}
#' \item{logLHR} {TOREVIEW Score for LOH base on frequencies in population. Sum of
#' the log10 of the frequencies of the observe gegenotype minus the
#' the sum of the log10 of the higher frequent genotype.
#' (-100 when normal genotype are present)}
#' \item{LH1} {TOREVIEW If normal genotype present and heterozygot the probability to be heterozygote
#' base on the coverage of each allele}
#' \item{LM1} {TOREVIEW If normal genotype present and heterozygot the max probability the max probability
#' for the read coverage at the position}
#' \item{homoScore} {a \code{numeric} representing \code{LH1} - \code{LM1}}
#' \item{nbSNV} {a \code{integer} representing th number of SNVs in
#' the box}
#' \item{nbPruned} {a \code{integer} representing the number of pruned SNVs in
#' the box}
#' \item{nbNorm} {TOREVIEW a \code{integer} representing of genotype
#' heterozygote for the normal in the block}
#' \item{LOH} {TOREVIEW a \code{integer} representing a flag if 1 it mean
#' the block is satisfy the criteria to be LOH. The value is not assign
#' in this function they are all 0}
#' }
#'
#' @examples
#'
#' ## Required library for GDS
#' library(SNPRelate)
#'
#' ## Path to the demo Reference GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
#'
#' ## Open the Reference GDS file (demo version)
#' gds1KG <- snpgdsOpen(fileGDS)
#'
#' ## Chromosome length information for hg38
#' ## chr23 is chrX, chr24 is chrY and chrM is 25
#' chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'     181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
#'     135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
#'     83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
#'     156040895L, 57227415L,  16569L)
#'
#' ## Data frame with SNV information for the specified chromosome (chr 1)
#' snpInfo <- data.frame(cnt.tot=c(41, 17, 27, 15, 11, 37, 16, 32),
#'     cnt.ref=c(40, 17, 27, 15, 4, 14, 16, 32),
#'     cnt.alt=c(0, 0, 0, 0, 7, 23, 0, 0),
#'     snp.pos=c(3722256, 3722328, 3767522, 3868160, 3869467, 4712655,
#'         6085318, 6213145),
#'     snp.chr=c(rep(1, 8)),
#'     normal.geno=c(rep(3, 8)), pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE,
#'     TRUE, TRUE),
#'     pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, rep(TRUE, 3)),
#'     snp.index=c(160, 162, 204, 256, 259, 288, 366, 465),
#'     keep=rep(TRUE, 8), hetero=c(rep(FALSE, 4), TRUE, TRUE, rep(FALSE, 2)),
#'     homo=c(rep(TRUE, 4), FALSE, FALSE, TRUE, TRUE),
#'     stringAsFactor=FALSE)
#'
#' ## The function returns a data frame containing the information about the
#' ## LOH regions in the specified chromosome
#' result <- RAIDS:::computeLOHBlocksDNAChr(gdsReference=gds1KG,
#'     chrInfo=chrInfo, snpPos=snpInfo, chr=1L, genoN=0.0001)
#' head(result)
#'
#' ## Close Reference GDS file (important)
#' closefn.gds(gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats dbinom
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeLOHBlocksDNAChr <- function(gdsReference, chrInfo, snpPos, chr,
                                        genoN=0.0001) {

    genoN1 <- 1 - 2 * genoN

    chrEnd <- chrInfo[chr]
    listHetero <- snpPos[snpPos$hetero == TRUE, "snp.pos"]

    homoBlock <- data.frame(chr=rep(chr, length(listHetero) + 1),
        start=c(1, listHetero + 1), end=c(listHetero, chrEnd))

    z <- cbind(c(homoBlock$start, homoBlock$end,
            snpPos$snp.pos[which(snpPos$homo == TRUE)]),
            c(seq_len(length(homoBlock$start)),
                -1*seq_len(length(homoBlock$start)),
                rep(0, length(which(snpPos$homo == TRUE)))),
            c(rep(0, length(homoBlock$start)), rep(0, length(homoBlock$start)),
                seq_len(length(which(snpPos$homo == TRUE)))))

    z <- z[order(z[, 1]), ]

    blcSNV <- data.frame(block=cumsum(z[, 2])[z[, 2] == 0],
                            snv=z[z[, 2] == 0, 3])

    listAF <- read.gdsn(index.gdsn(gdsReference, "snp.AF"))

    # Compute if the block is LOH
    homoBlock$logLHR <- rep(0, nrow(homoBlock))
    homoBlock$LH1 <- rep(0, nrow(homoBlock))
    homoBlock$LM1 <- rep(0, nrow(homoBlock))
    # homoScore LH1 - LM1
    homoBlock$homoScore <- rep(0, nrow(homoBlock))
    homoBlock$nbSNV <- rep(0, nrow(homoBlock))
    homoBlock$nbPruned <- rep(0, nrow(homoBlock))
    homoBlock$nbNorm <- rep(0, nrow(homoBlock))
    # Include a field for LOH but will be fill elsewhere
    homoBlock$LOH <- rep(0, nrow(homoBlock))

    homoBlock <- lapply(seq_len(nrow(homoBlock)),
                    FUN=function(i, homoBlock, blcSNV,
                                 listAF,snpPos){
                        blcCur <- blcSNV[blcSNV$block == i, ]
                        snvH <- snpPos[blcCur$snv, ]
                        lH1 <- 0
                        lM1 <- 0
                        logLHR <- 0
                        homoBlock$nbSNV[i] <- nrow(blcCur)
                        homoBlock$nbPruned[i] <- length(which(snvH$pruned))
                        if (length(which(snvH$normal.geno != 3)) > 0) {
                            listCount <- snvH$cnt.tot[which(snvH$normal.geno == 1)]
                            homoBlock$nbNorm[i] <- length(listCount)
                            if(homoBlock$nbNorm[i] > 0){
                                lH1 <-sum(log10(apply(snvH[which(snvH$normal.geno == 1),
                                                           c("cnt.ref", "cnt.tot"), drop=FALSE],
                                                      1, FUN=function(x){
                                                          return(dbinom(x[1], x[2], 0.5))
                                                          # genoN1 * dbinom(x[1], x[2], 0.5) + genoN
                                                      })))

                                lM1 <- sum(log10(apply(snvH[which(snvH$normal.geno == 1),
                                                            c("cnt.ref", "cnt.tot"), drop=FALSE],
                                                       1, FUN=function(x){
                                                           return(dbinom((x[2] + x[2]%%2)/2, x[2], 0.5))
                                                           #genoN1 *dbinom((x[2] + x[2]%%2)/2, x[2], 0.5) + genoN
                                                       })))
                                logLHR <- -100
                            }

                        } else if (length(which(snvH$pruned)) > 2) {
                            afSNV <- listAF[snvH$snp.index[which(snvH$pruned)]]
                            afSNV <- apply(X=matrix(afSNV, ncol=1), MARGIN=1,
                                           FUN=function(x){max(x, 0.01)})
                            snvR <- snvH$cnt.ref[which(snvH$pruned)] >
                                snvH$cnt.alt[which(snvH$pruned)]

                            # Check if it is unlikely the genotype are homo by error
                            lH1 <- -100
                            # Freq of the more likely geno

                            tmp <- apply(matrix(afSNV, ncol=1), 1,
                                         FUN=function(x){max(max(x, 1-x)^2, 2* x *(1-x)) })
                            # log10 (prod(FreqAllele^2) / prod(freq of more likely genotype))
                            # snvR * 1 + (-1)^snvR * afSNV freq of the genotype
                            # (snvR = 1 homo ref
                            # and 0 if homo alt)
                            logLHR <- sum(2 * log10(snvR * 1 + (-1)^snvR * afSNV)) -
                                sum(log10(tmp))
                        }

                        homoBlock$logLHR[i] <- max(logLHR, -100)
                        homoBlock$LH1[i] <- lH1
                        homoBlock$LM1[i] <- lM1
                        homoBlock$homoScore[i] <- lH1 - lM1
                        return(homoBlock[i,])
                    }, homoBlock=homoBlock,
                    blcSNV=blcSNV, listAF=listAF,
                    snpPos=snpPos)
    homoBlock = do.call(rbind,homoBlock)

    return(homoBlock)
}

#' @title TOREVIEW For each imbalance segments compute the allelic fraction
#'
#' @description TOREVIEW Compute the allelic fraction for each segment different
#' than 0.5. The allelic fraction of the segment can be decomposed in
#' sub-segments.
#'
#' @param snpPos a \code{data.frame} containing the genotype information for
#' a SNV dataset.
#'
#' @param w a single positive \code{numeric} representing the size of the
#' window to compute the allelic fraction.
#' Default: \code{10}.
#'
#' @param cutOff a \code{numeric} representing the cut-off for considering
#' a region imbalanced when comparing likelihood to gave allelic fraction
#' change and likelihood not to have allelic fraction change.
#' Default: \code{-3}.
#'
#' @return a \code{matrix} of \code{numeric} with 3 columns where each
#' row represent a segment
#' of imbalanced SNVs. The first column represents the position, in
#' \code{snpPos}, of the first
#' SNV in the segment. The second column represents the position, in the
#' \code{snpPos}, of the last SNV in the segment. The third column represents
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
#' RAIDS:::computeAlleleFraction(snpPos=snpInfo, w=10, cutOff=-3)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom stats median
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeAlleleFraction <- function(snpPos, w=10, cutOff=-3) {

    listBlockAR <- list()

    tmp <- as.integer(snpPos$imbAR == 1)
    z <- cbind(c(tmp[1], tmp[-1] - tmp[seq_len(length(tmp) -1)]),
               c(tmp[-1] - tmp[seq_len(length(tmp) -1)], tmp[length(tmp)] * -1))

    ## Split SNVs by segment of continuous imbalanced SNVs
    ## There must be at least one segment with imbalanced SNVs to go one
    if(length(which(z[,1] == 1)) > 0) {
        ## Find segmentsof imbalanced SNVs
        segImb <- data.frame(start=seq_len(nrow(snpPos))[which(z[,1] > 0)],
                             end=seq_len(nrow(snpPos))[which(z[,2] < 0)])

        listBlockAR <- lapply(seq_len(nrow(segImb)),
                        FUN=function(i, segImb, snpPos, w,
                            cutOff){
                                listBlockAR <- list()
                                j<-1
                                listSeg <- (segImb$start[i]):(segImb$end[i])
                                # index hetero segment
                                listHetero  <- listSeg[snpPos[listSeg,"hetero"] == TRUE]
                                # SNP hetero for the segment
                                snp.hetero <- snpPos[listHetero,]

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
                                listBlockAR <- do.call(rbind, listBlockAR)
                                return(listBlockAR)
                            },
                        segImb=segImb, snpPos=snpPos,
                        w=w, cutOff=cutOff)
    }

    # note NULL if length(listBlockAR) == 0
    listBlockAR <- do.call(rbind, listBlockAR)
    # print(all.equal(listBlockAR, listBlockAR1))
    return(listBlockAR)
}


#' @title Estimate the allelic fraction of the pruned SNVs for a specific
#' DNA-seq profile
#'
#' @description The function creates a \code{data.frame} containing the
#' allelic fraction for the pruned SNV dataset specific to a DNA-seq profile
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Reference GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
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
#' the length of the chromosomes.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#'
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
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
#' the window used to compute an empty box. Default: \code{9L}.
#'
#' @param verbose a \code{logicial} indicating if the function should print
#' message when running.
#'
#' @return a \code{data.frame} containing the allelic information for the
#' pruned SNV dataset with coverage > \code{minCov}. The \code{data.frame}
#' contains those columns:
#' \itemize{
#' \item{cnt.tot} {a \code{integer} representing the total allele count}
#' \item{cnt.ref} {a \code{integer} representing the reference allele count}
#' \item{cnt.alt} {a \code{integer} representing the alternative allele count}
#' \item{snp.pos} {a \code{integer} representing the position on the chromosome}
#' \item{snp.chr} {a \code{integer} representing the chromosome}
#' \item{normal.geno} {a \code{integer} representing the genotype
#' (0=wild-type reference; 1=heterozygote; 2=homozygote alternative; 3=unkown)}
#' \item{pruned} {a \code{logical} indicating if the SNV is retained after
#' pruning}
#' \item{snp.index} {a \code{integer} representing the index position of the
#' SNV in the Reference GDS file that contains all SNVs}
#' \item{keep} {a \code{logical} indicating if the genotype exists for the SNV}
#' \item{hetero} {a \code{logical} indicating if the SNV is heterozygote}
#' \item{homo} {a \code{logical} indicating if the SNV is homozygote}
#' \item{lap} {a \code{numeric} indicating lower allelic fraction}
#' \item{LOH} {a \code{integer} indicating if the SNV is in an LOH region
#' (0=not LOH, 1=in LOH)}
#' \item{imbAR} {a \code{integer} indicating if the SNV is in an imbalanced
#' region (-1=not classified as imbalanced or LOH, 0=in LOH; 1=tested
#' positive for imbalance in at least 1 window)}
#' }
#'
#' @examples
#'
#' ## Required library for GDS
#' library(SNPRelate)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
#'
#' ## Temporary Profile GDS file for one profile
#' fileProfile <- file.path(getwd(), "ex1.gds")
#'
#' ## Example can only be run if the current directory is in writing mode
#' if (file.access(getwd()) == 0 && !file.exists(fileProfile)) {
#'
#'     ## Copy the Profile GDS file demo that has been pruned and annotated
#'     file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                  fileProfile)
#'
#'     ## Open the reference GDS file (demo version)
#'     gds1KG <- snpgdsOpen(fileGDS)
#'
#'     ## Open Profile GDS file for one profile
#'     profileGDS <- openfn.gds(fileProfile)
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'         181538259L, 170805979L, 159345973L, 145138636L, 138394717L,
#'         133797422L, 135086622L, 133275309L, 114364328L, 107043718L,
#'         101991189L, 90338345L, 83257441L,  80373285L,  58617616L,
#'         64444167L,  46709983L, 50818468L, 156040895L, 57227415L,  16569L)
#'
#'     ## The function returns a data frame containing the allelic fraction info
#'     result <- RAIDS:::computeAllelicFractionDNA(gdsReference=gds1KG,
#'         gdsSample=profileGDS,
#'         currentProfile="ex1", studyID="MYDATA", chrInfo=chrInfo, minCov=10L,
#'         minProb=0.999, eProb=0.001, cutOffLOH=-5,
#'         cutOffHomoScore=-3, wAR=9L, verbose=FALSE)
#'     head(result)
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
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeAllelicFractionDNA <- function(gdsReference, gdsSample, currentProfile,
                                studyID, chrInfo, minCov=10L, minProb=0.999,
                                eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                                wAR=9L, verbose) {

    ## Extract the genotype information for a SNV dataset using
    ## the Profile GDS file and the Reference GDS file
    snpPos <- getTableSNV(gdsReference=gdsReference, gdsSample=gdsSample,
        currentProfile=currentProfile, studyID=studyID, minCov=minCov,
        minProb=minProb, eProb=eProb, verbose=verbose)

    snpPos$lap <- rep(-1, nrow(snpPos))
    snpPos$LOH <- rep(0, nrow(snpPos))
    snpPos$imbAR <- rep(-1, nrow(snpPos))



    snpPos <- lapply(unique(snpPos$snp.chr),
            FUN=function(chr,snpPos){
                if (verbose) {
                    message("chr ", chr)
                    message("Step 1 ", Sys.time())
                }
                homoBlock <- list()
                listChr <- which(snpPos$snp.chr == chr)
                ## Identify LOH regions
                homoBlock[[chr]] <- computeLOHBlocksDNAChr(
                                        gdsReference=gdsReference,
                                        chrInfo=chrInfo,
                                        snpPos=snpPos[listChr,],
                                        chr=chr)

                if (verbose) { message("Step 2 ", Sys.time()) }

                homoBlock[[chr]]$LOH <- as.integer(homoBlock[[chr]]$logLHR <=
                        cutOffLOH & homoBlock[[chr]]$homoScore <= cutOffHomoScore)

                z <- cbind(c(homoBlock[[chr]]$start, homoBlock[[chr]]$end,
                            snpPos[listChr, "snp.pos"]),
                        c(rep(0,  2* nrow(homoBlock[[chr]])),
                            rep(1, length(listChr))),
                        c(homoBlock[[chr]]$LOH, -1 * homoBlock[[chr]]$LOH,
                            rep(0, length(listChr)) ),
                        c(rep(0, 2 * nrow(homoBlock[[chr]])),
                            seq_len(length(listChr))))

                z <- z[order(z[,1], z[,2]), ]
                pos <- z[cumsum(z[,3]) > 0 & z[,4] > 0, 4]
                snpPos[listChr[pos], "lap"] <- 0
                snpPos[listChr[pos], "LOH"] <- 1

                if (verbose) { message("Step 3 ", Sys.time()) }

                ## Identify imbalanced SNVs in specified chromosome
                snpPos[listChr, "imbAR"] <-
                    computeAllelicImbDNAChr(snpPos=snpPos[listChr, ], chr=chr,
                        wAR=10, cutOffEmptyBox=-3)

                if (verbose) { message("Step 4 ", Sys.time()) }

                ## Compute allelic fraction for SNVs in specified chromosome
                blockAF <- computeAlleleFraction(snpPos=snpPos[listChr, ],
                                                w=10, cutOff=-3)

                if (verbose) { message("Step 5 ", Sys.time()) }

                if(! is.null(blockAF)) {
                    for(i in seq_len(nrow(blockAF))) {
                        snpPos[listChr[blockAF[i, 1]:blockAF[i, 2]], "lap"] <-
                            blockAF[i, 3]
                    }
                }
                return(snpPos[listChr,])
            },
            snpPos=snpPos)
    snpPos <- do.call(rbind, snpPos)
    snpPos[which(snpPos[, "lap"] == -1), "lap"] <- 0.5

    return(snpPos)
}


#' @title Estimate the allelic fraction of the pruned SNVs for a specific
#' RNA-seq sample
#'
#' @description The function creates a \code{data.frame} containing the
#' allelic fraction for the pruned SNV dataset specific to a RNA-seq sample.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Reference GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Profile GDS file.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Reference SNV Annotation GDS file.
#'
#' @param currentProfile a \code{character} string corresponding to
#' the sample identifier as used in \code{\link{pruningSample}} function.
#'
#' @param studyID a \code{character} string corresponding to the name of
#' the study as
#' used in \code{\link{pruningSample}} function.
#'
#' @param blockID a \code{character} string corresponding to the field gene
#' block in the GDS \code{gdsRefAnnot} to use split by gene.
#'
#' @param chrInfo a \code{vector} of \code{integer} values representing
#' the length of the chromosomes.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#'
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#'
#' @param eProb a single \code{numeric} between 0 and 1 representing the
#' probability of sequencing error. Default: \code{0.001}.
#'
#' @param cutOffLOH a single \code{numeric} log of the score to be LOH.
#' Default: \code{-5}.
#'
#' @param cutOffAR a single \code{numeric} representing the cutoff, in
#' log score, to tag SNVs located in a gene has having an allelic fraction
#' different 0.5
#' Default: \code{3}.
#'
#' @param verbose a \code{logicial} indicating if the function should print
#' message when running.
#'
#' @return a \code{data.frame} containing the allelic information for the
#' pruned SNV dataset with coverage > \code{minCov}. The \code{data.frame}
#' contains those columns:
#' \itemize{
#' \item{cnt.tot} {a \code{integer} representing the total allele count}
#' \item{cnt.ref} {a \code{integer} representing the reference allele count}
#' \item{cnt.alt} {a \code{integer} representing the alternative allele count}
#' \item{snp.pos} {a \code{integer} representing the position on the chromosome}
#' \item{snp.chr} {a \code{integer} representing the chromosome}
#' \item{normal.geno} {a \code{integer} representing the genotype
#' (0=wild-type reference; 1=heterozygote; 2=homozygote alternative; 3=unkown)}
#' \item{pruned} {a \code{logical} indicating if the SNV is retained after
#' pruning}
#' \item{snp.index} {a \code{integer} representing the index position of the
#' SNV in the Reference GDS file that contains all SNVs}
#' \item{keep} {a \code{logical} indicating if the genotype exists for the SNV}
#' \item{hetero} {a \code{logical} indicating if the SNV is heterozygote}
#' \item{homo} {a \code{logical} indicating if the SNV is homozygote}
#' \item{block.id} {TOREVIEW a \code{integer} indicating the block.id in gdsRefAnnot the
#' vairant is in}
#' \item{phase} {TOREVIEW a \code{integer} indicating the phase of the variant
#' if known, 3 if not known}
#' \item{lap} {a \code{numeric} indicating lower allelic fraction}
#' \item{LOH} {a \code{integer} indicating if the SNV is in an LOH region
#' (0=not LOH, 1=in LOH)}
#' \item{imbAR} {a \code{integer} indicating if the SNV is in an imbalanced
#' region (-1=not classified as imbalanced or LOH, 0=in LOH; 1=tested
#' positive for imbalance in at least 1 window)}
#' \item{freq} {a \code{numeric} indicating the frequency of the variant
#' in the the reference}
#' }
#'
#' @examples
#'
#' ## Required library for GDS
#' library(SNPRelate)
#'
#' #' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
#' fileAnnotGDS <- file.path(dataDir, "ex1_good_small_1KG_Annot_GDS.gds")
#'
#' ## Temporary Profile GDS file for one profile
#' fileProfile <- file.path(getwd(), "ex1.gds")
#'
#' ## Example can only be run if the current directory is in writing mode
#' if (file.access(getwd()) == 0 && !file.exists(fileProfile)) {
#'
#'     ## Copy the Profile GDS file demo that has been pruned and annotated
#'     file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                  fileProfile)
#'
#'     ## Open the reference GDS file (demo version)
#'     gds1KG <- snpgdsOpen(fileGDS)
#'    gdsRefAnnot <- openfn.gds(fileAnnotGDS)
#'     ## Open Profile GDS file for one profile
#'     profileGDS <- openfn.gds(fileProfile)
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'         181538259L, 170805979L, 159345973L, 145138636L, 138394717L,
#'         133797422L, 135086622L, 133275309L, 114364328L, 107043718L,
#'         101991189L, 90338345L, 83257441L,  80373285L,  58617616L,
#'         64444167L,  46709983L, 50818468L, 156040895L, 57227415L,  16569L)
#'
#'     ## The function returns a data frame containing the allelic fraction info
#'     result <- RAIDS:::computeAllelicFractionRNA(gdsReference=gds1KG,
#'         gdsSample=profileGDS, gdsRefAnnot=gdsRefAnnot,
#'         currentProfile="ex1", studyID="MYDATA",
#'         blockID="GeneS.Ensembl.Hsapiens.v86",
#'         chrInfo=chrInfo, minCov=10L,
#'         minProb=0.999, eProb=0.001, cutOffLOH=-5,
#'         cutOffAR=3, verbose=FALSE)
#'     head(result)
#'
#'     ## Close both GDS files (important)
#'     closefn.gds(profileGDS)
#'     closefn.gds(gds1KG)
#'     closefn.gds(gdsRefAnnot)
#'
#'     ## Remove Profile GDS file (created for demo purpose)
#'     unlink(fileProfile, force=TRUE)
#'
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn ls.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeAllelicFractionRNA <- function(gdsReference, gdsSample, gdsRefAnnot,
                    currentProfile, studyID, blockID, chrInfo, minCov=10L,
                    minProb=0.999, eProb=0.001, cutOffLOH=-5,
                    cutOffAR=3, verbose) {

    ## Extract the genotype information for a SNV dataset using
    ## the GDS Sample file and the Reference GDS file
    snpPos <- getTableSNV(gdsReference, gdsSample, currentProfile, studyID,
                                minCov, minProb, eProb)
    # Keep only SNV in GDS ref because to reduce SNV artefact from RNA
    snpPos <- snpPos[which(snpPos$snp.index > 0),]

    # Get the block structure base on genes from gdsRefAnnot
    snpPos$block.id <- getBlockIDs(gdsRefAnnot=gdsRefAnnot,
                            snpIndex=snpPos$snp.index, blockTypeID=blockID)

    snpPos$phase <- rep(3, nrow(snpPos))
    if ("phase" %in% ls.gdsn(node=gdsSample)) {
        snpPos$phase <- read.gdsn(index.gdsn(gdsSample,
                                            "phase"))[snpPos$snp.index]
    }
    snpPos$lap <- rep(-1, nrow(snpPos))
    snpPos$LOH <- rep(0, nrow(snpPos))
    snpPos$imbAR <- rep(-1, nrow(snpPos))
    snpPos$freq <- read.gdsn(index.gdsn(gdsReference,
                                            "snp.AF"))[snpPos$snp.index]
    # for each chromosome
    # listBlock <- list()

    listBlock <- lapply(unique(snpPos$snp.chr), FUN = function(x, snpPos, verbose){
            if (verbose) {
                message("chr ", x)
                message("Step 1 ", Sys.time())
            }
            listChr <- which(snpPos$snp.chr == x)
            blockAF <- tableBlockAF(snpPos=snpPos[listChr,])
            blockAF$aRF[blockAF$lRhomo <= cutOffLOH] <- 0
            blockAF$aRF[blockAF$lR >= cutOffAR] <- blockAF$aFraction[blockAF$lR
                                                                     >= cutOffAR]
            blockAF$aRF[blockAF$lR < cutOffAR & blockAF$nbHetero > 1] <- 0.5

            #listBlock[[x]] <- blockAF

            if (verbose) {
                message("Step 1 done ", Sys.time())
            }
            return(blockAF)
        }, snpPos=snpPos,
        verbose=verbose)
    blockAF <- do.call(rbind, listBlock)
    listMissing <- which(abs(blockAF$aRF + 1) < 1e-6)

    if(length(listMissing) > 0){
        if(length(blockAF$aRF[-1*listMissing]) == 0){
            blockAF[listMissing, "aRF"] <- 0.5
        }else{
            blockAF[listMissing, "aRF"] <- sample(blockAF$aRF[-1*listMissing],
                                length(listMissing), replace=TRUE)
        }

    }

    for(b in seq_len(nrow(blockAF))) {
        snpPos$lap[snpPos$block.id == blockAF$block[b]] <- blockAF$aRF[b]
    }

    return(snpPos)
}


#' @title Verify if SNVs are in an imbalance region
#'
#' @description The function verifies, for each SNV present in the data frame,
#' if the SNV is in an imbalance region.
#'
#' @param snpPos a \code{data.frame} containing the SNV information for the
#' chromosome specified by the \code{chr} argument. The \code{data.frame} must
#' contain:
#' \itemize{
#' \item{cnt.tot} {a single \code{integer} representing the total coverage for
#' the SNV.}
#' \item{cnt.ref} {a single \code{integer} representing the coverage for
#' the reference allele.}
#' \item{cnt.alt} {a single \code{integer} representing the coverage for
#' the alternative allele.}
#' \item{snp.pos} {a single \code{integer} representing the SNV position.}
#' \item{snp.chr} {a single \code{integer} representing the SNV chromosome.}
#' \item{normal.geno} {a single \code{numeric} indicating the genotype of the
#' SNV. The possibles are: \code{0} (wild-type homozygote), \code{1}
#' (heterozygote), \code{2} (altenative homozygote), \code{3} indicating that
#' the normal genotype is unknown.}
#' \item{pruned} {a \code{logical} indicating if the SNV is retained after
#' pruning}
#' \item{snp.index} {a \code{integer} representing the index position of the
#' SNV in the Reference GDS file that contains all SNVs}
#' \item{keep} {a \code{logical} indicating if the genotype exists for the SNV}
#' \item{hetero} {a \code{logical} indicating if the SNV is heterozygote}
#' \item{homo} {a \code{logical} indicating if the SNV is homozygote}
#' \item{lap} {a \code{numeric} indicating lower allelic fraction}
#' \item{LOH} {a \code{integer} indicating if the SNV is in an LOH region
#' (0=not LOH, 1=in LOH)}
#' }
#'
#' @param chr a single positive \code{integer} for the current chromosome.
#'
#' @param wAR a single positive \code{integer} representing the size-1 of
#' the window used to compute an empty box. Default: \code{10}.
#'
#' @param cutOffEmptyBox a \code{numeric} representing the cut-off for
#' considering a region imbalanced when comparing likelihood to be imbalanced
#' and likelihood not to be imbalanced. Default: \code{-3}.
#'
#' @return a \code{vector} of \code{integer} indicating if the SNV is in an
#' imbalanced region (-1=not classified as imbalanced or LOH, 0=in LOH; 1=tested
#' positive for imbalance in at least 1 window). The vector as an entry for
#' each SNV present in the
#' input \code{snpPos}.
#'
#' @examples
#'
#' ## Required library for GDS
#' library(gdsfmt)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
#'
#' ## Open the reference GDS file (demo version)
#' gds1KG <- snpgdsOpen(fileGDS)
#'
#' ## Chromosome length information
#' ## chr23 is chrX, chr24 is chrY and chrM is 25
#' chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'     181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
#'     135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
#'     83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
#'     156040895L, 57227415L,  16569L)
#'
#' ## Data frame with SNV information for the specified chromosome (chr 1)
#' snpInfo <- data.frame(cnt.tot=c(41, 17, 27, 15, 11, 37, 16, 32),
#'     cnt.ref=c(40, 17, 27, 15, 4, 14, 16, 32),
#'     cnt.alt=c(0, 0, 0, 0, 7, 23, 0, 0),
#'     snp.pos=c(3722256, 3722328, 3767522, 3868160, 3869467, 4712655,
#'         6085318, 6213145),
#'     snp.chr=c(rep(1, 8)),
#'     normal.geno=c(rep(1, 8)), pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE,
#'     TRUE, TRUE),
#'     pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, rep(TRUE, 3)),
#'     snp.index=c(160, 162, 204, 256, 259, 288, 366, 465),
#'     keep=rep(TRUE, 8), hetero=c(rep(FALSE, 4), TRUE, TRUE, rep(FALSE, 2)),
#'     homo=c(rep(TRUE, 4), FALSE, FALSE, TRUE, TRUE),
#'     lap=rep(-1, 8), LOH=rep(0, 8), imbAR=rep(-1, 8),
#'     stringAsFactor=FALSE)
#'
#' ## The function returns a data frame containing the information about the
#' ## LOH regions in the specified chromosome
#' result <- RAIDS:::computeAllelicImbDNAChr(snpPos=snpInfo, chr=1, wAR=10,
#'     cutOffEmptyBox=-3)
#' head(result)
#'
#' ## Close GDS file (important)
#' closefn.gds(gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeAllelicImbDNAChr <- function(snpPos, chr, wAR=10, cutOffEmptyBox=-3) {

    # We use wAR - 1 because
    # process the window ex: 1 to 1+wAR
    wAR <- wAR - 1
    listHetero <- NULL
    if(length(which(snpPos$normal.geno != 3) > 0)) {
        listHetero <- which(snpPos$keep == TRUE & snpPos$normal.geno == 1)
    } else{
        listHetero <- which(snpPos$hetero == TRUE)
    }

    heteroSNV <- snpPos[listHetero, ]

    if(nrow(heteroSNV) > wAR) {
        for(i in seq_len(nrow(heteroSNV)-wAR)) {
            if(sum(snpPos[listHetero[i]:listHetero[(i+wAR-1)], "LOH"]) == 0 ) {
                ## Check for imbalance regions for heterozygote SNVs
                cur <- testEmptyBox(heteroSNV[i:(i+wAR), c
                                ("cnt.alt", "cnt.ref")], cutOffEmptyBox)
                if(cur$pCut == 1) {
                    # Set all snv from tmpA (include homozygotes)
                    # in the window  to 1
                    snpPos[listHetero[i]:listHetero[(i+wAR)], "imbAR"] <- 1
                }
            }
        }
    }
    snpPos$imbAR[which(snpPos$LOH == 1)] <- 0

    return(snpPos$imbAR)
}


#' @title Calculate presence of allelic fraction change using specified SNVs
#' separately and together
#'
#' @description The function tests allelic fraction change using all
#' specified SNVs
#' separately and together. The function reports the associated results,
#' including statistic for the region represented by all the SNVs.
#'
#' @param matCov a \code{data.frame} containing only heterozygote SNVs. The
#' \code{data.frame} must contain those columns:
#' \itemize{
#' \item{cnt.ref} {a single \code{integer} representing the coverage for
#' the reference allele.}
#' \item{cnt.alt} {a single \code{integer} representing the coverage for
#' the alternative allele.}
#' }
#'
#' @param pCutOff a \code{numeric} representing the cut-off for considering
#' a region imbalanced when comparing likelihood to gave allelic fraction
#' change and likelihood not to have allelic fraction change.
#' Default: \code{-3}.
#'
#' @param vMean a positive \code{numeric} representing the current ratio
#' (minor allele/(minor allele + second allele)) that is going to be used as
#' reference to see if there is a allelic fraction change.
#'
#' @return a \code{list} containing 4 entries:
#' \itemize{
#' \item{pWin}{ a \code{vector} of \code{numeric} representing the
#' probability (x2) of obtaining the current
#' alternative/(alternative+reference) ratio from a reference distribution
#' specified by user.}
#' \item{p}{a \code{integer} indicating if all SNVs tested
#' positive (1=TRUE, 0=FALSE). The cut-off is 0.5. }
#' \item{pCut}{a \code{integer} indicating if all SNVs tested
#' positive (1=TRUE, 0-FALSE). }
#' \item{pCut1}{a \code{integer} indicating if the region tested
#' positive (1=TRUE, 0=FALSE) for allelic ratio change.}
#' }
#'
#' @examples
#'
#' ## Data frame with SNV information
#' snpInfo <- data.frame(
#'     cnt.ref=c(40, 17, 27, 15, 4, 14, 16, 32),
#'     cnt.alt=c(2, 4, 5, 10, 7, 23, 0, 0))
#'
#' RAIDS:::testAlleleFractionChange(matCov=snpInfo, pCutOff=-3, vMean=0.5)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom stats pbinom
#' @encoding UTF-8
#' @keywords internal
testAlleleFractionChange <- function(matCov, pCutOff=-3, vMean) {
    p <- 0
    pO <- 0

    matCov$pWin <- rep(1, nrow(matCov))

    matTmp <- apply(matCov[, c("cnt.alt", "cnt.ref")], 1,
                    FUN=function(x, vMean){
                        vCur <- ifelse(x[1] <= x[2],
                                       x[1], x[2])

                        diff2Mean <- abs(vMean * (x[1] +
                                                x[2]) - vCur)
                        pCur1 <- pbinom(round(vMean * (x[1] +
                                                    x[2]) - diff2Mean),
                                            size=x[2] + x[1], vMean)
                        pCur2 <- 1 - pbinom(round(vMean * (x[1] +
                                                    x[2]) + diff2Mean),
                                            size=x[2] + x[1], vMean)

                        pCur <- pCur1 + pCur2

                        pCurO <- max(1 - max(pCur, 0.01), 0.01)
                        pCurMax <- max(pCur, 0.01)
                        return(c(pCur, pCurMax, pCurO))
                    }, vMean=vMean)
    matCov$pWin <- matTmp[1, ]

    p <- sum(log10(matTmp[2,]))
    p0 <- sum(log10(matTmp[3,]))


    pCut1 <- as.integer((sum(matCov$pWin < 0.5) >= nrow(matCov)-1) &
                            matCov$pWin[1] < 0.5 &
                            (matCov$pWin[nrow(matCov)] < 0.5)  &
                            ((p-pO) <= pCutOff))
    res <- list(pWin=matCov$pWin, p=p,
                pCut=as.integer(sum(matCov$pWin < 0.5) == nrow(matCov)),
                pCut1=pCut1)
    return(res)
}


#' @title Calculate imbalance region using specified heterozygote SNVs
#' separately and together
#'
#' @description The function tests imbalance region using all specified SNVs
#' separately and together. The function reports the associated results,
#' including statistic for the region.
#'
#' @param matCov a \code{data.frame} containing only heterozygote SNVs. The
#' \code{data.frame} must contain those columns:
#' \itemize{
#' \item{cnt.ref} {a single \code{integer} representing the coverage for
#' the reference allele.}
#' \item{cnt.alt} {a single \code{integer} representing the coverage for
#' the alternative allele.}
#' }
#'
#' @param pCutOff a \code{numeric} representing the cut-off for considering
#' a region imbalanced when comparing likelihood to be imbalanced and
#' likelihood not to be imbalanced. Default: \code{-3}.
#'
#' @return a \code{list} containing 4 entries:
#' \itemize{
#' \item{pWin}{ a \code{vector} of \code{numeric} representing the
#' probability (x2) of obtaining the current
#' alternative/(alternative+reference) ratio from a 0.5 distribution.}
#' \item{p}{ a \code{numeric} representing the likelihood for the region}
#' \item{pCut}{ a \code{integer} indicating if all SNVs tested
#' positive (1=TRUE, 0=FALSE). The cut-off is 0.5. }
#' \item{pCut1}{ a \code{integer} indicating if the window tested
#' positive (1=TRUE, 0=FALSE) for imbalance.}
#' }
#'
#' @examples
#'
#' ## Data frame with SNV information for a list of heterozygote SNVs
#' snpInfo <- data.frame(
#'     cnt.ref=c(40, 17, 27, 15, 4, 14, 16, 32),
#'     cnt.alt=c(2, 4, 5, 10, 7, 23, 0, 0))
#'
#' ## Calculate imbalance for the region represented by the SNVs
#' RAIDS:::testEmptyBox(matCov=snpInfo, pCutOff=-3)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom stats pbinom
#' @encoding UTF-8
#' @keywords internal
testEmptyBox <- function(matCov, pCutOff=-3) {

    p <- 0
    pO <- 0
    vMean <- 0.5
    matCov$pWin <- rep(1, nrow(matCov))

    matTmp <- apply(matCov[, c("cnt.alt", "cnt.ref")], 1,
                FUN=function(x){
                    vCur1 <- ifelse(x[1] <= x[2],
                        x[1], x[2])

                    pCur <- pbinom(q=vCur1, size=x[2] + x[1],
                        prob=vMean)

                    pCurO <- max(1 - max(2 * pCur, 0.01), 0.01)

                    return(c(vCur1, pCur, pCurO))
                })
    matCov$pWin <- matTmp[2, ] * 2
    tmp <- which(matTmp[2,] > 0.01)
    p <- sum(log10(matTmp[2,tmp])) + (ncol(matTmp) - length(tmp)) * log10(0.01)
    p0 <- sum(log10(matTmp[3,]))


    ## Calculate a global statistic using all SNVs
    ## The region is imbalance or not
    pCut1 <- as.integer((sum(matCov$pWin < 0.5) >= nrow(matCov)-1) &
                matCov$pWin[1] < 0.5 & (matCov$pWin[nrow(matCov)] < 0.5) &
                            ((p-pO) <= pCutOff))

    ## Return a list
    res <- list(pWin=matCov$pWin, p=p,
                pCut=as.integer(sum(matCov$pWin < 0.5) == nrow(matCov)),
                pCut1=pCut1)
    return(res)
}


###############################################
# RNA-seq section
###############################################


#' @title TOREVIEW Compute the log likelihood ratio base on the coverage (read depth)
#' of each allele in block (gene in the case of RNA-seq)
#'
#' @description TOREVIEW For the block sum the log of read depth of the lowest depth
#' divide by the total depth of the position minus of likelhood of the allelic
#' fraction of 0.5. If the phase is known, the variant varaint in the same
#' haplotype are group.
#'
#' @param snpPosHetero For a specific gene (block) a \code{data.frame}
#' containing the SNV information.
#' The \code{data.frame} must contain those columns:
#' \itemize{
#' \item{cnt.ref} {a single \code{integer} representing the coverage for
#' the reference allele.}
#' \item{cnt.alt} {a single \code{integer} representing the coverage for
#' the alternative allele.}
#' \item{phase} {TOREVIEW a \code{integer} indicating the phase of the variant
#' if known, 3 if not known}
#' }
#'
#' @return TOREVIEW a \code{list}  for the block with the information of
#' relative to the heterozygotes.
#' The \code{list} contains:
#' \itemize{
#' \item{lR} {TOREVIEW a single \code{numeric} representing sum the log of read depth of the lowest depth
#' divide by the total depth of the position minus of likelhood of the allelic
#' fraction of 0.5.}
#' \item{aFraction} {TOREVIEW a single \code{numeric} representing the allele
#' fraction estimation.}
#' \item{sumAlleleLow} {TOREVIEW a \code{integer} representing the
#' sum of the allele read depth
#' of the lowest read alelle depth}
#' \item{sumAlleleHigh} {TOREVIEW a \code{integer} representing the
#' sum of the allele read depth
#' of the highsest read alelle depth}
#' }
#'
#' @examples
#'
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' snpPos <- readRDS(file.path(dataDir, "demoAllelicFraction", "demSnpPos.rds"))
#'
#' result <- RAIDS:::calcAFMLRNA(snpPos[which(
#'                 snpPos$block.id == 2750 &
#'                 snpPos$hetero), c("cnt.ref",
#'                     "cnt.alt", "phase")])
#' head(result)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
calcAFMLRNA <- function(snpPosHetero) {

    listPhase <- which(snpPosHetero$phase < 2)
    m <- data.frame(aL=rep(0, nrow(snpPosHetero)),
                        aH=rep(0, nrow(snpPosHetero)))
    # For the vairants phase
    # sum the coverage for each haplotype
    #
    if (length(listPhase) > 0) {
        mPhase <- data.frame(a1=rep(0, length(listPhase)),
                                a2=rep(0, length(listPhase)))
        if(length(which(snpPosHetero$phase == 0)) > 0){
            mPhase[which(snpPosHetero$phase == 0), "a1"] <-
                snpPosHetero[which(snpPosHetero$phase == 0),"cnt.ref"]
            mPhase[which(snpPosHetero$phase == 0), "a2"] <-
                snpPosHetero[which(snpPosHetero$phase == 0),"cnt.alt"]
        }
        if(length(which(snpPosHetero$phase == 1)) > 0){
            mPhase[which(snpPosHetero$phase == 1), "a2"] <-
                snpPosHetero[which(snpPosHetero$phase == 1),"cnt.ref"]
            mPhase[which(snpPosHetero$phase == 1), "a1"] <-
                snpPosHetero[which(snpPosHetero$phase == 1),"cnt.alt"]
        }

        m1 <- sum(mPhase[,"a1"])
        m2 <- sum(mPhase[,"a2"])
        minPhase <- which.min(c(m1,m2))
        m[listPhase, "aL"] <- mPhase[, minPhase]
        m[listPhase, "aH"] <- mPhase[, (minPhase+1)%%2]
    }

    listUnphase <- which(snpPosHetero$phase > 1)

    if(length(listUnphase) > 0){
        minUnphase <- apply(snpPosHetero[,c("cnt.ref", "cnt.alt")], 1,
                                FUN=min)
        maxUnphase <- apply(snpPosHetero[,c("cnt.ref", "cnt.alt")], 1,
                                FUN=max)
        m[listUnphase, "aL"] <- minUnphase
        m[listUnphase, "aH"] <- maxUnphase
    }
    d <- sum(rowSums(snpPosHetero[,c("cnt.ref", "cnt.alt")]))

    aF <- sum(m[,"aL"]) / d
    lM <- log10(aF) * sum(m[,"aL"]) + log10(1- aF) * sum(m[,"aH"])
    lR <- lM - log10(0.5) * d
    res <- list(lR = lR, aFraction=aF, nPhase = length(listPhase),
                    sumAlleleLow = sum(m[,"aL"]),
                    sumAlleleHigh = sum(m[,"aH"]))
    return(res)
}


#' @title TOREVIEW Compile the information about the variants
#' for each bloc.
#'
#' @description TOREVIEW For each block, the function evaluates a score
#' about lost of heterizygocity and allelic fration. It generates some
#' information about the variants in the block like number of homozigte or
#' heterozygote.
#' In the case of RNA-seq the blocks are genes.
#'
#' @param snpPos For a specific chromosome a \code{data.frame} with lap for
#' the SNV dataset with
#' coverage > \code{minCov}.
#'
#' @return TOREVIEW resBlock a \code{data.frame} containing only heterozygote
#' SNV information. The
#' \code{data.frame} contain those columns:
#' \itemize{
#' \item{block} {a single \code{integer} representing the unique identifier
#' of the block.}
#' \item{aRF} {a single \code{numeric} representing the final allelic fraction not compute
#' here -1 for all entries.}
#' \item{aFraction} {a single \code{integer} representing the possible allelic
#' fraction in absence of LOH.}
#' \item{lR} {a single \code{integer} representing the coverage for
#' the alternative allele.}
#' \item{nPhase} {a single \code{integer} representin the number of SNV
#' phases.}
#' \item{sumAlleleLow} {a single \code{integer} TOREVIEW sum of the allele with
#' the less coverage.}
#' \item{sumAlleleHigh} {a single \code{integer} TOREVIEW sum of the allele
#' with more coverage.}
#' \item{lH} {a single \code{numeric} for the homozygotes log10 of the product
#' frequencies of the allele not found in the profile (not a probability).}
#' \item{lM} {a single \code{numeric} log10 product frequency allele
#' in population.}
#' \item{lRhomo} {a single \code{numeric} score lH - lM.}
#' \item{nbHomo} {a single \code{integer} representing the number of
#' homozygote SNVs per block.}
#' \item{nbKeep} {a single \code{integer} representing the number of
#' SNVs retained per block.}
#' \item{nbHetero} {a single \code{integer} representing the number of
#' heterozygote SNVs per block.}
#' }
#'
#' @examples
#'
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' snpPos <- readRDS(file.path(dataDir, "demoAllelicFraction", "demSnpPos.rds"))
#'
#' result <- RAIDS:::tableBlockAF(snpPos[which(snpPos$snp.chr == 1),])
#' head(result)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors aggregate
#' @encoding UTF-8
#' @keywords internal
tableBlockAF <- function(snpPos) {

    listBlocks <- unique(snpPos$block.id)

    resBlock <- data.frame(block=listBlocks)


    # Number of homozygotes per block
    tmp <- aggregate(snpPos[, c( "homo"), drop=FALSE],
                        by = list(block=snpPos$block.id), sum)
    row.names(tmp) <- as.character(tmp[, 1])
    resBlock$nbHomo <- tmp[as.character(listBlocks), 2]

    # Number of SNVs keep per block
    tmp <- aggregate(snpPos[, c( "keep"), drop=FALSE],
                        by = list(block=snpPos$block.id), sum)
    row.names(tmp) <- as.character(tmp[, 1])
    resBlock$nbKeep <- tmp[as.character(listBlocks), 2]

    # Number of heterozygotes per block
    tmp <- aggregate(snpPos[, c( "hetero"), drop=FALSE],
                        by = list(block=snpPos$block.id), sum)
    row.names(tmp) <- as.character(tmp[, 1])
    resBlock$nbHetero <- tmp[as.character(listBlocks), 2]

    resBlock <- apply(resBlock, 1, FUN=function(x, snpPos) {
        resBlock <- data.frame(block=x[1],
                               nbHomo=x[2],
                               nbKeep=x[3],
                               nbHetero=x[4],
                               aRF=-1,
                               aFraction=-1,
                               lR=-1,
                               nPhase=-1,
                               sumAlleleLow=-1,
                               sumAlleleHigh=-1,
                               lH=-1,
                               lM=-1,
                               lRhomo=1)


            lH <- 1
            lM <- 1
            if (resBlock[1, "nbKeep"] > 0 &
                (resBlock[1, "nbKeep"] == resBlock[1, "nbHomo"] |
                 (resBlock[1, "nbHomo"] > 0 & resBlock[1, "nbHetero"] == 1))) {

                # Check if 1 hetero with allelic fraction (<=0.05)
                # it is considered as all homozygote
                flag <- TRUE
                if (resBlock[1, "nbHetero"] == 1) {
                    tmp <- min(snpPos[snpPos$block.id == resBlock$block[1] &
                                    snpPos$hetero, c("cnt.ref" , "cnt.alt")])/
                        sum(snpPos[snpPos$block.id == resBlock$block[1] &
                                    snpPos$hetero, c("cnt.ref" , "cnt.alt")])

                    # flag is true if allelic fraction <= 0.05
                    flag <- ifelse(tmp > 0.05, FALSE, TRUE)
                }

                if(flag) {
                    # List homozygote ref
                    listRef <- which(snpPos$block.id == resBlock$block[1] &
                                         snpPos$homo &
                                         snpPos$cnt.ref > snpPos$cnt.alt)
                    # list homozygote alt
                    listAlt <- which(snpPos$block.id == resBlock$block[1] &
                                         snpPos$homo &
                                         snpPos$cnt.ref < snpPos$cnt.alt)
                    # freq of the Ref allele in population of listRef
                    tmp <- snpPos$freq[listRef]
                    ## min freq is 0.01
                    tmp[which(tmp < 0.01)] <- 0.01
                    ## log10 of the product of the frequency of the alternative
                    ## allele in pop for listRef
                    lH <- ifelse(length(listRef) > 0, sum(log10(1-tmp)*2), 0)
                    ## freq of the Ref allele in population of listAlt
                    tmp <- snpPos$freq[listAlt]
                    tmp[which(tmp < 0.01)] <- 0.01
                    ## log10 of the product of the frequency of the
                    ## alternative allele in pop for listRef
                    ## plus log10 of the product of the frequency of
                    ##the reference allele in pop for listAlt
                    lH <- lH + ifelse(length(listAlt) > 0,
                                        sum(log10(tmp)*2), 0)

                    lM <- sum(log10(apply(snpPos[which(snpPos$block.id ==
                                            resBlock$block[1] & snpPos$homo),
                                                  "freq", drop=FALSE], 1,
                                FUN=function(x) {
                                    return(max(x^2, 2*(x * (1-x)), (1-x)^2))
                                          })))
                    resBlock$sumAlleleLow[1] <- 0
                    resBlock$sumAlleleHigh[1] <- sum(snpPos[listRef, "cnt.ref"]) +
                        sum(snpPos[listAlt, "cnt.alt"])
                }
            }

            ## compute the score of the homozygote on the block
            ## if heterozygote present lH = lM = 1 and lRhomo = 0
            resBlock[1, c("lH", "lM", "lRhomo")] <- c(lH, lM, lH-lM)

            ## get hetero and compute AF nbHetero > 1
            if (resBlock[1, "nbKeep"] > 0 & resBlock[1, "nbHetero"] > 1) {

                resML <- calcAFMLRNA(snpPos[which(snpPos$block.id ==
                                        resBlock$block[1] & snpPos$hetero),])

                resBlock$aFraction[1] <- resML$aFraction
                resBlock$lR[1] <- resML$lR
                resBlock$nPhase[1] <- resML$nPhase
                resBlock$sumAlleleLow[1] <- resML$sumAlleleLow
                resBlock$sumAlleleHigh[1] <- resML$sumAlleleHigh
            }
            return(resBlock)
        }, snpPos=snpPos)
    resBlock <- do.call(rbind, resBlock)

    return(resBlock)
}


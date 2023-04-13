#' @title TODO
#'
#' @description TODO
#'
#' @param gdsReference an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param chrInfo a \code{vector} of \code{integer} representing the length of
#' the chromosomes. As an example, the information ca be obtained from
#' package 'BSgenome.Hsapiens.UCSC.hg38'.
#'
#' @param snp.pos a \code{data.frame} containing TODO.
#'
#' @param chr a single positive \code{integer} for the current chromosome. The
#' \code{chrInfo} parameter must contain the value for the specified
#' chromosome.
#'
#' @param  genoN a single \code{numeric} between 0 and 1 representing TODO.
#' Default: \code{0.0001}.
#'
#' @return a \code{data.frame} containing:
#' \itemize{
#' \item{chr} {TODO}
#' \item{start} {TODO}
#' \item{end} {TODO}
#' \item{logLHR} {TODO}
#' \item{LH1} {TODO}
#' \item{LM1} {TODO}
#' \item{homoScore} {TODO}
#' \item{nbSNV} {TODO}
#' \item{nbPruned} {TODO}
#' \item{nbNorm} {TODO}
#' \item{LOH} {TODO}
#' }
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Chromosome length information
#' ## chr23 is chrX, chr24 is chrY and chrM is 25
#' chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'     181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
#'     135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
#'     83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
#'     156040895L, 57227415L,  16569L)
#'
#' ## A formal way to get the chormosome length information
#' ## library(BSgenome.Hsapiens.UCSC.hg38)
#' ## chrInfo <- integer(25L)
#' ## for(i in seq_len(22L)){ chrInfo[i] <-
#' ##                          length(Hsapiens[[paste0("chr", i)]])}
#' ## chrInfo[23] <- length(Hsapiens[["chrX"]])
#' ## chrInfo[24] <- length(Hsapiens[["chrY"]])
#' ## chrInfo[25] <- length(Hsapiens[["chrM"]])
#'
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats dbinom
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeLOHBlocksDNAChr <- function(gdsReference, chrInfo, snp.pos, chr,
                                    genoN=0.0001) {

    ## The chr parameter must be a single integer value
    if (!isSingleNumber(chr))  {
        stop("The \'chr\' must be a single integer value representing ",
            "a chromosome")
    }

    ## The genoN parameter must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(genoN) && (genoN >= 0.0) && (genoN <= 1.0)))  {
        stop("The \'genoN\' must be a single numeric positive ",
                "value between 0 and 1.")
    }

    ## The specified chromosome is not included in the chrInfo parameter
    if (is.na(chrInfo[chr])) {
        stop("The \'chr\' must be present in the \'chrInfo\' parameter.")
    }

    ## The snp.pos must be a data.frame
    if (!is.data.frame(snp.pos)) {
        stop("The \'snp.pos\' must be a data.frame.")
    }

    genoN1 <- 1 - 2 * genoN

    chrEnd <- chrInfo[chr]
    listHetero <- snp.pos[snp.pos$hetero == TRUE, "snp.pos"]

    homoBlock <- data.frame(chr=rep(chr, length(listHetero) + 1),
                                start=c(1, listHetero + 1),
                                end=c(listHetero, chrEnd))

    z <- cbind(c(homoBlock$start, homoBlock$end,
                        snp.pos$snp.pos[which(snp.pos$homo == TRUE)]),
                c(seq_len(length(homoBlock$start)),
                        -1*seq_len(length(homoBlock$start)),
                        rep(0, length(which(snp.pos$homo == TRUE)))),
                c(rep(0, length(homoBlock$start)),
                        rep(0, length(homoBlock$start)),
                        seq_len(length(which(snp.pos$homo == TRUE)))))

    z <- z[order(z[, 1]), ]

    blcSNV <- data.frame(block = cumsum(z[, 2])[z[, 2] == 0],
                snv = z[z[, 2] == 0, 3])

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

    for (i in seq_len(nrow(homoBlock))) {
        blcCur <- blcSNV[blcSNV$block == i, ]
        snvH <- snp.pos[blcCur$snv, ]
        lH1 <- 0
        lM1 <- 0
        logLHR <- 0
        homoBlock$nbSNV[i] <- nrow(blcCur)
        homoBlock$nbPruned[i] <- length(which(snvH$pruned))
        if (length(which(snvH$normal.geno != 3)) > 0) {

            listCount <- snvH$cnt.tot[which(snvH$normal.geno == 1)]
            homoBlock$nbNorm[i] <- length(listCount)

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
    } # end for each block

    return(homoBlock)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param matCov TODO
#'
#' @param pCutOff TODO, Default: \code{-3}.
#'
#'
#' @return a \code{list} containing 4 entries:
#' \itemize{
#' \item{pWin}{TODO}
#' \item{p}{TODO}
#' \item{pCut}{TODO}
#' \item{pCut1}{TODO}
#' }
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats pbinom
#' @encoding UTF-8
#' @keywords internal
testEmptyBox <- function(matCov, pCutOff=-3) {

    p <- 0
    pO <- 0
    vMean <- 0.5
    matCov$pWin <- rep(1, nrow(matCov))

    for (i in seq_len(nrow(matCov))) {

        vCur1 <- ifelse(matCov$cnt.alt[i] <= matCov$cnt.ref[i],
                            matCov$cnt.alt[i], matCov$cnt.ref[i])

        pCur <- pbinom(q=vCur1, size=matCov$cnt.ref[i] + matCov$cnt.alt[i],
                        prob=vMean)

        pCurO <- max(1 - max(2 * pCur,0.01),0.01)

        matCov$pWin[i] <- pCur * 2

        p <- p + log10(max(pCur,0.01))
        pO <- pO + log10(pCurO)
    }

    pCut1 <- as.integer((sum(matCov$pWin < 0.5) >= nrow(matCov)-1) &
                                matCov$pWin[1] < 0.5 &
                                (matCov$pWin[nrow(matCov)] < 0.5) &
                                ((p-pO) <= pCutOff))
    res <- list(pWin=matCov$pWin, p=p,
                    pCut=as.integer(sum(matCov$pWin < 0.5) == nrow(matCov)),
                    pCut1=pCut1)
    return(res)
}



#' @title TODO
#'
#' @description TODO
#'
#' @param snp.pos a \code{data.frame} containing the genotype information for
#' a SNV dataset. TODO
#'
#' @param chr a single positive \code{integer} for the chromosome.
#'
#' @param w a single positive \code{numeric} representing the size of the
#' window to compute the allelic fraction.
#' Default: \code{10}.
#'
#' @param cutOff a single \code{numeric} representing TODO. Default: \code{-3}.
#'
#' @return TODO
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats median
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeAlleleFraction <- function(snp.pos, chr, w=10, cutOff=-3) {

    listBlockAR <- list()
    j <- 1
    tmp <- as.integer(snp.pos$imbAR == 1)
    z <- cbind(c(tmp[1], tmp[-1] - tmp[seq_len(length(tmp) -1)]),
            c(tmp[-1] - tmp[seq_len(length(tmp) -1)], tmp[length(tmp)] * -1))

    if(length(which(z[,1] == 1)) > 0) {

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
                # I am here
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


###############################################
# RNA-seq section
###############################################


#' @title TODO
#'
#' @description TODO
#'
#' @param snp.pos.Hetero For a specific gene (block) a \code{data.frame} with
#' lap for the SNV heterozygote dataset with
#' coverage > \code{minCov}. The \code{data.frame} must contain those columns:
#' 'phase', 'cnt.ref', 'cnt.alt'. TODO
#'
#' @return TODO a \code{list} of \code{numeric} for the gene lR the score
#' for aFraction different than 0.5
#' aFraction allele estimation, nPhase number of SNV phase,
#' sumAlleleLow number of read overlapping the allele low
#' sumAlleleHigh number of read overlapping the allele high TODO
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
calcAF.MLRNA <- function(snp.pos.Hetero) {

    listPhase <- which(snp.pos.Hetero$phase < 2)
    m <- data.frame(aL = rep(0, nrow(snp.pos.Hetero)),
                    aH = rep(0, nrow(snp.pos.Hetero)))
    if (length(listPhase) > 0) {
        mPhase <- data.frame(a1 = rep(0, length(listPhase)),
                                a2 = rep(0, length(listPhase)))
        if(length(which(snp.pos.Hetero$phase == 0)) > 0){
            mPhase[which(snp.pos.Hetero$phase == 0), "a1"] <-
                    snp.pos.Hetero[which(snp.pos.Hetero$phase == 0),"cnt.ref"]
            mPhase[which(snp.pos.Hetero$phase == 0), "a2"] <-
                    snp.pos.Hetero[which(snp.pos.Hetero$phase == 0),"cnt.alt"]
        }
        if(length(which(snp.pos.Hetero$phase == 1)) > 0){
            mPhase[which(snp.pos.Hetero$phase == 1), "a2"] <-
                    snp.pos.Hetero[which(snp.pos.Hetero$phase == 1),"cnt.ref"]
            mPhase[which(snp.pos.Hetero$phase == 1), "a1"] <-
                    snp.pos.Hetero[which(snp.pos.Hetero$phase == 1),"cnt.alt"]
        }

        m1 <- sum(mPhase[,"a1"])
        m2 <- sum(mPhase[,"a2"])
        minPhase <- which.min(c(m1,m2))
        m[listPhase, "aL"] <- mPhase[, minPhase]
        m[listPhase, "aH"] <- mPhase[, (minPhase+1)%%2]
    }

    listUnphase <- which(snp.pos.Hetero$phase > 1)

    if(length(listUnphase) > 0){
        minUnphase <- apply(snp.pos.Hetero[,c("cnt.ref", "cnt.alt")], 1,
                            FUN=min)
        maxUnphase <- apply(snp.pos.Hetero[,c("cnt.ref", "cnt.alt")], 1,
                            FUN=max)
        m[listUnphase, "aL"] <- minUnphase
        m[listUnphase, "aH"] <- maxUnphase
    }
    d <- sum(rowSums(snp.pos.Hetero[,c("cnt.ref", "cnt.alt")]))

    aF <- sum(m[,"aL"]) / d
    lM <- log10(aF) * sum(m[,"aL"]) + log10(1- aF) * sum(m[,"aH"])
    lR <- lM - log10(0.5) * d
    res <- list(lR = lR, aFraction=aF, nPhase = length(listPhase),
                sumAlleleLow = sum(m[,"aL"]),
                sumAlleleHigh = sum(m[,"aH"]))
    return(res)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param snp.pos For a specific chromosome a \code{data.frame} with lap for
#' the SNV dataset with
#' coverage > \code{minCov}.
#'
#' @return TODO a \code{data.frame} with the information related to allelic
#' fraction for each block gene
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors aggregate
#' @encoding UTF-8
#' @keywords internal
tableBlockAF <- function(snp.pos) {

    listBlocks <- unique(snp.pos$block.id)

    resBlock <- data.frame(block = listBlocks,
                            aRF = rep(-1, length(listBlocks)),
                            aFraction = rep(-1, length(listBlocks)),
                            lR = rep(-1, length(listBlocks)),
                            nPhase = rep(-1, length(listBlocks)),
                            sumAlleleLow = rep(-1, length(listBlocks)),
                            sumAlleleHigh = rep(-1, length(listBlocks)),
                            lH = rep(-1, length(listBlocks)),
                            lM = rep(-1, length(listBlocks)),
                            lRhomo = rep(1, length(listBlocks)))

    tmp <- aggregate(snp.pos[, c( "homo"), drop = FALSE],
                        by = list(block=snp.pos$block.id) ,sum)
    row.names(tmp) <- as.character(tmp[,1])
    resBlock$nbHomo <- tmp[as.character(listBlocks),2]
    tmp <- aggregate(snp.pos[, c( "keep"), drop = FALSE],
                        by = list(block=snp.pos$block.id) ,sum)
    row.names(tmp) <- as.character(tmp[,1])
    resBlock$nbKeep <- tmp[as.character(listBlocks),2]

    tmp <- aggregate(snp.pos[, c( "hetero"), drop = FALSE],
                        by = list(block=snp.pos$block.id) ,sum)
    row.names(tmp) <- as.character(tmp[,1])
    resBlock$nbHetero <- tmp[as.character(listBlocks),2]


    for (i in seq_len(length(listBlocks))) {
        # start with LOH

        lH <- 1
        lM <- 1
        if (resBlock[i, "nbKeep"] > 0 &
            (resBlock[i, "nbKeep"] == resBlock[i, "nbHomo"] |
            (resBlock[i, "nbHomo"] > 0 & resBlock[i, "nbHetero"] == 1)) ) {

            # Check if 1 hetero with allelic fraction (<=0.05)
            # it is considered as all homozygote
            flag <- TRUE
            if (resBlock[i, "nbHetero"] == 1) {
                tmp <- min(snp.pos[snp.pos$block.id == resBlock$block[i] &
                                snp.pos$hetero, c("cnt.ref" , "cnt.alt")])/
                            sum(snp.pos[snp.pos$block.id == resBlock$block[i] &
                                snp.pos$hetero, c("cnt.ref" , "cnt.alt")])
                flag <- ifelse(tmp > 0.05, FALSE,TRUE)
            }
            if(flag){
                listRef <- which(snp.pos$block.id == resBlock$block[i] &
                                        snp.pos$homo &
                                        snp.pos$cnt.ref > snp.pos$cnt.alt)
                listAlt <- which(snp.pos$block.id == resBlock$block[i] &
                                        snp.pos$homo &
                                        snp.pos$cnt.ref < snp.pos$cnt.alt)
                tmp <- snp.pos$freq[listRef]
                tmp[which(tmp < 0.01)] <- 0.01
                lH <- ifelse(length(listRef) > 0, sum(log10(1-tmp)*2), 0)

                tmp <- snp.pos$freq[listAlt]
                tmp[which(tmp < 0.01)] <- 0.01
                lH <- lH + ifelse(length(listAlt) > 0, sum(log10(tmp)*2), 0)

                lM <- sum(log10(apply(snp.pos[which(snp.pos$block.id ==
                                            resBlock$block[i] & snp.pos$homo),
                                            "freq", drop=FALSE], 1,
                                FUN = function(x) {
                                    return(max(x^2, 2*(x * (1-x)), (1-x)^2))
                                })))
                resBlock$sumAlleleLow[i] <- 0
                resBlock$sumAlleleHigh[i] <- sum(snp.pos[listRef, "cnt.ref"]) +
                    sum(snp.pos[listAlt, "cnt.alt"])
            }
        }

        resBlock[i, c("lH", "lM", "lRhomo")] <- c(lH, lM, lH - lM)

        # get hetero and compute AF
        if (resBlock[i, "nbKeep"] > 0 & resBlock[i, "nbHetero"] > 1) {

            resML <- calcAF.MLRNA(snp.pos[which(snp.pos$block.id ==
                                                        resBlock$block[i] &
                                                        snp.pos$hetero),])

            resBlock$aFraction[i] <- resML$aFraction
            resBlock$lR[i] <- resML$lR
            resBlock$nPhase[i] <- resML$nPhase
            resBlock$sumAlleleLow[i] <- resML$sumAlleleLow
            resBlock$sumAlleleHigh[i] <- resML$sumAlleleHigh

        }
    }

    return(resBlock)
}



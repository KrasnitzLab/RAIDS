#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param gdsSample TODO
#'
#' @param sampleCurrent A \code{character} string corresponding to
#' the sample.id used in LDpruning.
#'
#' @param study.id A \code{string} corresponding to the study
#' use in LDpruning.
#'
#' @param minCov a single positive \code{integer} representing the
#' minimum coverage needed to retain TODO. Default: \code{10}.
#'
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing TODO.
#' Default: \code{0.999}.
#'
#' @param eProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error. Default: \code{0.001}.
#'
#' @return a \code{data.frame} containing:
#' \itemize{
#' \item{cnt.tot} {a single \code{integer} representing the total coverage for
#' the SNV.}
#' \item{cnt.ref} {a single \code{integer} representing the coverage for
#' the reference allele.}
#' \item{cnt.alt} {a single \code{integer} representing the coverage for
#' the alternative allele.}
#' \item{snp.pos} {a single \code{integer} representing the SNV position.}
#' \item{snp.chr} {a single \code{integer} representing the SNV chromosome.}
#' \item{normal.geno} {a single \code{numeric} \code{3} indicating that
#' the normal genotype is unknown. TODO}
#' \item{snp.index} {The \code{boolean} \code{FALSE} indicating TODO}
#' }
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @encoding UTF-8
#' @export
getTableSNV <- function(gds, gdsSample, sampleCurrent, study.id, minCov=10,
                            minProb=0.999, eProb=0.001) {

    ## The minCov must be a single positive number
    if (!(isSingleNumber(minCov) && (minCov >= 0)))  {
        stop("The \'minCov\' must be a single numeric positive value.")
    }

    ## The minProb must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(minProb) && (minProb >= 0.0) && (minProb <= 1.0)))  {
        stop("The \'minProb\' must be a single numeric positive ",
                "value between 0 and 1.")
    }

    ## The eProb must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(eProb) && (eProb >= 0.0) && (eProb <= 1.0)))  {
        stop("The \'eProb\' must be a single numeric positive ",
             "value between 0 and 1.")
    }

    ## Extract study information from GDS
    study.annot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    ## Retain the specified sample in the specified study
    posCur <- which(study.annot$data.id == sampleCurrent &
                        study.annot$study.id == study.id)

    # g <- read.gdsn(index.gdsn(gdsSample, "geno.ref"), start=c(1, posCur), count = c(-1,1))[listSNP]
    # g <- read.gdsn(index.gdsn(gds, "genotype"), start=c(1,i), count = c(-1,1))[listSNP]


    ## Extract SNV coverage from GDS file
    cnt.total <- read.gdsn(node=index.gdsn(gdsSample, "Total.count"),
                            start=c(1, posCur), count=c(-1,1))
    #read.gdsn(index.gdsn(gdsSample, "Total.count"))

    ## Only retained the SNV with the minimum required coverage
    listKeep <- cnt.total@i[which(cnt.total@x >= minCov)] + 1

    ## Create the data.frame with the required information
    snp.pos <- data.frame(cnt.tot=cnt.total[listKeep],
                    cnt.ref=read.gdsn(index.gdsn(gdsSample, "Ref.count"),
                                            start=c(1, posCur),
                                            count=c(-1,1))[listKeep],
                    cnt.alt=read.gdsn(index.gdsn(gdsSample, "Alt.count"),
                                            start=c(1, posCur),
                                            count=c(-1,1))[listKeep],
                    snp.pos=read.gdsn(index.gdsn(gds,
                                            "snp.position"))[listKeep],
                    snp.chr=read.gdsn(index.gdsn(gds,
                                            "snp.chromosome"))[listKeep],
                    normal.geno=rep(3,length(listKeep)), # Suppose the normal genotype unkown
                    pruned=rep(FALSE, length(listKeep)),#bit(length(listKeep)),
                    snp.index=listKeep,
                    stringsAsFactors=FALSE)

    snp.pruned <- read.gdsn(index.gdsn(gdsSample, "snp.index"))

    listKeepPruned <- which(listKeep %in% snp.pruned)
    snp.pos$pruned[listKeepPruned] <- TRUE

    rm(cnt.total, snp.pruned, listKeepPruned)


    if("normal.geno" %in% ls.gdsn(node=gdsSample)) { # if normal.geno exist mean there is count not in the ref
        # I have other genotype than 1KG
        print("Aye1")
        cnt.total <- read.gdsn(index.gdsn(gdsSample, "Total.count.o"))
        listKeep.o <- which(cnt.total >= minCov)

        snp.pos.o <- data.frame(cnt.tot=cnt.total[listKeep.o],
                        cnt.ref=read.gdsn(index.gdsn(gdsSample,
                                            "Ref.count.o"))[listKeep.o],
                        cnt.alt=read.gdsn(index.gdsn(gdsSample,
                                            "Alt.count.o"))[listKeep.o],
                        snp.pos=read.gdsn(index.gdsn(gds,
                                            "snp.position.o"))[listKeep.o],
                        snp.chr=read.gdsn(index.gdsn(gds,
                                            "snp.chromosome.o"))[listKeep.o],
                        normal.geno=read.gdsn(index.gdsn(gds,
                                            "normal.geno"))[listKeep.o],
                        pruned=rep(0, length(listKeep)),
                        snp.index=rep(0, length(listKeep.o)),
                        stringsAsFactors=FALSE)
        listChr <- unique(snp.pos.o$snp.chr)
        listUnion <- list()

        # if snp.pos.o intersect snp.pos and normal.geno != 3 (we know
        # the genotype of normal) change in snp.pos normal.geno
        z <- cbind(c(snp.pos.o$snp.chr, snp.pos$snp.chr, snp.pos.o$snp.chr),
                    c(snp.pos.o$snp.pos, snp.pos$snp.pos, snp.pos.o$snp.pos),
                    c(seq_len(nrow(snp.pos.o)), 0, -1*seq_len(nrow(snp.pos.o))),
                    c(rep(0, nrow(snp.pos.o)), seq_len(nrow(snp.pos)),
                            rep(0, nrow(snp.pos.o))))
        z <- z[order(z[,1], z[,2], z[,3]), ]
        vCum <- cumsum(z[,3])

        snp.pos[z[ vCum < 0 & z[,3] == 0,4],
                "normal.geno"] <- snp.pos.o[vCum[vCum < 0 & z[,3] == 0],
                                                    "normal.geno"]
        rm(z)

        # Keep the snp.pos.o not in snp.pos
        z <- cbind(c(snp.pos$snp.chr, snp.pos.o$snp.chr, snp.pos$snp.chr),
                    c(snp.pos$snp.pos, snp.pos.o$snp.pos, snp.pos$snp.pos),
                    c(seq_len(nrow(snp.pos)), 0, -1*seq_len(nrow(snp.pos))),
                    c(rep(0, nrow(snp.pos)), seq_len(nrow(snp.pos.o)),
                            rep(0, nrow(snp.pos))))
        z <- z[order(z[,1], z[,2], z[,3]), ]
        snp.pos <- rbind(snp.pos,
                            snp.pos.o[z[cumsum(z[,3] == 0 & z[,3] == 0),4],])

    }
    listCnt <- unique(snp.pos$cnt.tot)
    listCnt <- listCnt[order(listCnt)]

    cutOffA <- data.frame(count = unlist(vapply(listCnt,
                    FUN=function(x, minProb, eProb){
                            return(max(2,qbinom(minProb, x,eProb)))},
                        FUN.VALUE = numeric(1), minProb=minProb, eProb=eProb)),
                    allele = unlist(vapply(listCnt,
                        FUN=function(x, minProb, eProb){
                                return(max(2,qbinom(minProb, x,eProb)))},
                        FUN.VALUE = numeric(1), minProb=minProb, eProb=eProb)))
    row.names(cutOffA) <- as.character(listCnt)

    snp.pos$keep <- rowSums(snp.pos[, c("cnt.ref", "cnt.alt")]) >=
            snp.pos$cnt.tot - cutOffA[as.character(snp.pos$cnt.tot), "count"]

    snp.pos$hetero <- snp.pos$keep == TRUE &
            rowSums(snp.pos[, c("cnt.ref", "cnt.alt")] >=
                        cutOffA[as.character(snp.pos$cnt.tot), "allele"]) == 2

    # We set to homo if 2th allele can be explain by error
    # can switch low allelic fraction to LOH which is less a problem
    # then reduce the allelic ratio by seq error

    snp.pos$homo <- snp.pos$keep == TRUE &
            rowSums(snp.pos[, c("cnt.ref", "cnt.alt")] >=
                        cutOffA[as.character(snp.pos$cnt.tot), "allele"]) == 1

    ## If we know the normal is hetero then we call hetero
    ## if the cnt.alt and cnt.ref > 0
    listHeteroN <- which(snp.pos$homo == TRUE &
                        rowSums(snp.pos[, c("cnt.ref", "cnt.alt")] > 0) == 2 &
                        snp.pos$normal.geno == 1)

    if(length(listHeteroN) > 0) {
        snp.pos$hetero[listHeteroN] <- TRUE
        snp.pos$homo <- FALSE
    }

    return(snp.pos)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param chrInfo a vector chrInfo[i] = length(Hsapiens[[paste0("chr", i)]])
#'         Hsapiens library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' @param snp.pos a \code{data.frame} containing TODO.
#'
#' @param chr a single \code{integer} for the chromosome.
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
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats dbinom
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
computeLOHBlocksDNAChr <- function(gds, chrInfo, snp.pos, chr, genoN=0.0001) {


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

    # snp.pos <- getTableSNV(gds, gdsSample, minCov, minProb, eProb)
    # dfHetero <- snp.pos[snp.pos$hetero == TRUE,]
    # homoBlock <- list()
    # for(chr in unique(snp.pos$snp.chr)){
    #     # Define the end of chr
    #
    #
    # }
    # Cutoff genotype 0.001

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

    z <- z[order(z[,1]),]

    blcSNV <- data.frame(block = cumsum(z[,2])[z[,2] == 0],
                snv = z[z[,2] == 0, 3])
    listAF <- read.gdsn(index.gdsn(gds,"snp.AF"))

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
    for(i in seq_len(nrow(homoBlock))){
        blcCur <- blcSNV[blcSNV$block == i,]
        snvH <- snp.pos[blcCur$snv,]
        lH1 <- 0
        lM1 <- 0
        logLHR <- 0
        homoBlock$nbSNV[i] <- nrow(blcCur)
        homoBlock$nbPruned[i] <- length(which(snvH$pruned))
        if(length(which(snvH$normal.geno != 3)) > 0) {

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

        } else if(length(which(snvH$pruned)) > 2) {

            afSNV <- listAF[snvH$snp.index[which(snvH$pruned)]]
            afSNV <- apply(matrix(afSNV, ncol=1),
                           1,
                           FUN=function(x){max(x, 0.01) })
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

        homoBlock$logLHR[i] <- max(logLHR,-100)
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
#' @param pCutOff TODO
#'
#'
#' @return a list TODO.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats pbinom
#' @encoding UTF-8
#' @export
testEmptyBox <- function(matCov, pCutOff = -3) {

    p <- 0
    pO <- 0
    vMean <- 0.5
    matCov$pWin <- rep(1, nrow(matCov))

    for(i in seq_len(nrow(matCov))){

        vCur1 <- ifelse(matCov$cnt.alt[i] <= matCov$cnt.ref[i],
                            matCov$cnt.alt[i], matCov$cnt.ref[i])
        #vCur2 <- ifelse(matCov$cnt.alt[i] > matCov$cnt.ref[i],
        #                   matCov$cnt.alt[i], matCov$cnt.ref[i])


        pCur <- pbinom(vCur1,
                       size = matCov$cnt.ref[i] + matCov$cnt.alt[i],
                       vMean)
        #print(paste0("pCur ", pCur, " vCur1 ", vCur1, " size ",
        #    matCov$cnt.ref[i] + matCov$cnt.alt[i]))
        pCurO <- max(1 - max(2 * pCur,0.01),0.01)

        matCov$pWin[i] <- pCur * 2
        #print(paste0("Part ",pCur))
        p <- p + log10(max(pCur,0.01))
        pO <- pO + log10(pCurO)
    }
    pCut1 <- as.integer((sum(matCov$pWin < 0.5) >= nrow(matCov)-1) &
                                matCov$pWin[1] < 0.5 &
                                (matCov$pWin[nrow(matCov)] < 0.5) &
                                ((p-pO) <= pCutOff))
    res <- list(pWin = matCov$pWin, p=p,
                    pCut = as.integer(sum(matCov$pWin < 0.5) == nrow(matCov)),
                    pCut1 = pCut1)
    return(res)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param matCov TODO
#'
#' @param pCutOff TODO Default: \code{-3}.
#'
#' @param vMean TODO
#'
#' @return a \code{list} containing:
#' \itemize{
#' \item{pWin} {TODO.}
#' \item{pCut} {TODO.}
#' \item{pCut1} {TODO.}
#' }
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats pbinom
#' @encoding UTF-8
#' @export
testAlleleFractionChange <- function(matCov, pCutOff = -3, vMean){
    p <- 0
    pO <- 0

    matCov$pWin <- rep(1, nrow(matCov))

    for(i in seq_len(nrow(matCov))){

        vCur <- ifelse(matCov$cnt.alt[i] <= matCov$cnt.ref[i],
                        matCov$cnt.alt[i], matCov$cnt.ref[i])

        diff2Mean <- abs(vMean * (matCov$cnt.alt[i] +
                                        matCov$cnt.ref[i]) - vCur)
        pCur1 <- pbinom(round(vMean * (matCov$cnt.alt[i] +
                                            matCov$cnt.ref[i]) - diff2Mean),
                    size = matCov$cnt.ref[i] + matCov$cnt.alt[i], vMean)
        pCur2 <- 1 - pbinom(round(vMean * (matCov$cnt.alt[i] +
                                matCov$cnt.ref[i]) + diff2Mean),
                        size = matCov$cnt.ref[i] + matCov$cnt.alt[i], vMean)

        pCur <- pCur1 + pCur2

        matCov$pWin[i] <- pCur


        pCurO <- max(1 - max(pCur,0.01),0.01)

        p <- p + log10(max(pCur,0.01))
        pO <- pO + log10(pCurO)
    }
    pCut1 <- as.integer((sum(matCov$pWin < 0.5) >= nrow(matCov)-1) &
                                matCov$pWin[1] < 0.5 &
                                (matCov$pWin[nrow(matCov)] < 0.5)  &
                                ((p-pO) <= pCutOff))
    res <- list(pWin = matCov$pWin, p=p,
                    pCut = as.integer(sum(matCov$pWin < 0.5) == nrow(matCov)),
                    pCut1 = pCut1)
    return(res)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param snp.pos TODO
#'
#' @param chr a single positive \code{integer} for the chromosome.
#'
#' @param wAR a single positive \code{integer} representing the size-1 of
#' the window used to compute an empty box. Default: \code{10}.
#'
#' @param cutOffEmptyBox TODO Default: \code{-3}.
#'
#' @return a \code{vector} of TODO representing the imbAR TODO.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
computeAllelicImbDNAChr <- function(snp.pos, chr, wAR=10,
                                        cutOffEmptyBox=-3) {

    ## The chr parameter must be a single integer value
    if (!isSingleNumber(chr))  {
        stop("The \'chr\' must be a single integer value representing ",
             "a chromosome")
    }

    ## The wAR parameter must be a single positive numeric superior to 1
    if (!(isSingleNumber(wAR) && (wAR >= 1)))  {
        stop("The \'wAR\' must be a single numeric positive value.")
    }

    # We use wAR - 1 because
    # process the window ex: 1 to 1+wAR
    wAR <- wAR - 1
    listHetero <- NULL
    if(length(which(snp.pos$normal.geno != 3) > 0)){
        listHetero <- which(snp.pos$keep == TRUE & snp.pos$normal.geno == 1)
    } else{
        listHetero <- which(snp.pos$hetero == TRUE)
    }

    heteroSNV <- snp.pos[listHetero,]

    if(nrow(heteroSNV) > wAR) {
        for(i in seq_len(nrow(heteroSNV)-wAR)) {
            if(sum(snp.pos[listHetero[i]:listHetero[(i+wAR-1)], "LOH"]) == 0 ){
                cur <- testEmptyBox(heteroSNV[i:(i+wAR), c
                                    ("cnt.alt", "cnt.ref")], cutOffEmptyBox)
                if(cur$pCut == 1){
                    # Set all snv from tmpA (include homozygotes)
                    # in the window  to 1
                    snp.pos[listHetero[i]:listHetero[(i+wAR)], "imbAR"] <- 1
                }
            }
        }
    }
    snp.pos$imbAR[which(snp.pos$LOH == 1)] <- 0
    return(snp.pos$imbAR)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param snp.pos TODO
#'
#' @param chr A integer for the chromosome TODO
#'
#' @param w a single positive \code{numeric} representing the size of the
#' window to compute the allelic fraction.
#' Default: \code{10}.
#'
#' @param cutOff TODO . Default: \code{-3}.
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats median
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
computeAlleleFraction <- function(snp.pos, chr, w=10, cutOff=-3) {

    ## The w parameter must be a single positive numeric superior to 1
    if (!(isSingleNumber(w) && (w >= 1)))  {
        stop("The \'w\' must be a single numeric positive value.")
    }

    listBlockAR <- list()
    j <- 1
    tmp <- as.integer(snp.pos$imbAR == 1)
    z <- cbind(c(tmp[1], tmp[-1] - tmp[seq_len(length(tmp) -1)]),
            c(tmp[-1] - tmp[seq_len(length(tmp) -1)], tmp[length(tmp)] * -1))

    if(length(which(z[,1] == 1)) > 0) {

        segImb <- data.frame(start = seq_len(nrow(snp.pos))[which(z[,1] > 0)],
                                end = seq_len(nrow(snp.pos))[which(z[,2] < 0)])

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


                        if(nrow(snp.hetero) - start < w){ # Close the segment
                            lapCur <- median(apply(snp.hetero[start:nrow(snp.hetero), c("cnt.ref", "cnt.alt")], 1, min) / (rowSums(snp.hetero[start:nrow(snp.hetero),c("cnt.ref", "cnt.alt")])))


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
                            lapCur <- median(apply(snp.hetero[start:nrow(snp.hetero), c("cnt.ref", "cnt.alt")], 1, min) / (rowSums(snp.hetero[start:nrow(snp.hetero),c("cnt.ref", "cnt.alt")])))

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
                lapCur <- median(apply(snp.hetero[, c("cnt.ref", "cnt.alt")], 1, min) / (rowSums(snp.hetero[,c("cnt.ref", "cnt.alt")])))

                listBlockAR[[j]] <- c(segImb$start[i], segImb$end[i], lapCur)

                j <- j + 1
            }
        }
    }

    # note NULL if length(listBlockAR) == 0
    listBlockAR <- do.call(rbind, listBlockAR)

    return(listBlockAR)
}


#' @title TODO
#'
#' @description Create a data.frame with the allelic fraction TODO
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param gdsSample TODO
#'
#' @param sampleCurrent A \code{string} corresponding to
#' the sample.id
#' use in LDpruning
#'
#' @param study.id A \code{character} string corresponding to the study
#' used in LDpruning
#'
#' @param chrInfo a vector chrInfo[i] = length(Hsapiens[[paste0("chr", i)]])
#'         Hsapiens library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' @param minCov an \code{integer} default 10
#'
#' @param minProb a single \code{numeric} between 0 and 1 representing TODO.
#' Default: \code{0.999}.
#'
#' @param eProb a single \code{numeric} between 0 and 1 representing the
#' probability of sequencing error. Default: \code{0.001}.
#'
#' @param cutOffLOH log of the score to be LOH . Default: \code{-5}.
#'
#' @param cutOffHomoScore TODO. Default: \code{-3}.
#'
#' @param wAR a single positive \code{integer} representing the size-1 of
#' the window used to compute an empty box. Default: \code{9}.
#'
#' @return a \code{data.frame} with lap for the snv with depth > minCov. TODO
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
computeAllelicFractionDNA <- function(gds, gdsSample,
                                        sampleCurrent, study.id,
                                        chrInfo, minCov=10,
                                        minProb=0.999,
                                        eProb=0.001,
                                        cutOffLOH=-5,
                                        cutOffHomoScore=-3,
                                        wAR=9) {

    ## The minProb parameter must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(minProb) && (minProb >= 0.0) && (minProb <= 1.0))) {
        stop("The \'minProb\' must be a single numeric positive ",
                "value between 0 and 1.")
    }

    ## The eProb parameter must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(eProb) && (eProb >= 0.0) && (eProb <= 1.0))) {
        stop("The \'eProb\' must be a single numeric positive ",
                "value between 0 and 1.")
    }

    ## The wAR parameter must be a single positive numeric superior to 1
    if (!(isSingleNumber(wAR) && (wAR >= 1)))  {
        stop("The \'wAR\' must be a single numeric positive value.")
    }

    snp.pos <- getTableSNV(gds, gdsSample, sampleCurrent, study.id,
                            minCov, minProb, eProb)
    snp.pos$lap <- rep(-1, nrow(snp.pos))
    snp.pos$LOH <- rep(0, nrow(snp.pos))
    snp.pos$imbAR <- rep(-1, nrow(snp.pos))

    homoBlock <- list()

    for(chr in unique(snp.pos$snp.chr)) {
        # Define the end of chr
        print(paste0("chr ", chr))
        print(paste0("Step 1 ", Sys.time()))


        #listHetero <- dfHetero[dfHetero$snp.chr == chr, "snp.pos"]
        listChr <- which(snp.pos$snp.chr == chr)
        # snp.pos.chr <- snp.pos[listChr,]


        homoBlock[[chr]] <- computeLOHBlocksDNAChr(gds=gds, chrInfo=chrInfo,
                                            snp.pos=snp.pos[listChr,], chr=chr)

        print(paste0("Step 2 ", Sys.time()))

        homoBlock[[chr]]$LOH <- as.integer(homoBlock[[chr]]$logLHR <=
                cutOffLOH & homoBlock[[chr]]$homoScore <= cutOffHomoScore)

        z <- cbind(c(homoBlock[[chr]]$start, homoBlock[[chr]]$end,
                        snp.pos[listChr, "snp.pos"]),
                    c(rep(0,  2* nrow(homoBlock[[chr]])),
                        rep(1, length(listChr))),
                    c(homoBlock[[chr]]$LOH,
                        -1 * homoBlock[[chr]]$LOH,
                        rep(0, length(listChr)) ),
                    c(rep(0, 2 * nrow(homoBlock[[chr]])),
                        seq_len(length(listChr))))

        z <- z[order(z[,1], z[,2]), ]
        pos <- z[cumsum(z[,3]) > 0 & z[,4] > 0, 4]
        snp.pos[listChr[pos], "lap"] <- 0
        snp.pos[listChr[pos], "LOH"] <- 1

        print(paste0("Step 3 ", Sys.time()))

        snp.pos[listChr, "imbAR"] <-
                computeAllelicImbDNAChr(snp.pos=snp.pos[listChr, ], chr=chr,
                                            wAR=10, cutOffEmptyBox=-3)

        print(paste0("Step 4 ", Sys.time()))
        blockAF <- computeAlleleFraction(snp.pos=snp.pos[listChr, ], chr=chr,
                                                    w=10, cutOff=-3)

        print(paste0("Step 5 ", Sys.time()))

        if(! is.null(blockAF)) {
            for(i in seq_len(nrow(blockAF))) {
                snp.pos[listChr[blockAF[i, 1]:blockAF[i, 2]], "lap"] <-
                                blockAF[i, 3]
            }
        }
    }

    snp.pos[which(snp.pos[, "lap"] == -1), "lap"] <- 0.5

    return(snp.pos)
}



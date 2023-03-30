#' @title Validate input parameters for getTableSNV() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{getTableSNV}} function.
#'
#' @param gdsRef an object of class \code{\link[gdsfmt]{gds.class}} (a GDS
#' file), the opened 1KG GDS file (reference).
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Sample GDS file.
#'
#' @param sampleCurrent a \code{character} string corresponding to
#' the sample identifier used in \code{\link{pruningSample}} function.
#'
#' @param study.id a \code{character} string corresponding to the study
#' identifier used in \code{\link{pruningSample}} function.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage.
#'
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#'
#' @param eProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error.
#'
#' @param verbose a \code{logicial} indicating if messages should be printed
#' when the function is running.
#'
#' @return The integer \code{0L} when the function is successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located (refence and sample files)
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened) (reference file)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The Sample GDS file (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                          "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validateGetTableSNV(gdsRef=gds1KG, gdsSample=gdsSample,
#'      sampleCurrent="A101TCGA", study.id="TCGA", minCov=10L,
#'      minProb=0.998, eProb=0.002, verbose=TRUE)
#'
#' ## All GDS files must be closed
#' closefn.gds(gdsfile=gds1KG)
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateGetTableSNV <- function(gdsRef, gdsSample, sampleCurrent, study.id,
                                    minCov, minProb, eProb, verbose) {

    ## The gdsRef must be an object of class "gds.class"
    validateGDSClass(gds=gdsRef, name="gds")

    ## The gdsSample must be an object of class "gds.class"
    validateGDSClass(gds=gdsSample, name="gdsSample")

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

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, "verbose")

    ## Successful
    return(0L)
}


#' @title Estimate the allelic fraction of the pruned SNVs for a specific
#' DNA-seq sample
#'
#' @description The function creates a \code{data.frame} containing the
#' allelic fraction for the pruned SNV dataset specific to a DNA-seq sample.
#'
#' @param gds an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the GDS Sample file.
#'
#' @param sampleCurrent a \code{character} string corresponding to
#' the sample identifier as used in \code{\link{pruningSample}} function.
#'
#' @param study.id a \code{character} string corresponding to the name of
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
#' message when running. Default: \code{FALSE}.
#'
#' @return a \code{data.frame} with lap for the pruned SNV dataset with
#' coverage > \code{minCov}. TODO
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
#' @keywords internal
computeAllelicFractionDNA <- function(gds, gdsSample, sampleCurrent, study.id,
                                chrInfo, minCov=10L, minProb=0.999,
                                eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                                wAR=9L, verbose=FALSE) {

    ## Extract the genotype information for a SNV dataset using
    ## the GDS Sample file and the 1KG GDS file
    snp.pos <- getTableSNV(gds, gdsSample, sampleCurrent, study.id,
                            minCov, minProb, eProb)

    snp.pos$lap <- rep(-1, nrow(snp.pos))
    snp.pos$LOH <- rep(0, nrow(snp.pos))
    snp.pos$imbAR <- rep(-1, nrow(snp.pos))

    homoBlock <- list()

    for(chr in unique(snp.pos$snp.chr)) {

        if (verbose) {
            message("chr ", chr)
            message("Step 1 ", Sys.time())
        }

        listChr <- which(snp.pos$snp.chr == chr)

        homoBlock[[chr]] <- computeLOHBlocksDNAChr(gds=gds, chrInfo=chrInfo,
                                        snp.pos=snp.pos[listChr,], chr=chr)

        if (verbose) { message("Step 2 ", Sys.time()) }

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

        if (verbose) { message("Step 3 ", Sys.time()) }

        snp.pos[listChr, "imbAR"] <-
            computeAllelicImbDNAChr(snp.pos=snp.pos[listChr, ], chr=chr,
                                        wAR=10, cutOffEmptyBox=-3)

        if (verbose) { message("Step 4 ", Sys.time()) }

        blockAF <- computeAlleleFraction(snp.pos=snp.pos[listChr, ], chr=chr,
                                            w=10, cutOff=-3)

        if (verbose) { message("Step 5 ", Sys.time()) }

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


#' @title Estimate the allelic fraction of the pruned SNVs for a specific
#' RNA-seq sample
#'
#' @description The function creates a \code{data.frame} containing the
#' allelic fraction for the pruned SNV dataset specific to a RNA-seq sample.
#' TODO
#'
#' @param gds an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the GDS Sample file.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file.
#'
#' @param sampleCurrent a \code{character} string corresponding to
#' the sample identifier as used in \code{\link{pruningSample}} function.
#'
#' @param study.id a \code{character} string corresponding to the name of
#' the study as
#' used in \code{\link{pruningSample}} function.
#'
#' @param block.id a \code{character} corresponding to the field gene block
#' in the GDS \code{gdsRefAnnot} to use split by gene.
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
#' @param cutOffLOH a single log of the score to be LOH TODO.
#' Default: \code{-5}.
#'
#' @param cutOffAR a single \code{numeric} representing the cutoff, in
#' log score, to tag SNVs located in a gene has having an allelic fraction
#' different 0.5
#' Default: \code{3}.
#'
#' @param verbose a \code{logicial} indicating if the function should print
#' message when running. Default: \code{FALSE}.
#'
#' @return a \code{data.frame} with lap for the pruned SNV dataset with
#' coverage > \code{minCov}. TODO
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn ls.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computeAllelicFractionRNA <- function(gds, gdsSample, gdsRefAnnot,
                    sampleCurrent, study.id, block.id, chrInfo, minCov=10L,
                    minProb=0.999, eProb=0.001, cutOffLOH=-5,
                    cutOffAR=3, verbose=FALSE) {

    ## Extract the genotype information for a SNV dataset using
    ## the GDS Sample file and the 1KG GDS file
    snp.pos <- getTableSNV(gds, gdsSample, sampleCurrent, study.id,
                                minCov, minProb, eProb)
    # Keep only SNV in GDS ref because to reduce SNV artefact from RNA
    snp.pos <- snp.pos[which(snp.pos$snp.index > 0),]

    # Get the block structure base on genes from gdsRefAnnot
    snp.pos$block.id <- get.Gene.Block(gdsRefAnnot, snp.pos$snp.index,
                                            block.id)

    snp.pos$phase <- rep(3, nrow(snp.pos))
    if ("phase" %in% ls.gdsn(node=gdsSample)) {
        snp.pos$phase <- read.gdsn(index.gdsn(gdsSample,
                                            "phase"))[snp.pos$snp.index]
    }
    snp.pos$lap <- rep(-1, nrow(snp.pos))
    snp.pos$LOH <- rep(0, nrow(snp.pos))
    snp.pos$imbAR <- rep(-1, nrow(snp.pos))
    snp.pos$freq <- read.gdsn(index.gdsn(gds, "snp.AF"))[snp.pos$snp.index]
    # for each chromosome
    listBlock <- list()
    for(chr in unique(snp.pos$snp.chr)) {

        if (verbose) {
            message("chr ", chr)
            message("Step 1 ", Sys.time())
        }

        #listHetero <- dfHetero[dfHetero$snp.chr == chr, "snp.pos"]
        listChr <- which(snp.pos$snp.chr == chr)
        # snp.pos.chr <- snp.pos[listChr,]

        blockAF <- tableBlockAF(snp.pos=snp.pos[listChr,])
        # LOH
        blockAF$aRF[blockAF$lRhomo <= cutOffLOH] <- 0
        blockAF$aRF[blockAF$lR >= cutOffAR] <- blockAF$aFraction[blockAF$lR
                                                            >= cutOffAR]
        blockAF$aRF[blockAF$lR < cutOffAR & blockAF$nbHetero > 1] <- 0.5

        listBlock[[chr]] <- blockAF

        if (verbose) {
            message("Step 1 done ", Sys.time())
        }
    }

    blockAF <- do.call(rbind, listBlock)
    listMissing <- which(abs(blockAF$aRF + 1) < 1e-6)
    blockAF[listMissing, "aRF"] <- sample(blockAF$aRF[-1*listMissing],
                                        length(listMissing), replace=TRUE)

    for(b in seq_len(nrow(blockAF))) {
        snp.pos$lap[snp.pos$block.id == blockAF$block[b]] <- blockAF$aRF[b]
    }

    return(snp.pos)
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
#' @keywords internal
computeAllelicImbDNAChr <- function(snp.pos, chr, wAR=10, cutOffEmptyBox=-3) {

    # We use wAR - 1 because
    # process the window ex: 1 to 1+wAR
    wAR <- wAR - 1
    listHetero <- NULL
    if(length(which(snp.pos$normal.geno != 3) > 0)) {
        listHetero <- which(snp.pos$keep == TRUE & snp.pos$normal.geno == 1)
    } else{
        listHetero <- which(snp.pos$hetero == TRUE)
    }

    heteroSNV <- snp.pos[listHetero, ]

    if(nrow(heteroSNV) > wAR) {
        for(i in seq_len(nrow(heteroSNV)-wAR)) {
            if(sum(snp.pos[listHetero[i]:listHetero[(i+wAR-1)], "LOH"]) == 0 ) {
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
#' @param matCov TODO
#'
#' @param pCutOff TODO Default: \code{-3}.
#'
#' @param vMean TODO
#'
#' @return a \code{list} containing 3 entries:
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
#' @keywords internal
testAlleleFractionChange <- function(matCov, pCutOff=-3, vMean) {
    p <- 0
    pO <- 0

    matCov$pWin <- rep(1, nrow(matCov))

    for(i in seq_len(nrow(matCov))) {

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

        pCurO <- max(1 - max(pCur, 0.01), 0.01)

        p <- p + log10(max(pCur, 0.01))
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


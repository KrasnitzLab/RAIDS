#' @title Extract the genotype information for a SNV dataset using
#' the Profile GDS file and the 1KG GDS file
#'
#' @description The function generates a \code{data.frame} containing the
#' genotype information from a initial list of SNVs associated to a specific
#' profile. The function uses the information present in the 1KG GDS file
#' (reference) and the Profile GDS file.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}} (a
#' GDS file), the opened 1KG GDS file.
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
#' \item{snp.pos} {a single \code{integer} representing the SNV position.}
#' \item{snp.chr} {a single \code{integer} representing the SNV chromosome.}
#' \item{normal.geno} {a single \code{numeric} indicating the genotype of the
#' SNV. The possibles are: \code{0} (wild-type homozygote), \code{1}
#' (heterozygote), \code{2} (altenative homozygote), \code{3} indicating that
#' the normal genotype is unknown.}
#' \item{snp.index} {a \code{vector} of \code{integer} representing the
#' position of the SNVs in the 1KG GDS file.}
#' }
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
#'
#' ## Copy the Profile GDS file demo that has been pruned and annotated
#' ## into a test directory (deleted after the example has been run)
#' dataDirAllelicFraction <- file.path(system.file("extdata", package="RAIDS"),
#'                  "demoAllelicFraction")
#' dir.create(dataDirAllelicFraction, showWarnings=FALSE,
#'                  recursive=FALSE, mode="0777")
#' file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                  file.path(dataDirAllelicFraction, "ex1.gds"))
#'
#' ## Open the reference GDS file (demo version)
#' gds1KG <- snpgdsOpen(fileGDS)
#'
#' ## Profile GDS file for one profile
#' fileProfile <- file.path(dataDirAllelicFraction, "ex1.gds")
#' profileGDS <- openfn.gds(fileProfile)
#'
#' ## The function returns a data frame containing the SNVs information
#' result <- RAIDS:::getTableSNV(gdsReference=gds1KG, gdsSample=profileGDS,
#'     currentProfile="ex1", studyID="MYDATA", minCov=10L, minProb=0.999,
#'     eProb=0.001, verbose=FALSE)
#' head(result)
#'
#' ## Close both GDS files (important)
#' closefn.gds(profileGDS)
#' closefn.gds(gds1KG)
#'
#' ## Unlink Profile GDS file (created for demo purpose)
#' unlink(file.path(dataDirAllelicFraction, "ex1.gds"))
#' unlink(dataDirAllelicFraction)
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

    ## Extract SNV coverage from Sample GDS file
    cnt.total <- read.gdsn(node=index.gdsn(gdsSample, "Total.count"),
                            start=c(1, posCur), count=c(-1, 1))

    ## Only retained the SNV with the minimum required coverage
    listKeep <- cnt.total@i[which(cnt.total@x >= minCov)] + 1

    ## Create the data.frame with the required information
    snp.pos <- data.frame(cnt.tot=cnt.total[listKeep],
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
    snp.pos$pruned[listKeepPruned] <- TRUE

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

        snp.pos.o <- data.frame(cnt.tot=cnt.total[listKeep.o],
            cnt.ref=read.gdsn(index.gdsn(gdsSample, "Ref.count.o"))[listKeep.o],
            cnt.alt=read.gdsn(index.gdsn(gdsSample, "Alt.count.o"))[listKeep.o],
            snp.pos=read.gdsn(index.gdsn(gdsReference,
                "snp.position.o"))[listKeep.o],
            snp.chr=read.gdsn(index.gdsn(gdsReference,
                "snp.chromosome.o"))[listKeep.o],
            normal.geno=read.gdsn(index.gdsn(gdsSample,
                "normal.geno"))[listKeep.o], pruned=rep(0, length(listKeep)),
            snp.index=rep(0, length(listKeep.o)), stringsAsFactors=FALSE)
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

        snp.pos[z[ vCum < 0 & z[,3] == 0, 4], "normal.geno"] <-
            snp.pos.o[vCum[vCum < 0 & z[, 3] == 0], "normal.geno"]
        rm(z)

        # Keep the snp.pos.o not in snp.pos
        z <- cbind(c(snp.pos$snp.chr, snp.pos.o$snp.chr, snp.pos$snp.chr),
                c(snp.pos$snp.pos, snp.pos.o$snp.pos, snp.pos$snp.pos),
                c(seq_len(nrow(snp.pos)), 0, -1*seq_len(nrow(snp.pos))),
                c(rep(0, nrow(snp.pos)), seq_len(nrow(snp.pos.o)),
                                rep(0, nrow(snp.pos))))
        z <- z[order(z[,1], z[,2], z[,3]), ]
        # merge snp.pos with snp.pos.o
        snp.pos <- rbind(snp.pos,
                            snp.pos.o[z[cumsum(z[,3] == 0 & z[,3] == 0),4],])

    }
    listCnt <- unique(snp.pos$cnt.tot)
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

    if (length(listHeteroN) > 0) {
        snp.pos$hetero[listHeteroN] <- TRUE
        snp.pos$homo <- FALSE
    }

    return(snp.pos)
}


#' @title Estimate the allelic fraction of the pruned SNVs for a specific
#' DNA-seq profile
#'
#' @description The function creates a \code{data.frame} containing the
#' allelic fraction for the pruned SNV dataset specific to a DNA-seq profile
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened 1KG GDS file.
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
#' SNV in the 1KG GDS file that contains all SNVs}
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
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
#'
#' ## Copy the Profile GDS file demo that has been pruned and annotated
#' ## into a test directory (deleted after the example has been run)
#' dataDirAllelicFraction <- file.path(system.file("extdata", package="RAIDS"),
#'                  "demoAllelicFraction")
#' dir.create(dataDirAllelicFraction, showWarnings=FALSE,
#'                  recursive=FALSE, mode="0777")
#' file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                  file.path(dataDirAllelicFraction, "ex1.gds"))
#'
#' ## Open the reference GDS file (demo version)
#' gds1KG <- snpgdsOpen(fileGDS)
#'
#' ## Profile GDS file for one profile
#' fileProfile <- file.path(dataDirAllelicFraction, "ex1.gds")
#' profileGDS <- openfn.gds(fileProfile)
#'
#' ## Chromosome length information
#' ## chr23 is chrX, chr24 is chrY and chrM is 25
#' chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'     181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
#'     135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
#'     83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
#'     156040895L, 57227415L,  16569L)
#'
#' ## The function returns a data frame containing the allelic fraction info
#' result <- RAIDS:::computeAllelicFractionDNA(gdsReference=gds1KG,
#'     gdsSample=profileGDS,
#'     currentProfile="ex1", studyID="MYDATA", chrInfo=chrInfo, minCov=10L,
#'     minProb=0.999, eProb=0.001, cutOffLOH=-5,
#'     cutOffHomoScore=-3, wAR=9L, verbose=FALSE)
#' head(result)
#'
#' ## Close both GDS files (important)
#' closefn.gds(profileGDS)
#' closefn.gds(gds1KG)
#'
#' ## Unlink Profile GDS file (created for demo purpose)
#' unlink(file.path(dataDirAllelicFraction, "ex1.gds"))
#' unlink(dataDirAllelicFraction)
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
    ## the Profile GDS file and the 1KG GDS file
    snp.pos <- getTableSNV(gdsReference=gdsReference, gdsSample=gdsSample,
        currentProfile=currentProfile, studyID=studyID, minCov=minCov,
        minProb=minProb, eProb=eProb, verbose=verbose)

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

        homoBlock[[chr]] <- computeLOHBlocksDNAChr(gdsReference=gdsReference,
                chrInfo=chrInfo, snp.pos=snp.pos[listChr,], chr=chr)

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
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the GDS Sample file.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file.
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
#' @param cutOffLOH a single log of the score to be LOH TODO.
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
#' @return a \code{data.frame} with lap for the pruned SNV dataset with
#' coverage > \code{minCov}. TODO
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
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
    ## the GDS Sample file and the 1KG GDS file
    snp.pos <- getTableSNV(gdsReference, gdsSample, currentProfile, studyID,
                                minCov, minProb, eProb)
    # Keep only SNV in GDS ref because to reduce SNV artefact from RNA
    snp.pos <- snp.pos[which(snp.pos$snp.index > 0),]

    # Get the block structure base on genes from gdsRefAnnot
    snp.pos$block.id <- getGeneBlock(gdsRefAnnot=gdsRefAnnot,
                            snp.index=snp.pos$snp.index, blockID=blockID)

    snp.pos$phase <- rep(3, nrow(snp.pos))
    if ("phase" %in% ls.gdsn(node=gdsSample)) {
        snp.pos$phase <- read.gdsn(index.gdsn(gdsSample,
                                            "phase"))[snp.pos$snp.index]
    }
    snp.pos$lap <- rep(-1, nrow(snp.pos))
    snp.pos$LOH <- rep(0, nrow(snp.pos))
    snp.pos$imbAR <- rep(-1, nrow(snp.pos))
    snp.pos$freq <- read.gdsn(index.gdsn(gdsReference,
                                            "snp.AF"))[snp.pos$snp.index]
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
#' dataDir <- system.file("extdata", package="RAIDS")
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
#' dataDir <- system.file("extdata", package="RAIDS")
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
                    size=matCov$cnt.ref[i] + matCov$cnt.alt[i], vMean)
        pCur2 <- 1 - pbinom(round(vMean * (matCov$cnt.alt[i] +
                    matCov$cnt.ref[i]) + diff2Mean),
                    size=matCov$cnt.ref[i] + matCov$cnt.alt[i], vMean)

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
    res <- list(pWin=matCov$pWin, p=p,
                pCut=as.integer(sum(matCov$pWin < 0.5) == nrow(matCov)),
                pCut1=pCut1)
    return(res)
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
                matCov$pWin[1] < 0.5 & (matCov$pWin[nrow(matCov)] < 0.5) &
                            ((p-pO) <= pCutOff))
    res <- list(pWin=matCov$pWin, p=p,
                pCut=as.integer(sum(matCov$pWin < 0.5) == nrow(matCov)),
                pCut1=pCut1)
    return(res)
}


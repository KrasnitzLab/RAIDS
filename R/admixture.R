

#' @title create a file tfam file for plink from the gdsReference file
#'
#' @description TODO
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param listSample a \code{array} with the sample to keep TODO
#'
#' @param pedOut TODO a PATH and file name to the output file
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils write.table
#' @encoding UTF-8
#' @keywords export
gds2tfam <- function(gdsReference, listSample, pedOut) {

    sampleGDS <- index.gdsn(gdsReference, "sample.id")
    sampleId <-read.gdsn(sampleGDS)
    listS <- which(sampleId %in% listSample)

    sampleGDS <- index.gdsn(gdsReference, "sample.annot")
    sampleANNO <-read.gdsn(sampleGDS)

    dfPed <- data.frame(famId=paste0("F", seq_len(length(listSample))),
                            id=sampleId[listS],
                            fa=rep("0",length(listSample)),
                            mo=rep("0",length(listSample)),
                            sex=sampleANNO$sex[listS],
                            pheno=rep(1,length(listSample)),
                            stringsAsFactors=FALSE)

    write.table(dfPed, pedOut,
                quote=FALSE, sep="\t",
                row.names=FALSE,
                col.names=FALSE)

}


#' @title create a file tfam file for plink from the gdsProfile file
#'
#' @description TODO
#'
#' @param gdsProfile an object of class \link[gdsfmt]{gds.class} (a GDS
#' file), the open Profile GDS file.
#'
#' @param listSample  a \code{array} with the sample to keep
#'
#' @param sampleANNO a \code{data.frame} with at least column sex and the name
#' must be sample.id
#'
#' @param pedOut TODO a PATH and file name to the output file
#'
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils write.table
#' @encoding UTF-8
#' @keywords export
gds2tfamSample <- function(gdsProfile, listSample, sampleANNO, pedOut) {

    sampleGDS <- index.gdsn(gdsProfile, "sample.id")
    sampleId <-read.gdsn(sampleGDS)
    listS <- which(sampleId %in% listSample)

    dfPed <- data.frame(famId=paste0("F", seq_len(length(listSample))),
                            id=sampleId[listS],
                            fa=rep("0",length(listSample)),
                            mo=rep("0",length(listSample)),
                            sex=sampleANNO[sampleId[listS], "sex"],
                            pheno=rep(1,length(listSample)),
                            stringsAsFactors=FALSE)

    write.table(dfPed, pedOut, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE)
}


#' @title create a file tped file for plink from the gds file
#'
#' @description TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param listSample  a \code{array} with the sample to keep
#'
#' @param listSNP  a \code{array} with the snp.id to keep
#'
#' @param pedOut TODO a PATH and file name to the output file
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils write.table
#' @encoding UTF-8
#' @keywords export
gds2tped <- function(gds, listSample, listSNP, pedOut) {

    sampleGDS <- index.gdsn(gds, "sample.id")
    sampleId <-read.gdsn(sampleGDS)
    listS <- which(sampleId %in% listSample)

    snpGDS <- index.gdsn(gds, "snp.id")
    snpId <- read.gdsn(snpGDS)
    listKeep <- which(snpId %in% listSNP)
    snpId <- snpId[listKeep]

    snpGDS <- index.gdsn(gds, "snp.chromosome")
    snpChr <- read.gdsn(snpGDS)
    snpChr <- snpChr[listKeep]

    snpGDS <- index.gdsn(gds, "snp.position")
    snpPos <- read.gdsn(snpGDS)
    snpPos <- snpPos[listKeep]

    tped <- list()
    tped[[1]] <- snpChr
    tped[[2]] <- snpId
    tped[[3]] <- rep(0,length(snpId))
    tped[[4]] <- snpPos
    k<-4
    geno.var <- index.gdsn(gds, "genotype")
    for(i in listS){
        k <- k + 1

        tmp <- read.gdsn(geno.var, start=c(1, i), count=c(-1,1))[listKeep]

        # 0 1 0 1 0 1
        tped[[k]] <- (tmp == 2) + 1
        k <- k + 1
        tped[[k]] <- (tmp > 0) + 1

    }

    write.table(tped, pedOut, quote=FALSE, sep="\t", row.names=FALSE,
                    col.names=FALSE)
}


#' @title Merge the pruning files by chromosome in one RDS file
#'
#' @description The function reads the information from all chromosomes. The
#' information is scattered in different files (one file per chromosome).
#' Once all information is loaded,  the function merges
#' the information and saves it into a RDS file.
#'
#' @param pathPrunedGDS a \code{character} string representing the path where
#' the pruned files for each chromosome are located.
#' The path must exists.
#'
#' @param filePref a \code{character} string representing the prefix used for
#' the pruned files for each chromosome. The prefix represent the complete
#' string that is before the chromosome number in the file names.
#'
#' @param fileOut a \code{character} string representing name of the output
#' file that will be created. The file will contain the information for all
#' pruned chromosomes. The file must have a ".rds" extension.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
groupChrPruning <- function(pathPrunedGDS, filePref, fileOut) {

    prunedList <- list()

    # Read the content of each file (one file per chromosome)
    for(i in seq_len(22)) {
        fileChr <- file.path(pathPrunedGDS, paste0(filePref, i, ".rds"))

        if(file.exists(fileChr)) {
            prunedList[[i]] <- readRDS(fileChr)
        } else {
            stop("Problem with the file: ", fileChr)
        }
    }

    ## Merge the content of all files
    pruned <- do.call(c, prunedList)

    ## Save all the information into one file
    saveRDS(pruned, fileChr <- file.path(pathPrunedGDS, fileOut))

    ## Successful
    return(0L)
}

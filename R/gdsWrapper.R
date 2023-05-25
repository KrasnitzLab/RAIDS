#' @title Append sample names into a GDS file
#'
#' @description This function append the sample identifiers into the
#' "samples.id" node of a GDS file.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param listSample a \code{vector} of \code{character} string representing
#' the sample identifiers to be added to GDS file.
#'
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(dataDir, "GDS_TEMP_04.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create "sample.id" node (the node must be present)
#' add.gdsn(node=GDS_file_tmp, name="sample.id", val=c("sample_01",
#'     "sample_02"))
#'
#' sync.gds(gdsfile=GDS_file_tmp)
#'
#' ## Add information about 2 samples to the GDS file
#' RAIDS:::appendGDSSampleOnly(gds=GDS_file_tmp,
#'     listSamples=c("sample_03", "sample_04"))
#'
#' ## Read sample identifier list
#' ## Only "sample_03" and "sample_04" should have been added
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.id"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn append.gdsn
#' @encoding UTF-8
#' @keywords internal
appendGDSSampleOnly <- function(gds, listSamples) {

    sampleGDS <- index.gdsn(gds, "sample.id")

    append.gdsn(sampleGDS, val=listSamples, check=TRUE)

    return(0L)
}


#' @title Add information related to SNVs into a Reference GDS file
#'
#' @description the function adds the SNV information into a Reference
#' GDS file.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Reference GDS file.
#'
#' @param fileFREQ a \code{character} string representing the path and file
#' name of the RDS file with the filtered SNP information.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required package
#' library(withr)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Temporary Reference GDS file
#' file1KG <- local_file(file.path(dataDir, "1KG_TEMP_002.gds"))
#' filenewGDS <- createfn.gds(file1KG)
#'
#' ## The RDS file containing the filtered SNP information
#' fileFilerterSNVs <- file.path(dataDir, "mapSNVSelected_Demo.rds")
#'
#' ## Add SNV information to Reference GDS
#' RAIDS:::generateGDSSNPinfo(gdsReference=filenewGDS,
#'     fileFREQ=fileFilerterSNVs, verbose=TRUE)
#'
#' ## Close GDS file (important)
#' closefn.gds(filenewGDS)
#'
#' ## Remove temporary files
#' deferred_run()
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
generateGDSSNPinfo <- function(gdsReference, fileFREQ, verbose) {

    mapSNVSel <- readRDS(file=fileFREQ)

    if(verbose) { message("Read mapSNVSel DONE ", Sys.time()) }

    add.gdsn(node=gdsReference, name="snp.id",
                    paste0("s", seq_len(nrow(mapSNVSel))))

    if(verbose) { message("SNP part snp.id DONE ", Sys.time()) }

    add.gdsn(node=gdsReference, name="snp.chromosome",
            as.integer(gsub("chr", "", mapSNVSel$CHROM)), storage = "uint16")

    if(verbose) { message("SNP part snp.chromosome DONE ", Sys.time()) }

    add.gdsn(node=gdsReference, name="snp.position", as.integer(mapSNVSel$POS),
                storage="int32")

    if(verbose) { message("SNP part snp.position DONE ", Sys.time()) }

    add.gdsn(node=gdsReference, name="snp.allele",
                paste0(mapSNVSel$REF, "/", mapSNVSel$ALT))

    if(verbose) { message("SNP part 1 DONE ", Sys.time()) }

    add.gdsn(node=gdsReference, name="snp.AF", as.numeric(mapSNVSel$AF),
                    storage="packedreal24")

    if(verbose) { message("SNP part AF DONE ", Sys.time()) }

    add.gdsn(node=gdsReference, name="snp.EAS_AF",
        val=as.numeric(mapSNVSel$EAS_AF), storage="packedreal24")
    add.gdsn(node=gdsReference, name="snp.EUR_AF",
        val=as.numeric(mapSNVSel$EUR_AF), storage="packedreal24")
    add.gdsn(node=gdsReference, name="snp.AFR_AF",
        val=as.numeric(mapSNVSel$AFR_AF), storage="packedreal24")
    add.gdsn(node=gdsReference, name="snp.AMR_AF",
        val=as.numeric(mapSNVSel$AMR_AF), storage="packedreal24")
    add.gdsn(node=gdsReference, name="snp.SAS_AF",
        val=as.numeric(mapSNVSel$SAS_AF), storage="packedreal24")

    ## Return OL when successful
    return(0L)
}


#' @title Add information related to profile genotype into a Reference GDS file
#'
#' @description This function adds the genotype fields with the associated
#' information into the Reference GDS file for the selected profiles.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Reference GDS file.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the 1K genotyping files for each sample are located. The name of the
#' genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#'
#' @param fileSNPsRDS a \code{character} string representing the path and file
#' name of the RDS file that contains the indexes of the retained SNPs. The
#' file must exist. The file must be a RDS file.
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to profiles (must be the profile identifiers) that will be
#' retained and added to the Reference GDS file.
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn write.gdsn
#' @importFrom utils read.csv2
#' @encoding UTF-8
#' @keywords internal
generateGDSgenotype <- function(gds, pathGeno, fileSNPsRDS, listSamples,
                                    verbose) {

    # File with the description of the SNP keep
    listMat1k <- dir(pathGeno, pattern=".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMat1k)

    listSNP <- readRDS(fileSNPsRDS)

    for(i in seq_len(length(listSamples))) {
        pos <- which(listSample1k == listSamples[i])

        if(verbose) { message(listSamples[i]) }

        if( length(pos) == 1) {
            matSample <- read.csv2(file.path(pathGeno, listMat1k[pos]),
                                        row.names=NULL)
            matSample <- matSample[listSNP,, drop=FALSE]
            if(i == 1) {
                var.geno <- add.gdsn(gds, "genotype",
                                valdim=c(nrow(matSample),
                                            length(listSamples)),
                                            storage="bit2")
            }

            # Easier to read
            matSample[matSample[,1] == "0|0", 1] <- 0
            matSample[matSample[,1] == "0|1" | matSample[,1] == "1|0", 1] <- 1
            matSample[matSample[,1] == "1|1", 1] <- 2

            g <- as.matrix(matSample)[,1]

            write.gdsn(var.geno, g, start=c(1, i), count=c(-1, 1))

            rm(matSample)

            if(verbose) { message(listMat1k[pos], " ", i) }
        }else{
            stop("Missing samples genotype in ", listSamples[i])
        }
    }

    return(0L)
}

#' @title This function append the field genotype in the gds file
#'
#' @description TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param pathGeno TODO a PATH to a directory with the a file for each
#' samples with the genotype.
#'
#' @param fileSNPsRDS TODO list of SNP to keep in the file genotype
#'
#' @param listSamples  a \code{array} with the sample to keep
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running. Default: \code{FALSE}.
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils read.csv2
#' @encoding UTF-8
#' @keywords internal
appendGDSgenotype <- function(gds, listSample, pathGeno, fileSNPsRDS,
                                verbose=FALSE) {

    # File with the description of the SNP keep
    listMat1k <- dir(pathGeno, pattern = ".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMat1k)

    listSNP <- readRDS(file=fileSNPsRDS)
    geno.var <- index.gdsn(gds, "genotype")
    g <- read.gdsn(node=geno.var, start=c(1, 1), count=c(1,-1))
    nbSample <- length(g)
    if(verbose) { message(nbSample) }
    for(i in seq_len(length(listSample))) {
        pos <- which(listSample1k == listSample[i])
        if( length(pos) == 1) {
            matSample <- read.csv2(file.path(pathGeno, listMat1k[pos]),
                                        row.names = NULL)
            matSample <- matSample[listSNP,, drop=FALSE]


            # Not faster but harder to read
            # matSample[,1] <- rowSums(t(matrix(as.numeric(unlist(strsplit(
            #                      matSample[,1], "\\|"))),nr=2)))
            # Easier to read
            matSample[matSample[,1] == "0|0",1] <- 0
            matSample[matSample[,1] == "0|1" | matSample[,1] == "1|0",1] <- 1
            matSample[matSample[,1] == "1|1",1] <- 2

            g <- as.matrix(matSample)[,1]
            append.gdsn(geno.var,g, check=TRUE)

            rm(matSample)
            if(verbose) { message(listMat1k[pos], " ", i) }
        }else {
            stop("Missing 1k samples ", listSample[i])
        }
    }

    return(0L)
}


#' @title create a file tfam file for plink from the gdsReference file
#'
#' @description TODO
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param listSample a \code{array} with the sample to keep TODO
#'
#' @param pedOUT TODO a PATH and file name to the output file
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
gds2tfam <- function(gdsReference, listSample, pedOUT) {

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

    write.table(dfPed, pedOUT,
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
#' @param pedOUT TODO a PATH and file name to the output file
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
gds2tfamSample <- function(gdsProfile, listSample, sampleANNO, pedOUT) {

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

    write.table(dfPed, pedOUT, quote=FALSE, sep="\t",
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
#' @param pedOUT TODO a PATH and file name to the output file
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
gds2tped <- function(gds, listSample, listSNP, pedOUT) {

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

    write.table(tped, pedOUT, quote=FALSE, sep="\t", row.names=FALSE,
                        col.names=FALSE)

}


#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened for the sample
#'
#' @param listBlock TODO
#'
#' @param blockName TODO
#'
#' @param blockDesc TODO
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn ls.gdsn compression.gdsn
#' @importFrom gdsfmt append.gdsn sync.gds
#' @encoding UTF-8
#' @keywords internal
addGDS1KGLDBlock <- function(gds, listBlock, blockName, blockDesc) {

    block.annot <- data.frame(block.id=blockName,
                                block.desc=blockDesc,
                                stringsAsFactors=FALSE)

    if(! ("block.annot" %in% ls.gdsn(gds))) {
        var.block.annot <- add.gdsn(gds, "block.annot", block.annot)
    }else {
        curAnnot <- index.gdsn(gds, "block.annot/block.id")
        append.gdsn(curAnnot,block.annot$block.id)
        curAnnot <- index.gdsn(gds, "block.annot/block.desc")
        append.gdsn(curAnnot, block.annot$block.desc)
    }

    var.block <- NULL
    if(! ("block" %in% ls.gdsn(gds))){
        var.block <- add.gdsn(gds, "block",
                                valdim=c(length(listBlock), 1),
                                listBlock, storage="int32",
                                compress = "LZ4_RA")
        readmode.gdsn(var.block)

    }else {
        if(is.null(var.block)) {
            var.block <- index.gdsn(gds, "block")
            var.block <- compression.gdsn(var.block, "")
        }
        append.gdsn(var.block, listBlock)
        var.block <- compression.gdsn(var.block, "LZ4_RA")
    }

    sync.gds(gds)

    return(0L)
}


#' @title Add information related to segments associated to the SNV
#' dataset for a specific sample into a GDS file
#'
#' @description The function adds the information related to segments
#' associated to the SNV dataset for a specific sample into a
#' GDS file, more specifically, in the "segment" node. If the "segment" node
#' already exists, the previous information is erased.
#'
#' @param gdsProfile an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), a GDS Sample file.
#'
#' @param snp.seg a \code{vector} of \code{integer} representing the segment
#' identifiers associated to each SNV selected for the specific sample. The
#' length of the \code{vector} should correspond to the number of SNVs
#' present in the "snp.id" entry of the GDS sample file.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(dataDir, "GDS_TEMP.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Vector of segment identifiers
#' segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L)
#'
#' ## Add segments to the GDS file
#' RAIDS:::addUpdateSegment(gdsProfile=GDS_file_tmp, snp.seg=segments)
#'
#' ## Read segments information from GDS file
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="segment"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @keywords internal
addUpdateSegment <- function(gdsProfile, snp.seg) {

    if("segment" %in% ls.gdsn(gdsProfile)) {
        snpLap <- write.gdsn(index.gdsn(gdsProfile, "segment"), snp.seg)
    } else{
        snpLap <- add.gdsn(gdsProfile, "segment", snp.seg, storage="uint32")
    }

    sync.gds(gdsProfile)

    ## Successful
    return(0L)
}


#' @title Get the block number for each SNV in snp.index
#'
#' @description TODO
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened 1KG SNV Annotation GDS file. RNA specific
#' Default: \code{NULL}.
#'
#' @param snp.index TODO
#'
#' @param blockID a \code{character} string corresponding to the block
#' identifier in \code{gdsRefAnnot}.  RNA specific
#' Default: \code{NULL}
#'
#' @return TODO a \code{vector} of \code{numeric} corresponding to the
#' block identifier
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @keywords internal
getGeneBlock <- function(gdsRefAnnot, snp.index, blockID) {

    block.annot <- read.gdsn(index.gdsn(gdsRefAnnot, "block.annot"))
    pos <- which(block.annot$block.id == blockID)

    if(length(pos) != 1) {
        stop("Try to get Gene.Block with blockID problematic ", blockID)
    }

    b <- read.gdsn(index.gdsn(gdsRefAnnot, "block"), start=c(1,pos),
                    count = c(-1,1))[snp.index]

    return(b)
}


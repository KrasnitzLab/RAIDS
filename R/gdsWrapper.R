#' @title Add information related to SNVs into a Population Reference GDS file
#'
#' @description The function adds the SNV information into a Population
#' Reference GDS file.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Reference GDS file.
#'
#' @param fileFreq a \code{character} string representing the path and file
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
#' library(gdsfmt)
#'
#' ## The RDS file containing the filtered SNP information
#' dataDir <- system.file("extdata", package="RAIDS")
#' fileFilerterSNVs <- file.path(dataDir, "mapSNVSelected_Demo.rds")
#'
#' ## Temporary Reference GDS file in temporary directory
#' file1KG <- file.path(tempdir(), "1KG_TEMP_002.gds")
#' filenewGDS <- createfn.gds(file1KG)
#'
#' ## Add SNV information to Reference GDS
#' RAIDS:::generateGDSSNPinfo(gdsReference=filenewGDS,
#'             fileFreq=fileFilerterSNVs, verbose=TRUE)
#'
#' ## Close GDS file (important)
#' closefn.gds(filenewGDS)
#'
#' ## Remove temporary 1KG_TEMP_002.gds file
#' unlink(file1KG, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
generateGDSSNPinfo <- function(gdsReference, fileFreq, verbose) {

    mapSNVSel <- readRDS(file=fileFreq)

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


#' @title Add information related to profile genotypes into a Population
#' Reference GDS file
#'
#' @description This function adds the genotype fields with the associated
#' information into the Population Reference GDS file for the selected
#' profiles.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Population Reference
#' GDS file.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the reference genotyping files for each sample are located. The name of the
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
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Path where the demo genotype CSV files are located
#' pathGeno <- file.path(dataDir, "demoProfileGenotypes")
#'
#' ## The RDS file containing the pedigree information
#' pedigreeFile <- file.path(dataDir, "PedigreeReferenceDemo.rds")
#'
#' ## The RDS file containing the indexes of the retained SNPs
#' snpIndexFile <- file.path(dataDir, "listSNPIndexes_Demo.rds")
#'
#' ## The RDS file containing the filtered SNP information
#' filterSNVFile <- file.path(dataDir, "mapSNVSelected_Demo.rds")
#'
#' ## Temporary Reference GDS file in temporary directory
#' tempRefGDS <- file.path(tempdir(), "Ref_TEMP01.gds")
#'
#' ## Create temporary Reference GDS file
#' newGDS <- createfn.gds(tempRefGDS)
#' put.attr.gdsn(newGDS$root, "FileFormat", "SNP_ARRAY")
#'
#' ## Read the pedigree file
#' ped1KG <- readRDS(pedigreeFile)
#'
#' ## Add information about samples to the Reference GDS file
#' listSampleGDS <- RAIDS:::generateGDSRefSample(gdsReference=newGDS,
#'                 dfPedReference=ped1KG, listSamples=NULL)
#'
#' ## Add SNV information to the Reference GDS
#' RAIDS:::generateGDSSNPinfo(gdsReference=newGDS, fileFreq=filterSNVFile,
#'                 verbose=FALSE)
#'
#' ## Add genotype information to the Reference GDS
#' RAIDS:::generateGDSgenotype(gds=newGDS, pathGeno=pathGeno,
#'         fileSNPsRDS=snpIndexFile, listSamples=listSampleGDS, verbose=FALSE)
#'
#' ## Close file
#' closefn.gds(newGDS)
#'
#' ## Remove temporary files
#' unlink(tempRefGDS, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn write.gdsn
#' @importFrom utils read.csv2
#' @encoding UTF-8
#' @keywords internal
generateGDSgenotype <- function(gds, pathGeno, fileSNPsRDS, listSamples,
                                    verbose) {

    # File with the description of the SNP keep
    listMatRef <- dir(pathGeno, pattern=".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMatRef)

    listSNP <- readRDS(fileSNPsRDS)

    for(i in seq_len(length(listSamples))) {
        pos <- which(listSample1k == listSamples[i])

        if(verbose) { message(listSamples[i]) }

        if( length(pos) == 1) {
            matSample <- read.csv2(file.path(pathGeno, listMatRef[pos]),
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

            if(verbose) { message(listMatRef[pos], " ", i) }
        }else{
            stop("Missing samples genotype in ", listSamples[i])
        }
    }

    return(0L)
}

#' @title Append information related to profile genotypes into a Population
#' Reference GDS file (associated node already present in the GDS)
#'
#' @description This function appends the genotype fields with the associated
#' information into the Population Reference GDS file for the selected
#' profiles. The associated node must already present in the GDS file.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Population Reference
#' GDS file.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the reference genotyping files for each sample are located. The name of the
#' genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#'
#' @param fileSNPsRDS a \code{character} string representing the path and file
#' name of the RDS file that contains the indexes of the retained SNPs. The
#' file must exist. The file must be a RDS file.
#'
#' @param listSamples a \code{character} string representing the path and file
#' name of the RDS file that contains the indexes of the retained SNPs. The
#' file must exist. The file must be a RDS file.
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running.
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Path where the demo genotype CSV files are located
#' pathGeno <- file.path(dataDir, "demoProfileGenotypes")
#'
#' ## The RDS file containing the pedigree information
#' pedigreeFile <- file.path(dataDir, "PedigreeReferenceDemo.rds")
#'
#' ## The RDS file containing the indexes of the retained SNPs
#' snpIndexFile <- file.path(dataDir, "listSNPIndexes_Demo.rds")
#'
#' ## The RDS file containing the filtered SNP information
#' filterSNVFile <- file.path(dataDir, "mapSNVSelected_Demo.rds")
#'
#' ## Temporary Reference GDS file
#' tempRefGDS <- file.path(tempdir(), "Ref_TEMP02.gds")
#'
#' ## Create temporary Reference GDS file
#' newGDS <- createfn.gds(tempRefGDS)
#' put.attr.gdsn(newGDS$root, "FileFormat", "SNP_ARRAY")
#'
#' ## Read the pedigree file
#' ped1KG <- readRDS(pedigreeFile)
#'
#' ## Add information about samples to the Reference GDS file
#' listSampleGDS <- RAIDS:::generateGDSRefSample(gdsReference=newGDS,
#'                 dfPedReference=ped1KG, listSamples=NULL)
#'
#' ## Add SNV information to the Reference GDS
#' RAIDS:::generateGDSSNPinfo(gdsReference=newGDS, fileFreq=filterSNVFile,
#'                 verbose=FALSE)
#'
#' ## Add genotype information to the Reference GDS for the 3 first samples
#' RAIDS:::generateGDSgenotype(gds=newGDS, pathGeno=pathGeno,
#'         fileSNPsRDS=snpIndexFile, listSamples=listSampleGDS[1:3],
#'         verbose=FALSE)
#'
#' ## Append genotype information to the Reference GDS for the other samples
#' RAIDS:::appendGDSgenotype(gds=newGDS, pathGeno=pathGeno,
#'         fileSNPsRDS=snpIndexFile,
#'         listSample=listSampleGDS[4:length(listSampleGDS)],
#'         verbose=FALSE)
#'
#' ## Close file
#' closefn.gds(newGDS)
#'
#' ## Remove temporary files
#' unlink(tempRefGDS, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils read.csv2
#' @encoding UTF-8
#' @keywords internal
appendGDSgenotype <- function(gds, listSample, pathGeno, fileSNPsRDS,
                                verbose) {

    # File with the description of the SNP keep
    listMatRef <- dir(pathGeno, pattern=".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMatRef)

    listSNP <- readRDS(file=fileSNPsRDS)
    geno.var <- index.gdsn(gds, "genotype")
    g <- read.gdsn(node=geno.var, start=c(1, 1), count=c(1,-1))
    nbSample <- length(g)
    if(verbose) { message(nbSample) }
    for(i in seq_len(length(listSample))) {
        pos <- which(listSample1k == listSample[i])
        if( length(pos) == 1) {
            matSample <- read.csv2(file.path(pathGeno, listMatRef[pos]),
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
            append.gdsn(geno.var, g, check=TRUE)

            rm(matSample)
            if(verbose) { message(listMatRef[pos], " ", i) }
        }else {
            stop("Missing reference samples ", listSample[i])
        }
    }

    return(0L)
}


#' @title Add block information in a Population Reference GDS Annotation file
#'
#' @description This function appends the information for one specific type
#' of blocks into a Population Reference GDS Annotation file. More
#' specifically, the node 'block.annot' is created if it does not exists. This
#' node contains a \code{data.frame} which will be append the description of
#' the current block. The node 'block' is also created if it does not exists.
#' This node is a \code{matrix} that will contain all the entries for the
#' current block. All the values for a specific block type are contained in a
#' single column that corresponds to the row number in the 'block.annot' node.
#'
#' @param gds an object of class \code{gds} opened in writing mode.
#'
#' @param listBlock a \code{array} of \code{integer} representing all the
#' entries for the current block.
#'
#' @param blockName a \code{character} string representing the unique
#' block name.
#'
#' @param blockDesc a \code{character} string representing the description of
#' the current block.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Temporary GDS Annotation file in current directory
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_Annot_14.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## One block
#' blockType <- "EAS.0.05.500k"
#'
#' ## The description of the block
#' blockDescription <- "EAS population blocks based on 500k windows"
#'
#' ## The values for each entry related to the block (integers)
#' blockEntries <- c(1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3)
#'
#' RAIDS:::addBlockInGDSAnnot(gds=GDS_file_tmp, listBlock=blockEntries,
#'     blockName=blockType, blockDesc=blockDescription)
#'
#' ## Read 'block.annot' node
#' read.gdsn(index.gdsn(GDS_file_tmp, "block.annot"))
#'
#' ## Read 'block' node
#' read.gdsn(index.gdsn(GDS_file_tmp, "block"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn ls.gdsn compression.gdsn
#' @importFrom gdsfmt append.gdsn sync.gds
#' @encoding UTF-8
#' @keywords internal
addBlockInGDSAnnot <- function(gds, listBlock, blockName, blockDesc) {

    blockAnnot <- data.frame(block.id=blockName,
                                block.desc=blockDesc,
                                stringsAsFactors=FALSE)

    if(! ("block.annot" %in% ls.gdsn(gds))) {
        ## Create 'block.annot' node when not existing
        varBlockAnnot <- add.gdsn(gds, "block.annot", blockAnnot)
    }else {
        ## Append content to 'block.annot' node when existing
        curAnnot <- index.gdsn(gds, "block.annot/block.id")
        append.gdsn(curAnnot, blockAnnot$block.id)
        curAnnot <- index.gdsn(gds, "block.annot/block.desc")
        append.gdsn(curAnnot, blockAnnot$block.desc)
    }

    varBlock <- NULL
    if(!("block" %in% ls.gdsn(gds))) {
        ## Create 'block' node that will contain a matrix of integers
        ## stored in compressed mode
        varBlock <- add.gdsn(node=gds, name="block",
                                valdim=c(length(listBlock), 1),
                                listBlock, storage="int32",
                                compress="LZ4_RA")
        readmode.gdsn(varBlock)

    } else {
        if(is.null(varBlock)) {
            varBlock <- index.gdsn(gds, "block")
            varBlock <- compression.gdsn(varBlock, "")
        }
        append.gdsn(varBlock, listBlock)
        ## Compressed data using LZ4_RA method
        varBlock <- compression.gdsn(varBlock, "LZ4_RA")
    }

    sync.gds(gds)

    return(0L)
}





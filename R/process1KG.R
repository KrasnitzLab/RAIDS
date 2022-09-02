#' @title Prepare the pedigree file using pedigree information from 1KG
#'
#' @description Using the pedigree file from 1KG, this function extracts
#' needed information and formats it into a \code{data.frame} so in can
#' be used in following steps of the ancestry inference process. The
#' function also requires that the genotyping files associated to each
#' sample be available in a specified directory.
#'
#' @param pedFile a \code{character} string representing the path and
#' file name of the pedigree file from 1KG. The file must exist.
#'
#' @param PATHGENO a \code{character} string representing the path where
#' the 1K genotyping files for each sample are located. Only the samples with
#' associated genotyping files are retained in the creation of the final
#' \code{data.frame}. The name of the genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#' Default: \code{"./data/sampleGeno"}.
#'
#' @param batch.v a\code{integer} that uniquely identifies the source of the
#' pedigree information. The 1KG is usually \code{0L}. Default: \code{0L}.
#'
#' @return A \code{data.frame} containing the needed pedigree information
#' from 1K.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo pedigree file
#' pedDemoFile <- file.path(data.dir, "PedigreeDemo.ped")
#'
#' ## Create a data.frame containing the information of the retained
#' ## samples (samples with existing genotyping files)
#' prepPed1KG(pedFile=pedDemoFile, PATHGENO=data.dir, batch.v=0L)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils read.delim
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
prepPed1KG <- function(pedFile, PATHGENO=file.path("data", "sampleGeno"),
                        batch.v=0L) {

    ## Validate that the batch is an integer
    if (! isSingleNumber(batch.v)) {
        stop("The batch.v must be an integer.")
    }

    ## Validate that the pedigree file exists
    if (! file.exists(pedFile)) {
        stop("The file \'", pedFile, "\' does not exist." )
    }

    ## Validate that the path for the genotyping files exists
    if (! file.exists(PATHGENO)) {
        stop("The path \'", PATHGENO, "\' does not exist." )
    }

    ## Read the pedigree file from 1KG
    ped1KG <- read.delim(pedFile)

    ## Create a data.frame containing the needed information
    pedAll <- data.frame(
                sample.id=c(ped1KG$Individual.ID),
                Name.ID=c(ped1KG$Individual.ID),
                sex=c(ped1KG$Gender),
                pop.group=c(ped1KG$Population),
                superPop=rep(NA, length(c(ped1KG$Population))),
                batch=c(rep(batch.v, nrow(ped1KG))),
                stringsAsFactors=FALSE
    )

    ## Create a list with all populations associated to each super-population
    ## TODO The population versus super.population is hard-coded
    ## TODO change to parameters
    listSuperPop1000G <- list()
    listSuperPop1000G[['EAS']] <- c("CHB", "JPT", "CHS", "CDX", "KHV")
    listSuperPop1000G[['EUR']] <- c("CEU", "TSI", "FIN", "GBR", "IBS")
    listSuperPop1000G[['AFR']] <- c("YRI", "LWK", "GWD", "MSL", "ESN",
                                        "ASW", "ACB")
    listSuperPop1000G[['AMR']] <- c("MXL", "PUR", "CLM", "PEL")
    listSuperPop1000G[['SAS']] <- c("GIH", "PJL", "BEB", "STU", "ITU")


    ## Identify the super-population associated to each sample in
    ## the data.frame
    listSuperPop <- c("EAS", "EUR", "AFR", "AMR", "SAS")
    for(sp in listSuperPop){
        pedAll[which(pedAll$pop.group %in% listSuperPop1000G[[sp]]),
                "superPop"] <- sp
    }

    ## Assign row names to the data.frame
    row.names(pedAll) <- pedAll$sample.id

    ## Change column format for Sex information
    ## TODO: could be done when the data.frame is created
    pedAll$sex <- as.character(pedAll$sex)


    ## Only retained samples with existing genotyping file
    listMat1k <- dir(PATHGENO, pattern = ".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMat1k)

    pedAll <- pedAll[listSample1k, ]

    return(pedAll)
}


#' @title Generate the filter SNP information file in RDS format
#'
#' @description The function applies a cut-off filter to the SNP information
#' file to retain only the SNP that have a frequency superior or equal to the
#' specified cut-off in at least one super population. The information about
#' the retained SNPs is saved in a RDS format file. A RDS file containing the
#' indexes of the retained SNP is also created.
#'
#' @param cutOff a single \code{numeric} value, the cut-off
#' for the frequency in at least one super population. Default: \code{0.01}.
#'
#' @param fileSNV  a \code{character} string representing the path and
#' file name of the bulk SNP information file from 1KG. The file must be in
#' text format. The file must exist.
#'
#' @param fileLSNP a \code{character} string representing the path and
#' file name of the RDS file that will contain the indexes of the retained
#' SNPs. The file extension must be '.rds'.
#'
#' @param fileFREQ a \code{character} string representing the path and
#' file name of the RDS file that will contain the filtered SNP information.
#' The file extension must be '.rds'.
#'
#' @return The integer \code{0} when successful.
#'
#' @details
#'
#' The filtered SNP information RDS file (parameter \code{fileFREQ}), contains
#' a \code{data.frame} with those columns:
#' \itemize{
#' \item{CHROM} {a \code{character} string representing the chromosome where
#' the SNV is located.}
#' \item{POS} {a \code{character} string representing the SNV position on the
#' chromosome.}
#' \item{REF} {a \code{character} string representing the reference DNA base
#' for the SNV.}
#' \item{ALT} {a \code{character} string representing the alternative DNA base
#' for the SNV.}\
#' \item{EAS_AF} {a \code{character} string representing the allele frequency
#' of the EAS super population.}
#' \item{AFR_AF} {a \code{character} string representing the allele frequency
#' of the AFR super population.}
#' \item{AMR_AF} {a \code{character} string representing the allele frequency
#' of the AMR super population.}
#' \item{SAS_AF} {a \code{character} string representing the allele frequency
#' of the SAS super population.}
#' }
#'
#' @examples
#'
#' ## Needed package
#' library(withr)
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo SNV information file used as input
#' snvFile <- file.path(data.dir, "matFreqSNV_Demo.txt.bz2")
#'
#' ## Temporary output files
#' ## The first file contains the indexes of the retained SNPs
#' ## The second file contains the filtered SNP information
#' snpIndexFile <- local_file(file.path(data.dir, "listSNP_TEMP.rds"))
#' filterSNVFile <- local_file(file.path(data.dir, "mapSNVSel_TEMP.rds"))
#'
#' ## Create a data.frame containing the information of the retained
#' ## samples (samples with existing genotyping files)
#' generateMapSnvSel(cutOff=0.01, fileSNV=snvFile,
#'     fileLSNP=snpIndexFile, fileFREQ=filterSNVFile)
#'
#' ## Remove temporary files
#' deferred_run()
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
generateMapSnvSel <- function(cutOff=0.01, fileSNV, fileLSNP, fileFREQ) {

    ## The cut-off must be a single number
    if (!isSingleNumber(cutOff)) {
        stop("The cutOff must be a single numeric value.")
    }

    ## Validate that the bulk SNP file exists
    if (! file.exists(fileSNV)) {
        stop("The file \'", fileSNV, "\' does not exist." )
    }

    ## Read the bulk SNP information file
    mapSNVSel <- read.csv2(fileSNV)

    ## Identify SNPs that have a frequency equal or superior to the cut-off
    ## in at least one super population
    listSNP <- which(rowSums(mapSNVSel[,c("EAS_AF",
                                            "EUR_AF",
                                            "AFR_AF",
                                            "AMR_AF",
                                            "SAS_AF")] >= (cutOff - 1e-10) &
                                mapSNVSel[,c("EAS_AF",
                                                "EUR_AF",
                                                "AFR_AF",
                                                "AMR_AF",
                                                "SAS_AF")] <= (1 - cutOff +
                                                                1e-10)) > 0)

    ## Filter the data.frame
    mapSNVSel <- mapSNVSel[listSNP,]

    ## Save the indexes of the retained SNPs in RDS
    saveRDS(listSNP, fileLSNP)
    ## Save the filtered SNP information file in RDS
    saveRDS(mapSNVSel, fileFREQ)

    return(0L)
}


#' @title Generate the GDS file that will contain the information from 1KG
#'
#' @description This function generates the GDS file that will contain the
#' information from 1KG. The function also add the samples information, the
#' SNP information and the genotyping information into the GDS file.
#'
#' @param PATHGENO a \code{character} string representing the path where
#' the 1K genotyping files for each sample are located. The name of the
#' genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#' Default: \code{"./data/sampleGeno"}.
#'
#' @param fileNamePED a \code{character} string representing the path and file
#' name of the RDS file that contains the pedigree information. The file must
#' exist. The file must be a RDS file.
#'
#' @param fileListSNP a \code{character} string representing the path and file
#' name of the RDS file that contains the indexes of the retained SNPs. The
#' file must exist. The file must be a RDS file.
#'
#' @param fileSNPSel a \code{character} string representing the path and file
#' name of the RDS file that contains the filtered SNP information. The
#' file must exist. The file must be a RDS file.
#'
#' @param fileNameGDS a \code{character} string representing the path and file
#' name of the GDS file that will be created. The GDS file will contain the
#' SNP information, the genotyping information and the pedigree information
#' from 1000 Genomes.
#' The extension of the file must be '.gds'.
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to samples (must be the sample.ids) that will be retained and added to the
#' GDS file. When \code{NULL}, all the samples are retained.
#' Default: \code{NULL}.
#'
#' @param verbose a \code{logical} indicating if the funciton must print
#' messages when running. Default: \code{FALSE}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @details
#'
#' More information about GDS file format can be found at the Bioconductor
#' gdsfmt website:
#' https://bioconductor.org/packages/gdsfmt/
#'
#' @examples
#'
#' ## Needed package
#' library(withr)
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The RDS file containing the pedigree information
#' pedigreeFile <- file.path(data.dir, "PedigreeDemo.rds")
#'
#' ## The RDS file containing the indexes of the retained SNPs
#' snpIndexFile <- file.path(data.dir, "listSNPIndexes_Demo.rds")
#'
#' ## The RDS file containing the filtered SNP information
#' filterSNVFile <- file.path(data.dir, "mapSNVSelected_Demo.rds")
#'
#' ## Temporary GDS file containing 1KG information
#' GDS_file <- local_file(file.path(data.dir, "1KG_TEMP.gds"))
#'
#' ## Create a temporary GDS file containing information from 1KG
#' generateGDS1KG(PATHGENO=data.dir, fileNamePED=pedigreeFile,
#'     fileListSNP=snpIndexFile, fileSNPSel=filterSNVFile,
#'     fileNameGDS=GDS_file, listSamples=NULL)
#'
#' ## Remove temporary files
#' deferred_run()
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#'
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds
#' @encoding UTF-8
#' @export
generateGDS1KG <- function(PATHGENO=file.path("data", "sampleGeno"),
                            fileNamePED, fileListSNP,
                            fileSNPSel, fileNameGDS,
                            listSamples=NULL, verbose=FALSE) {

    ## Validate that the pedigree file exists
    if (! file.exists(fileNamePED)) {
        stop("The file \'", fileNamePED, "\' does not exist." )
    }

    ## Validate that the path for the genotyping files exists
    if (! file.exists(PATHGENO)) {
        stop("The path \'", PATHGENO, "\' does not exist." )
    }

    ## Validate that the SNP indexes file exists
    if (! file.exists(fileListSNP)) {
        stop("The file \'", fileListSNP, "\' does not exist." )
    }

    ## Validate that the SNP information file exists
    if (! file.exists(fileSNPSel)) {
        stop("The file \'", fileSNPSel, "\' does not exist." )
    }

    ## Read the pedigree file
    ped1KG <- readRDS(fileNamePED)

    # list in the file genotype we keep from fileLSNP in generateMapSnvSel
    listKeep <- readRDS(fileListSNP)

    # Create the GDS file
    newGDS <- createfn.gds(fileNameGDS)
    put.attr.gdsn(newGDS$root, "FileFormat", "SNP_ARRAY")

    if(verbose) { message("Start ", Sys.time()) }

    listSampleGDS <- generateGDSSample(newGDS, ped1KG, listSamples)

    if(verbose) { message("Sample info DONE ", Sys.time()) }

    generateGDSSNPinfo(newGDS, fileSNPSel)

    if(verbose) { message("SNP info DONE ", Sys.time()) }

    generateGDSgenotype(newGDS, PATHGENO, fileListSNP, listSampleGDS)

    if(verbose) { message("Genotype DONE ", Sys.time()) }

    closefn.gds(newGDS)

    return(0L)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened
#'
#' @param gdsPhase TODO
#'
#' @param PATHGENO a \code{character} string representing the path where
#' the 1K genotyping files for each sample are located. The name of the
#' genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#' Default: \code{"./data/sampleGeno"}.
#'
#' @param fileLSNP TODO
#'
#' @param verbose a \code{logicial} indicating if the function should
#' print messages when running. Default: \code{FALSE}.
#'
#' @return TODO a \code{vector} of \code{string}
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn readmode.gdsn
#' @encoding  UTF-8
#' @export
generatePhase1KG2GDS <- function(gds, gdsPhase,
                            PATHGENO, fileLSNP, verbose=FALSE) {

    ## The verbose parameter must be a logical
    if (!(is.logical(verbose) && length(verbose) == 1)) {
        stop("The \'verbose\' parameters must be a single logical value ",
             "(TRUE or FALSE).")
    }

    sample.id <- read.gdsn(index.gdsn(gds,"sample.id"))
    listSNP <- readRDS(fileLSNP)

    var.phase <- NULL
    for(i in seq_len(length(sample.id))){

        if (verbose) { message("S ", i, " ", Sys.time()) }

        file1KG <- file.path(PATHGENO, paste0(sample.id[i],".csv.bz2"))
        matSample <- read.csv2( file1KG,
                                row.names = NULL)[listSNP,, drop=FALSE]
        matSample <- matrix(as.numeric(unlist(strsplit(matSample[, 1],
                                                      "\\|"))), nrow=2)[1,]

        if (verbose) { message("GDS ", i, " ", Sys.time()) }

        if(! ("phase" %in% ls.gdsn(gdsPhase))){
            var.phase <- add.gdsn(gdsPhase, "phase",
                                    valdim=c(length(listSNP), 1),
                                    matSample,
                                    storage="bit2",
                                    compress = "LZ4_RA")

        }else{
            if(is.null(var.phase)){
                var.phase <- index.gdsn(gdsPhase, "phase")
            }
            append.gdsn(var.phase, matSample)
        }
        rm(matSample)
    }

    readmode.gdsn(var.phase)

    return(0L)
}

#' @title Identify genetically unrelated patients in GDS 1KG file
#'
#' @description The function identify patients that are genetically related in
#' the 1KG file. It generates a first RDS file with the list of unrelated
#' patient. It also generates a second RDS file with the kinship coefficient
#' between the patients.
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, the
#' 1KG GDS file.
#'
#' @param maf a single \code{numeric} representing the threshold for the minor
#' allele frequency. Only the SNPs with ">= maf" will be used.
#' Default: \code{0.05}.
#'
#' @param thresh a single \code{numeric} representing the threshold value used
#' to decide if a pair of individuals is ancestrally divergent.
#' Default: \code{2^(-11/2)}.
#'
#' @param fileIBD a \code{character} string representing the path and file
#' name of the RDS file that will be created. The RDS file will contain the
#' kinship coefficient between the patients.
#' The extension of the file must be '.rds'.
#'
#' @param filePart a \code{character} string representing the path and file
#' name of the RDS file that will be created. The RDS file will contain the
#' information about the 1KG patients that are unrelated. The file will
#' contains two lists: the \code{list} of related samples, called \code{rels}
#' and the list of unrelated samples, called \code{unrels}.
#' The extension of the file must be '.rds'.
#'
#' @return \code{NULL} invisibly.
#'
#' @examples
#'
#' ## Needed packages
#' library(withr)
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## Open existing 1K GDS file
#' GDS_file <- file.path(data.dir, "1KG_Demo.gds")
#' gdsFile <- snpgdsOpen(GDS_file)
#'
#' ## Temporary output files
#' ## The first RDS file will contain the list of unrelated patients
#' ## The second RDS file will contain the kinship information between patients
#' patientTmpFile <- local_file(file.path(data.dir,
#'     "unrelatedPatients_TEMP.rds"))
#' ibdTmpFile <- local_file(file.path(data.dir,"ibd_TEMP.rds"))
#'
#' ## Identify unrelated patients in 1KG GDS file
#' identifyRelative(gds=gdsFile, maf=0.05, thresh=2^(-11/2),
#'     fileIBD=ibdTmpFile, filePart=patientTmpFile)
#'
#' ## Close 1K GDS file
#' closefn.gds(gdsFile)
#'
#' ## Remove temporary files
#' deferred_run()
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#'
#' @importFrom GENESIS pcairPartition
#' @importFrom S4Vectors isSingleNumber
#' @importFrom methods is
#' @encoding UTF-8
#' @export
identifyRelative <- function(gds, maf=0.05, thresh=2^(-11/2),
                                fileIBD, filePart) {

    ## Validate that the maf is a single numeric value
    if (! isSingleNumber(maf)) {
        stop("The \'maf\' parameter must be a single numeric value." )
    }

    ## Validate that the thresh is a single numeric value
    if (! isSingleNumber(thresh)) {
        stop("The \'thresh\' parameter must be a single numeric value." )
    }

    ## Validate that gds is an object of class SNPGDSFileClass
    if (! is(gds, "SNPGDSFileClass")) {
        stop("The \'gds\' parameter must be an object of ",
                    "class \'SNPGDSFileClass\'.")
    }

    ibd.robust <- runIBDKING(gds=gds, maf=maf)

    matKING <- ibd.robust$kinship
    colnames(matKING) <- ibd.robust$sample.id
    row.names(matKING) <- ibd.robust$sample.id

    part <- pcairPartition(kinobj=matKING, divobj=matKING,
                                kin.thresh=thresh, div.thresh=-1*thresh)

    saveRDS(ibd.robust, fileIBD)
    saveRDS(part, filePart)
}


#' @title Add the information about the unrelated patients to the GDS 1KG file
#'
#' @description This function adds the information about the unrelated patients
#' to the GDS 1KG file. More specificaly, it creates the field
#' \code{sample.ref} which as the value \code{1} when the sample
#' is unrelated and the value \code{0} otherwise.
#' The \code{sample.ref} is filled based on the information present in the
#' input RDS file.
#'
#' @param fileNameGDS a \code{character} string representing the path and file
#' name of the GDS file that contains the 1KG information. The GDS file must
#' contain the SNP information, the genotyping information and
#' the pedigree information from 1000 Genomes.
#' The extension of the file must be '.gds'.
#'
#' @param filePart a \code{character} string representing the path and file
#' name of the RDS file that contains the
#' information about the 1KG patients that are unrelated.
#' The extension of the file must be '.rds'. The file must exists.
#'
#' @return \code{None}
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## 1K GDS file
#' GDS_file <- file.path(data.dir, "1KG_Demo.gds")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @export
addRef2GDS1KG <- function(fileNameGDS, filePart) {

    ## Validate that the GDS file exists
    if (! file.exists(fileNameGDS)) {
        stop("The file \'", fileNameGDS, "\' does not exist." )
    }

    ## Validate that the RDS file exists
    if (! file.exists(filePart)) {
        stop("The file \'", filePart, "\' does not exist." )
    }

    gds <- snpgdsOpen(fileNameGDS, readonly=FALSE)

    addGDSRef(gds, filePart)

    closefn.gds(gds)
}


#' @title Find the pruned snv in 1KG by chr
#'
#' @description TODO
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param method a \code{character string} TODO . Default: \code{"corr"}.
#'
#' @param listSamples TODO
#'
#' @param slide.max.bp.v TODO
#'
#' @param ld.threshold.v TODO
#'
#' @param np TODO . Default: \code{NULL}.
#'
#' @param verbose.v a \code{logical} specifying if the function must provide
#' more information about the process. Default: \code{FALSE}.
#'
#' @param chr TODO
#'
#' @param minAF TODO
#'
#' @param outPref TODO
#'
#' @param keepObj a \code{logical} specifying if the function must save the
#' the processed information into a RDS object. Default: \code{FALSE}.
#'
#' @return \code{NULL} invisibly.
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
#' @encoding UTF-8
#' @keywords internal
pruning1KG.Chr <- function(gds, method="corr",
                            listSamples=NULL,
                            slide.max.bp.v=5e5,
                            ld.threshold.v=sqrt(0.1),
                            np=1, verbose.v=FALSE,
                            chr=NULL,
                            minAF=NULL,
                            outPref="pruned_1KG",
                            keepObj=FALSE) {

    filePruned <- file.path(paste0(outPref, ".rds"))
    fileObj <- file.path(paste0(outPref, "Obj.rds"))
    snpGDS <- index.gdsn(gds, "snp.id")
    listKeep <- NULL
    if(is.null(minAF)){
        if(!is.null(chr)){
            snpGDS <- index.gdsn(gds, "snp.id")
            snpID <- read.gdsn(snpGDS)

            chrGDS <- index.gdsn(gds, "snp.chromosome")
            snpCHR <- read.gdsn(chrGDS)

            listKeep <- snpID[which(snpCHR == chr)]
        }
    } else{
        snpGDS <- index.gdsn(gds, "snp.id")
        snpID <- read.gdsn(snpGDS)
        afGDS <- index.gdsn(gds, "snp.AF")
        snpAF <- read.gdsn(afGDS)

        if(is.null(chr)){
            listKeep <- snpID[which(snpAF >= minAF & snpAF <= 1-minAF)]
        } else{
            chrGDS <- index.gdsn(gds, "snp.chromosome")
            snpCHR <- read.gdsn(chrGDS)

            listKeep <- snpID[which(snpCHR == chr & snpAF >= minAF &
                                        snpAF <= 1-minAF)]
        }
    }

    snpset <- runLDPruning(gds,
                            method,
                            listSamples=listSamples,
                            listKeep=listKeep,
                            slide.max.bp.v = slide.max.bp.v,
                            ld.threshold.v=ld.threshold.v)

    pruned <- unlist(snpset, use.names=FALSE)
    saveRDS(pruned, filePruned)
    if(keepObj){
        saveRDS(snpset, fileObj)
    }
}



#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param listSample.Ref  A \code{vector} of \code{string} corresponding to
#' the sample.ids
#'
#' @param listSNP the list of snp.id keep
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1}.
#'
#' @return TODO a \code{list}  with with two objects
#' pca.unrel -> \code{snpgdsPCAClass}
#' and a snp.load -> \code{snpgdsPCASNPLoading}
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASNPLoading
#' @encoding UTF-8
#' @export
basePCASample <- function(gds, listSample.Ref=NULL, listSNP=NULL, np=1) {

    listPCA <- list()

    listPCA[["SNP"]] <- listSNP

    listPCA[["pca.unrel"]] <- snpgdsPCA(gds, sample.id=listSample.Ref,
                                            snp.id=listSNP,
                                            num.thread=np,
                                            verbose=TRUE)

    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                                    gdsobj=gds,
                                                    num.thread=np,
                                                    verbose=TRUE)

    return(listPCA)
}

#' @title TODO contain the information from 1KG
#'
#' @description TODO
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), TODO
#'
#' @param gdsOut an object of class \code{gds} in writing
#'
#' @param PATHBLOCK TODO
#'
#' @param superPop TODO
#'
#' @param blockName TODO
#'
#' @param blockDesc TODO
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{TRUE}.
#'
#' @return \code{OL} when the function is successful.
#'
#' @details
#'
#' More information about GDS file format can be found at the Bioconductor
#' gdsfmt website:
#' https://bioconductor.org/packages/gdsfmt/
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#'
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds
#' @encoding UTF-8
#' @export
addBlockFromPlink2GDS <- function(gds, gdsOut, PATHBLOCK,
                                    superPop, blockName,
                                    blockDesc, verbose=FALSE) {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'")
    }

    ## The verbose must be a logical
    if (!(is.logical(verbose) && length(verbose) == 1)) {
        stop("The \'verbose\' parameter must be a logical (TRUE or FALSE).")
    }

    ## Extract the SNP chromosomes and positions
    snp.chromosome <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
    snp.position <- read.gdsn(index.gdsn(gds, "snp.position"))

    listChr <- unique(snp.chromosome)

    listChr <- listChr[order(listChr)]
    listChr <- seq_len(22)
    listBlock <- list()
    for(chr in listChr) {
        if(verbose) { message("chr", chr, " ",Sys.time()) }

        snp.keep <- snp.position[snp.chromosome == chr]

        listBlock[[chr]] <- processBlockChr(snp.keep, PATHBLOCK, superPop, chr)
        if(chr > 1) {
            vMax <- max(listBlock[[chr-1]])
            vMin <- min(listBlock[[chr-1]])
            listBlock[[chr]][listBlock[[chr]] > 0] <-
                    listBlock[[chr]][listBlock[[chr]] > 0] + vMax
            if(vMin < 0) {
                listBlock[[chr]][listBlock[[chr]] < 0] <-
                        listBlock[[chr]][listBlock[[chr]] < 0] + vMin
            }
        }
    }
    listBlock <- do.call(c, listBlock)

    ## Save the information into the GDS file
    addGDS1KGLDBlock(gdsOut, listBlock, blockName, blockDesc)

    ## Success
    return(0L)
}

#' @title Extract the specified column from the "sample.ref" node present in
#' the 1KG GDS file
#'
#' @description The function extract the specified column for the "sample.ref"
#' node present in the 1KG GDS file. The column must be present in the
#' \code{data.frame} saved in the "sample.ref" node.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param popName a \code{character} string representing the name of the column
#' that will be fetched in the \code{data.frame} present in the 1KG GDS
#' "sample.ref" node. The column must be present in the \code{data.frame}.
#'  Default: \code{"superPop"}.
#'
#' @return \code{data.frame} TODO study.annot with study.annot == study.id and
#' with the column popName.
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats rmultinom
#' @encoding UTF-8
#' @export
getRef1KGPop <- function(gds, popName="superPop") {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'")
    }

    ## The popName is a character string
    if (!is.character(popName)) {
        stop("The \'popName\' parameter must be a single character string.")
    }

    sample.ref <- read.gdsn(index.gdsn(gds, "sample.ref"))
    dataRef <- read.gdsn(index.gdsn(gds,
                            "sample.annot"))[which(sample.ref == TRUE),]

    if(! popName %in% colnames(dataRef)) {
        stop("The population ", popName, " is not supported ",
                "(not found in the 1KG GDS file).")
    }

    dataRef <- dataRef[, popName]
    names(dataRef) <- read.gdsn(index.gdsn(node=gds,
                            "sample.id"))[which(sample.ref == TRUE)]

    return(dataRef)
}

#' @title Generate two indexes based on gene annotation for gdsAnnot1KG block
#'
#' @description TODO
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param winSize a single positive \code{integer} representing the
#' size of the window to use to group the SNVs when the SNVs are in a
#' non-coding region. Default: \code{10000}.
#'
#' @param EnsDb An object with the ensembl genome annotation
#' Default: \code{EnsDb.Hsapiens.v86}.
#'
#' @return  a \code{data.frame} with those columns:
#' "chr", "pos", "snp.allele", "Exon", "GName", "Gene", "GeneS"
#' Example for GName and the two indexes "Gene", "GeneS"
#' GName Gene GeneS
#' 470                                 ENSG00000230021   17  3820
#' 471                                 ENSG00000230021   17  3820
#' 472                 ENSG00000230021:ENSG00000228794   17  3825
#' 473                 ENSG00000230021:ENSG00000228794   17  3825
#' 481 ENSG00000230021:ENSG00000228794:ENSG00000225880   17  3826
#' 482 ENSG00000230021:ENSG00000228794:ENSG00000225880   17  3826
#' 483 ENSG00000230021:ENSG00000228794:ENSG00000225880   17  3826
#' 492                 ENSG00000230021:ENSG00000228794   17  3825
#' 493                 ENSG00000230021:ENSG00000228794   17  3825
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom S4Vectors Rle
#' @importFrom BSgenome strand
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom IRanges IRanges
#' @importFrom AnnotationDbi select
#' @importFrom ensembldb exonsBy toSAF genes
#' @importFrom AnnotationFilter GeneIdFilter
#' @encoding UTF-8
#' @keywords internal
generateGeneBlock <- function(gds, winSize=10000, EnsDb) {

    edb <- EnsDb
    listEnsId <- unique(names(genes(edb)))

    cols <- c("GENEID", "SYMBOL", "GENENAME", "GENESEQSTART",
                "GENESEQEND", "SEQNAME")

    annot <- select(edb, keys=listEnsId, columns=cols, keytype="GENEID")
    annot <- annot[which(annot$SEQNAME %in% c(1:22, "X")),]

    # All the genes
    grGene <- GRanges(
        seqnames = annot$SEQNAME,
        ranges = IRanges(annot$GENESEQSTART, end = annot$GENESEQEND),
        strand = Rle(strand(rep("+", nrow(annot)))),
        mcols = annot[,c("GENEID", "GENENAME")])

    # Data frame of the all genes
    dfGenneAll <- as.data.frame(grGene)

    # group the overlapping gene
    grGeneReduce <- reduce(grGene)
    # data.frame version of grGeneReduce
    dfGene <- as.data.frame(grGeneReduce)

    # All exon
    allExon <- exonsBy(edb, by = "gene", filter = GeneIdFilter(listEnsId))
    # Transforming the GRangesList into a data.frame in SAF format
    dfExon <- toSAF(allExon)
    # remove the duplicate
    dfExon <- unique(dfExon)
    # Group the overlap
    exonReduce <- reduce(allExon)
    # Transforming the GRangesList into a data.frame in SAF format
    dfExonReduce <- toSAF(exonReduce)
    listMat <- list()

    matFreqAll <- data.frame(chr=read.gdsn(index.gdsn(gds, "snp.chromosome")),
                        pos=read.gdsn(index.gdsn(gds, "snp.position")),
                        snp.allele=read.gdsn(index.gdsn(gds, "snp.allele")),
                        stringsAsFactors=FALSE)
    offsetGene <- 0
    offsetGeneS <- 0
    offsetGene.O <- 0

    for(chr in seq_len(22))
    {
        dfExonChr <- dfExonReduce[which(dfExonReduce$Chr == chr),]
        dfGenneAllChr <- dfGenneAll[which(dfGenneAll$seqnames == chr),]
        dfGeneChr <- dfGene[which(dfGene$seqnames == chr),]
        # matFreq <- NULL
        #    matFreq <- read.csv2(fileSNV,
        #                         header=FALSE)


        # colnames(matFreq) <- c("chr", "pos", "ref", "alt", "af", "EAS_AF",
        #                        "EUR_AF","AFR_AF", "AMR_AF", "SAS_AF")
        print(system.time({
            # SNV in the GDS
            matFreq <- matFreqAll[which(matFreqAll$chr == chr),]
            # create two vector (one for the exon and one for the gene) of char
            # with 1 entry for each SNV in the GDS
            # I will keep the name of the gene and exon at this position
            listSNVExons <- character(nrow(matFreq))
            listSNVGenes <- character(nrow(matFreq))

            listPos <- seq_len(nrow(matFreq))
            listPos <- listPos[order(matFreq$pos)]
            # Create an index to accelerate the process
            startIndex <- seq(1, nrow(matFreq), 1000)
            # Add if the last entry is not the last position
            # is not the nb row of matFreq add the the last
            #position
            if(startIndex[length(startIndex)] < nrow(matFreq)){
                startIndex <- c(startIndex, nrow(matFreq))
            }
            # For gene in the chr
            # slow but acceptable
            #    user  system elapsed
            #    26.116   0.074  26.201
            # see blockAnnotation.R for slower alternatives
            for(genePos in seq_len(nrow(dfGenneAllChr))){
                # the gene is where SNV exists
                if(dfGenneAllChr$end[genePos] >= matFreq$pos[listPos[1]] &
                   dfGenneAllChr$start[genePos] <=  matFreq$pos[nrow(matFreq)]){
                    # In which partitions from the index the gene is located
                    vStart <- max(c(which(matFreq$pos[startIndex] <=
                                            dfGenneAllChr$start[genePos]), 1))
                    vEnd <- min(c(which(matFreq$pos[startIndex] >=
                                            dfGenneAllChr$end[genePos]),
                                            length(startIndex)))
                    # List of SNV in the gene
                    listP <- which(matFreq$pos[listPos[startIndex[vStart]:startIndex[vEnd]]] >= dfGenneAllChr$start[genePos] &
                                matFreq$pos[listPos[startIndex[vStart]:startIndex[vEnd]]] <= dfGenneAllChr$end[genePos])

                    # if SNV in the gene
                    if(length(listP) > 0 ) {
                        # listPos in the gene
                        listP <-
                            listPos[startIndex[vStart]:startIndex[vEnd]][listP]

                        # Add the name of the gene of SNVs
                        listSNVGenes[listP] <- paste0(listSNVGenes[listP], ":",
                            dfGenneAllChr$mcols.GENEID[genePos])

                        # Allow run on all without check if the SNV have
                        # already gene name
                        listSNVGenes[listP] <- gsub("^:", "",
                                                        listSNVGenes[listP])

                        # Exon of the gene
                        dfExon <- dfExonChr[which(dfExonChr$GeneID ==
                                        dfGenneAllChr$mcols.GENEID[genePos]),]
                        k <- 1

                        listE <- list()
                        for(pos in listP) {
                            if(length(which(dfExon$Start <= matFreq$pos[pos] &
                                        dfExon$End >= matFreq$pos[pos])) > 0) {
                                listE[[k]] <- pos
                                k <- k + 1
                            }
                        }

                        if(length(listE) > 0) {
                            listE <- do.call(c, listE)
                            listSNVExons[listE] <- paste0(listSNVExons[listE],
                                    ":", dfGenneAllChr$mcols.GENEID[genePos])
                            listSNVExons[listE] <- gsub("^:", "",
                                    listSNVExons[listE])
                        }
                    }
                }
            }
        }))


        # add the column Exon with the list of gene with an exon with the SNV
        matFreq$Exon <- listSNVExons
        # add the column GName with the list of gene with which include the SNV
        matFreq$GName <- listSNVGenes

        # dfGeneChr are reduced (merge all the overlap interval)
        z <- cbind(c(dfGeneChr$start, dfGeneChr$end, as.integer(matFreq$pos)),
                   c(seq_len(nrow(dfGeneChr)), -1 * seq_len(nrow(dfGeneChr)),
                                rep(0, nrow(matFreq))))
        z <- z[order(z[,1], -1 * z[,2]),]

        # group by interval which in overlap a gene
        matFreq$Gene[listPos] <- cumsum(z[,2])[z[,2] == 0]
        matFreq$Gene[matFreq$Gene > 0] <- matFreq$Gene[matFreq$Gene > 0] +
                                                                    offsetGene
        offsetGene <- max(offsetGene, max(matFreq$Gene))

        listD <- which(matFreq$Gene > 0)

        tmp <- paste0(matFreq[listD, "GName"], "_", matFreq[listD, "Gene"])
        listO <- order(tmp)


        # Create an index for each gene different if
        # two gene overlap the order don't have meaning.
        # gene ex: ENSG00000238009:ENSG00000239945_6 and ENSG00000238009_6
        # have a different number.
        # Note the order is base on the name not on position
        # Ex:
        #         GeneN                         indexNew
        # 135                 ENSG00000230021  4089
        # 136                 ENSG00000230021  4089
        # 148 ENSG00000230021:ENSG00000237973  4094
        # 149 ENSG00000230021:ENSG00000237973  4094
        # 159 ENSG00000229344:ENSG00000230021  4036
        # 160 ENSG00000229344:ENSG00000230021  4036
        # 161 ENSG00000230021:ENSG00000248527  4095
        # 162 ENSG00000198744:ENSG00000230021  3168
        # 163                 ENSG00000230021  4089
        # 164                 ENSG00000230021  4089
        # 165                 ENSG00000230021  4089
        #
        indexNew <- cumsum(!(duplicated(tmp[listO])))

        matFreq$GeneS <- rep(0, nrow(matFreq))
        matFreq$GeneS[listD][listO] <- indexNew + offsetGeneS
        offsetGeneS <- max(offsetGeneS, max(matFreq$GeneS))

        matFreq$GeneS[matFreq$GeneS < 0] <- 0
        matFreq$GeneS[matFreq$Gene < 0] <- 0
        listOrph <- which(matFreq$GeneS == 0)
        flag <- TRUE
        v <- offsetGene.O - 1
        i <- 1
        curZone <- "GeneS"
        curZone1 <- "Gene"
        winSize <- 10000

        if(length(listOrph) > 0){
            # Very slow can do better
            # but just run 1 time so less priority
            #    user  system elapsed
            # 517.595   7.035 524.658
            #    user  system elapsed
            # 558.526   2.274 561.043
            #
            while(flag){
                #use the index
                vStart <- min(c(which(matFreq$pos[startIndex] >
                                (matFreq[listOrph[i], "pos"] + winSize)),
                                length(startIndex)))

                preList <- listOrph[i]:startIndex[vStart]
                listWin <- which(matFreq[preList, "pos"] >
                                    (matFreq[listOrph[i], "pos"] + winSize) |
                                    (matFreq[preList, "pos"] >
                                        matFreq[listOrph[i], "pos"] &
                                        matFreq[preList,"GeneS"] > 0))

                j <- ifelse(length(listWin) > 0, preList[listWin[1]] - 1,
                                listOrph[i])

                matFreq[listOrph[i]:j, curZone] <- v
                matFreq[listOrph[i]:j, curZone1] <- v
                v <- v - 1
                i <- which(listOrph == j) + 1
                flag <- ifelse(i <= length(listOrph), TRUE, FALSE)

            }
            offsetGene.O <- min(offsetGene.O, min(matFreq$Gene))
        }

        listMat[[chr]] <- matFreq

        # save the matrix for each chr
        # create the space at the begining
    }

    matGene.Block <- do.call(rbind, listMat)
    rm(listMat)
    return(matGene.Block)
}


#' @title Generate two indexes based on gene annotation for gdsAnnot1KG
#' block and add the indexes into the
#' gdsAnnot1KG
#'
#' @description TODO
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param file.gdsRefAnnot the filename corresponding to an object of
#' class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file. The function will
#' open it in write mode and close it after.
#'
#' @param winSize a single positive \code{integer} representing the
#' size of the window to use to group the SNVs when the SNVs are in a
#' non-coding region. Default: \code{10000L}.
#'
#' @param EnsDb An object with the ensembl genome annotation
#' Default: \code{EnsDb.Hsapiens.v86}.
#'
#' @param suffixe.blockName TODO ex Ensembl.Hsapiens.v86
#'
#' @return \code{OL} when the function is successful.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt openfn.gds closefn.gds
#' @encoding UTF-8
#' @export
addGeneBlockGDSRefAnnot <- function(gds, file.gdsRefAnnot, winSize=10000,
                                            EnsDb, suffixe.blockName) {

    dfGeneBlock <- generateGeneBlock(gds, winSize, EnsDb)

    ## Opne GDS 1KG Annotation file in writting mode
    gdsRefAnnot <- openfn.gds(file.gdsRefAnnot, readonly=FALSE)


    blockName <- paste0("Gene.", suffixe.blockName)
    blockDesc <- paste0("List of blocks including overlapping genes ",
                                suffixe.blockName)
    addGDS1KGLDBlock(gdsRefAnnot, dfGeneBlock$Gene, blockName, blockDesc)
    blockName <- paste0("GeneS.", suffixe.blockName)
    blockDesc <- paste0("List of blocks of split by genes ", suffixe.blockName)
    addGDS1KGLDBlock(gdsRefAnnot, dfGeneBlock$GeneS, blockName, blockDesc)

    ## Close GDS 1KG annotation file
    closefn.gds(gdsRefAnnot)

    ## Success
    return(0L)
}


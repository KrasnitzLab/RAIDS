#' @title Prepare the pedigree file using pedigree information from Reference
#'
#' @description Using the pedigree file from Reference, this function extracts
#' needed information and formats it into a \code{data.frame} so in can
#' be used in following steps of the ancestry inference process. The
#' function also requires that the genotyping files associated to each
#' sample be available in a specified directory.
#'
#' @param filePed a \code{character} string representing the path and
#' file name of the pedigree file (PED file) that contains the information
#' related to the profiles present in the Reference GDS file. The PED file must
#' exist.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the Reference genotyping files for each profile are located. Only the
#' profiles with associated genotyping files are retained in the creation of
#' the final \code{data.frame}. The name of the genotyping files must
#' correspond to the individual identification (Individual.ID) in the
#' pedigree file (PED file).
#' Default: \code{"./data/sampleGeno"}.
#'
#' @param batch a\code{integer} that uniquely identifies the source of the
#' pedigree information. The Reference is usually \code{0L}.
#' Default: \code{0L}.
#'
#' @return a \code{data.frame} containing the needed pedigree information
#' from Reference. The \code{data.frame} contains those columns:
#' \itemize{
#' \item{sample.id}{a \code{character} string representing the profile unique
#' ID.}
#' \item{Name.ID}{a \code{character} string representing the profile name.}
#' \item{sex}{a \code{character} string representing the sex of the profile.}
#' \item{pop.group}{a \code{character} string representing the
#' sub-continental ancestry of the profile.}
#' \item{superPop }{a \code{character} string representing the continental
#' ancestry of the profile.}
#' \item{superPop }{a \code{integer} representing the batch of the profile.}
#' }
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo pedigree file
#' pedDemoFile <- file.path(dataDir, "PedigreeDemo.ped")
#'
#' ## Create a data.frame containing the information of the retained
#' ## samples (samples with existing genotyping files)
#' prepPed1KG(filePed=pedDemoFile, pathGeno=dataDir, batch=0L)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils read.delim
#' @encoding UTF-8
#' @export
prepPed1KG <- function(filePed, pathGeno=file.path("data", "sampleGeno"),
                        batch=0L) {

    ## Validate parameters
    validatePrepPed1KG(filePed=filePed, pathGeno=pathGeno, batch=batch)

    ## Read the pedigree file from Reference
    ped1KG <- read.delim(filePed)

    ## Create a data.frame containing the needed information
    pedAll <- data.frame(sample.id=c(ped1KG$Individual.ID),
                Name.ID=c(ped1KG$Individual.ID), sex=c(ped1KG$Gender),
                pop.group=c(ped1KG$Population),
                superPop=rep(NA, length(c(ped1KG$Population))),
                batch=c(rep(batch, nrow(ped1KG))),
                stringsAsFactors=FALSE)

    ## Create a list with all populations associated to each super-population
    ## NOTE The population versus super.population is hard-coded
    ## NOTE change to parameters
    listSuperPop1000G <- list()
    listSuperPop1000G[['EAS']] <- c("CHB", "JPT", "CHS", "CDX", "KHV")
    listSuperPop1000G[['EUR']] <- c("CEU", "TSI", "FIN", "GBR", "IBS")
    listSuperPop1000G[['AFR']] <- c("YRI", "LWK", "GWD", "MSL", "ESN",
                                        "ASW", "ACB")
    listSuperPop1000G[['AMR']] <- c("MXL", "PUR", "CLM", "PEL")
    listSuperPop1000G[['SAS']] <- c("GIH", "PJL", "BEB", "STU", "ITU")

    ## Identify the super-population associated to each sample in data.frame
    listSuperPop <- c("EAS", "EUR", "AFR", "AMR", "SAS")
    for(sp in listSuperPop){
        pedAll[which(pedAll$pop.group %in% listSuperPop1000G[[sp]]),
                "superPop"] <- sp
    }

    ## Assign row names to the data.frame
    row.names(pedAll) <- pedAll$sample.id

    ## Change column format for Sex information
    ## NOTE: could be done when the data.frame is created
    pedAll$sex <- as.character(pedAll$sex)

    ## Only retained samples with existing genotyping file
    listMat1k <- dir(pathGeno, pattern = ".+.csv.bz2")
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
#' file name of the bulk SNP information file from Reference. The file must
#' be in text format. The file must exist.
#'
#' @param fileSNPsRDS a \code{character} string representing the path and
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
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo SNV information file used as input
#' snvFile <- file.path(dataDir, "matFreqSNV_Demo.txt.bz2")
#'
#' ## Temporary output files
#' ## The first file contains the indexes of the retained SNPs
#' ## The second file contains the filtered SNP information
#' snpIndexFile <- file.path(getwd(), "listSNP_TEMP.rds")
#' filterSNVFile <- file.path(getwd(), "mapSNVSel_TEMP.rds")
#'
#' ## Run only if directory in writing mode
#' if (file.access(getwd()) == 0 && !file.exists(snpIndexFile) &&
#'         !file.exists(filterSNVFile)) {
#'
#'     ## Create a data.frame containing the information of the retained
#'     ## samples (samples with existing genotyping files)
#'     generateMapSnvSel(cutOff=0.01, fileSNV=snvFile,
#'         fileSNPsRDS=snpIndexFile, fileFREQ=filterSNVFile)
#'
#'     ## Remove temporary files
#'     unlink(snpIndexFile, force=TRUE)
#'     unlink(filterSNVFile, force=TRUE)
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
generateMapSnvSel <- function(cutOff=0.01, fileSNV, fileSNPsRDS, fileFREQ) {

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
    saveRDS(listSNP, fileSNPsRDS)
    ## Save the filtered SNP information file in RDS
    saveRDS(mapSNVSel, fileFREQ)

    return(0L)
}


#' @title Generate the GDS file that will contain the information from
#' Reference data set (reference data set)
#'
#' @description This function generates the GDS file that will contain the
#' information from Reference. The function also add the samples information,
#' the SNP information and the genotyping information into the GDS file.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the 1K genotyping files for each sample are located. The name of the
#' genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#' Default: \code{"./data/sampleGeno"}.
#'
#' @param filePedRDS a \code{character} string representing the path and file
#' name of the RDS file that contains the pedigree information. The file must
#' exist. The file must be a RDS file.
#'
#' @param fileSNVIndex a \code{character} string representing the path and file
#' name of the RDS file that contains the indexes of the retained SNPs. The
#' file must exist. The file must be a RDS file.
#'
#' @param fileSNVSelected a \code{character} string representing the path and
#' file name of the RDS file that contains the filtered SNP information. The
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
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The RDS file containing the pedigree information
#' pedigreeFile <- file.path(dataDir, "PedigreeDemo.rds")
#'
#' ## The RDS file containing the indexes of the retained SNPs
#' snpIndexFile <- file.path(dataDir, "listSNPIndexes_Demo.rds")
#'
#' ## The RDS file containing the filtered SNP information
#' filterSNVFile <- file.path(dataDir, "mapSNVSelected_Demo.rds")
#'
#' ## Temporary Reference GDS file
#' tempRefGDS <- file.path(getwd(), "1KG_TEMP.gds")
#'
#' ## Only run example if the directory is writable
#' if (file.access(getwd()) == 0 && !file.exists(tempRefGDS)) {
#'
#'     ## Different code depending of the withr package availability
#'     if (requireNamespace("withr", quietly=TRUE)) {
#'
#'         ## Temporary Reference GDS file
#'         gdsFile <- withr::local_file(tempRefGDS)
#'
#'         ## Create a temporary Reference GDS file containing
#'         ## information from reference file
#'         generateGDS1KG(pathGeno=dataDir, filePedRDS=pedigreeFile,
#'             fileSNVIndex=snpIndexFile, fileSNVSelected=filterSNVFile,
#'             fileNameGDS=gdsFile, listSamples=NULL)
#'
#'         ## Remove temporary files
#'         withr::deferred_run()
#'
#'     } else {
#'
#'         ## Create a temporary Reference GDS file
#'         generateGDS1KG(pathGeno=dataDir, filePedRDS=pedigreeFile,
#'             fileSNVIndex=snpIndexFile, fileSNVSelected=filterSNVFile,
#'             fileNameGDS=tempRefGDS, listSamples=NULL)
#'
#'         ## Remove temporary files
#'         unlink(tempRefGDS, force=TRUE)
#'
#'     }
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#'
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds
#' @encoding UTF-8
#' @export
generateGDS1KG <- function(pathGeno=file.path("data", "sampleGeno"),
                            filePedRDS, fileSNVIndex,
                            fileSNVSelected, fileNameGDS,
                            listSamples=NULL, verbose=FALSE) {

    validateGenerateGDS1KG(pathGeno=pathGeno, filePedRDS=filePedRDS,
        fileSNVIndex=fileSNVIndex, fileSNVSelected=fileSNVSelected,
        fileNameGDS=fileNameGDS, listSamples=listSamples, verbose=verbose)

    ## Read the pedigree file
    ped1KG <- readRDS(filePedRDS)

    # list in the file genotype we keep from fileSNPsRDS in generateMapSnvSel
    listKeep <- readRDS(fileSNVIndex)

    # Create the GDS file
    newGDS <- createfn.gds(fileNameGDS)
    put.attr.gdsn(newGDS$root, "FileFormat", "SNP_ARRAY")

    if(verbose) { message("Start ", Sys.time()) }

    listSampleGDS <- generateGDSRefSample(gdsReference=newGDS,
        dfPedReference=ped1KG, listSamples=listSamples)
    if(verbose) { message("Sample info DONE ", Sys.time()) }

    generateGDSSNPinfo(gdsReference=newGDS, fileFreq=fileSNVSelected,
                        verbose=verbose)
    if(verbose) { message("SNP info DONE ", Sys.time()) }

    generateGDSgenotype(gds=newGDS, pathGeno=pathGeno, fileSNPsRDS=fileSNVIndex,
        listSamples=listSampleGDS, verbose=verbose)
    if(verbose) { message("Genotype DONE ", Sys.time()) }

    closefn.gds(newGDS)
    return(0L)
}

#' @title Adding the phase information into the Reference GDS file
#'
#' @description The function is adding the phase information into the
#' Reference Phase GDS file. The phase information is extracted from a Reference
#' GDS file and is added into a Reference Phase GDS file. An entry called
#' 'phase' is added to the Reference Phase GDS file.
#'
#' @param gdsReference an object of class \link[gdsfmt]{gds.class} (GDS file),
#' an opened Reference GDS file.
#'
#' @param gdsReferencePhase an object of class \link[gdsfmt]{gds.class}
#' (GDS file), an opened Reference Phase GDS file.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the 1K genotyping files for each sample are located. The name of the
#' genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#' Default: \code{"./data/sampleGeno"}.
#'
#' @param fileSNPsRDS a \code{character} string representing the path and file
#' name of the RDS file that contains the indexes of the retained SNPs. The
#' file must exist. The file must be a RDS file.
#'
#' @param verbose a \code{logicial} indicating if the function should
#' print messages when running. Default: \code{FALSE}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Required package
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The RDS file containing the pedigree information
#' pedigreeFile <- file.path(dataDir, "PedigreeDemo.rds")
#'
#' ## The RDS file containing the indexes of the retained SNPs
#' snpIndexFile <- file.path(dataDir, "listSNPIndexes_Demo.rds")
#'
#' ## The RDS file containing the filtered SNP information
#' filterSNVFile <- file.path(dataDir, "mapSNVSelected_Demo.rds")
#'
#' ## Temporary Reference GDS file containing reference information
#' fileReferenceGDS <- "1KG_TEMP_02.gds"
#'
#' ## Only run example if the directory is writable
#' if (file.access(getwd()) == 0 && !file.exists(fileReferenceGDS)) {
#'
#'     ## Create a temporary Reference GDS file containing information from 1KG
#'     generateGDS1KG(pathGeno=dataDir, filePedRDS=pedigreeFile,
#'         fileSNVIndex=snpIndexFile, fileSNVSelected=filterSNVFile,
#'         fileNameGDS=fileReferenceGDS, listSamples=NULL)
#'
#'     ## Temporary Phase GDS file that will contain the 1KG Phase information
#'     fileRefPhaseGDS <- "1KG_TEMP_Phase_02.gds"
#'
#'     ## Create Reference Phase GDS file
#'     gdsPhase <- createfn.gds(fileRefPhaseGDS)
#'
#'     ## Open Reference GDS file
#'     gdsRef <- openfn.gds(fileReferenceGDS)
#'
#'     ## Fill temporary Reference Phase GDS file
#'     if (FALSE) {
#'         generatePhase1KG2GDS(gdsReference=gdsRef,
#'             gdsReferencePhase=gdsPhase,
#'             pathGeno=dataDir, fileSNPsRDS=filterSNVFile,
#'             verbose=FALSE)
#'     }
#'
#'     ## Close Reference Phase information file
#'     closefn.gds(gdsPhase)
#'
#'     ## Close Reference information file
#'     closefn.gds(gdsRef)
#'
#'     ## Remove temporary files
#'     unlink(fileReferenceGDS, force=TRUE)
#'     unlink(fileRefPhaseGDS, force=TRUE)
#'
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn readmode.gdsn
#' @encoding  UTF-8
#' @export
generatePhase1KG2GDS <- function(gdsReference, gdsReferencePhase,
                                    pathGeno, fileSNPsRDS, verbose=FALSE) {

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, "verbose")

    sample.id <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
    listSNP <- readRDS(fileSNPsRDS)

    var.phase <- NULL
    for(i in seq_len(length(sample.id))){

        if (verbose) { message("S ", i, " ", Sys.time()) }

        file1KG <- file.path(pathGeno, paste0(sample.id[i],".csv.bz2"))
        matSample <- read.csv2( file1KG,
                                row.names = NULL)[listSNP,, drop=FALSE]
        matSample <- matrix(as.numeric(unlist(strsplit(matSample[, 1],
                                                        "\\|"))), nrow=2)[1,]

        if (verbose) { message("GDS ", i, " ", Sys.time()) }

        if(! ("phase" %in% ls.gdsn(gdsReferencePhase))) {
            var.phase <- add.gdsn(gdsReferencePhase, "phase",
                                    valdim=c(length(listSNP), 1),
                                    matSample, storage="bit2",
                                    compress="LZ4_RA")
        }else {
            if(is.null(var.phase)) {
                var.phase <- index.gdsn(gdsReferencePhase, "phase")
            }
            append.gdsn(var.phase, matSample)
        }
        rm(matSample)
    }

    readmode.gdsn(var.phase)

    return(0L)
}

#' @title Identify genetically unrelated patients in GDS Reference file
#'
#' @description The function identify patients that are genetically related in
#' the Reference file. It generates a first RDS file with the list of unrelated
#' patient. It also generates a second RDS file with the kinship coefficient
#' between the patients.
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, the
#' Reference GDS file.
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
#' information about the Reference patients that are unrelated. The file will
#' contains two lists: the \code{list} of related samples, called \code{rels}
#' and the list of unrelated samples, called \code{unrels}.
#' The extension of the file must be '.rds'.
#'
#' @return \code{NULL} invisibly.
#'
#' @examples
#'
#' ## Required package
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Open existing demo Reference GDS file
#' fileGDS <- file.path(dataDir, "1KG_Demo.gds")
#' tmpGDS <- snpgdsOpen(fileGDS)
#'
#' ## Temporary output files
#' ## The first RDS file will contain the list of unrelated patients
#' ## The second RDS file will contain the kinship information between patients
#' patientTmpFile <-  "unrelatedPatients_TEMP.rds"
#' ibdTmpFile <- "ibd_TEMP.rds"
#'
#' ## Different code depending of the withr package availability
#' if (requireNamespace("withr", quietly=TRUE)) {
#'
#'     ## Temporary output files
#'     ## The first RDS file will contain the list of unrelated patients
#'     ## The second RDS file will contain the kinship information
#'     ## between patients
#'     patientTmpFileLocal <- withr::local_file(patientTmpFile)
#'     ibdTmpFileLocal <- withr::local_file(ibdTmpFile)
#'
#'     ## Identify unrelated patients in demo Reference GDS file
#'     identifyRelative(gds=tmpGDS, maf=0.05, thresh=2^(-11/2),
#'         fileIBD=ibdTmpFileLocal, filePart=patientTmpFileLocal)
#'
#'     ## Close demo Reference GDS file
#'     closefn.gds(tmpGDS)
#'
#'     ## Remove temporary files
#'     withr::deferred_run()
#'
#' } else {
#'
#'     ## Identify unrelated patients in demo Reference GDS file
#'     identifyRelative(gds=tmpGDS, maf=0.05, thresh=2^(-11/2),
#'         fileIBD=ibdTmpFile, filePart=patientTmpFile)
#'
#'     ## Close demo Reference GDS file
#'     closefn.gds(tmpGDS)
#'
#'     ## Remove temporary files
#'     unlink(patientTmpFile, force=TRUE)
#'     unlink(ibdTmpFile, force=TRUE)
#' }
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

    ibd.robust <- runIBDKING(gds=gds, maf=maf, verbose=FALSE)

    matKING <- ibd.robust$kinship
    colnames(matKING) <- ibd.robust$sample.id
    row.names(matKING) <- ibd.robust$sample.id

    part <- pcairPartition(kinobj=matKING, divobj=matKING,
                                kin.thresh=thresh, div.thresh=-1*thresh)

    saveRDS(ibd.robust, fileIBD)
    saveRDS(part, filePart)
}


#' @title Add the information about the unrelated patients to the Reference
#' GDS file
#'
#' @description This function adds the information about the unrelated patients
#' to the Reference GDS file. More specifically, it creates the field
#' \code{sample.ref} which as the value \code{1} when the sample
#' is unrelated and the value \code{0} otherwise.
#' The \code{sample.ref} is filled based on the information present in the
#' input RDS file.
#'
#' @param fileNameGDS a \code{character} string representing the path and file
#' name of the GDS file that contains the Reference information. The
#' Reference GDS file must contain the SNP information, the genotyping
#' information and the pedigree information from Reference dataset.
#' The extension of the file must be '.gds'.
#'
#' @param filePart a \code{character} string representing the path and file
#' name of the RDS file that contains the
#' information about the Reference patients that are unrelated.
#' The extension of the file must be '.rds'. The file must exists.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Locate RDS with unrelated/related status for 1KG samples
#' dataDir <- system.file("extdata", package="RAIDS")
#' rdsFilePath <- file.path(dataDir, "unrelatedPatientsInfo_Demo.rds")
#'
#' ## Create a temporary GDS file in an test directory
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(dataDir, "GDS_TEMP_201.gds")
#'
#' ## Create and open the GDS file
#' tmpGDS  <- createfn.gds(filename=gdsFilePath)
#
#' ## Create "sample.id" node (the node must be present)
#' sampleIDs <- c("HG00104", "HG00109", "HG00110")
#' add.gdsn(node=tmpGDS, name="sample.id", val=sampleIDs)
#'
#' ## Create "snp.id" node (the node must be present)
#' snpIDs <- c("s1", "s2", "s3", "s4", "s5", "s6")
#' add.gdsn(node=tmpGDS, name="snp.id", val=snpIDs)
#'
#' ## Create "snp.position" node (the node must be present)
#' snpPositions <- c(16102, 51478, 51897, 51927, 54489, 54707)
#' add.gdsn(node=tmpGDS, name="snp.position", val=snpPositions)
#'
#' ## Create "snp.chromosome" node (the node must be present)
#' snpPositions <- c(1, 1, 1, 1, 1, 1)
#' add.gdsn(node=tmpGDS, name="snp.chromosome", val=snpPositions)
#'
#' ## Create "genotype" node (the node must be present)
#' genotype <- matrix(rep(1, 18), ncol = 3)
#' add.gdsn(node=tmpGDS, name="genotype", val=genotype)
#'
#' ## Close GDS file
#' closefn.gds(tmpGDS)
#'
#' ## Create  "sample.ref" node in GDS file using RDS information
#' addRef2GDS1KG(fileNameGDS=gdsFilePath, filePart=rdsFilePath)
#'
#' ## Read sample reference data.frame
#' fileGDS <- openfn.gds(gdsFilePath, readonly=TRUE)
#' read.gdsn(index.gdsn(node=fileGDS, path="sample.ref"))
#' closefn.gds(gdsfile=fileGDS)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
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

    ## Open GDS file in write mode
    gdsReference <- snpgdsOpen(fileNameGDS, readonly=FALSE)

    ## Add information to the GDS file
    addGDSRef(gdsReference=gdsReference, filePart)

    ## Close GDS file
    closefn.gds(gdsReference)

    ## Success
    return(0L)
}


#' @title Compute principal component axes (PCA) on SNV data using the
#' reference samples
#'
#' @description The function runs a Principal Component Analysis (PCA) on
#' the SNv genotype data. The function also loads SNVs into the PCA to
#' calculate the SNV eigenvectors. Those 2 steps are done with the
#'  \code{\link[SNPRelate]{snpgdsPCA}} and
#'  \code{\link[SNPRelate]{snpgdsPCASNPLoading}}
#' functions.
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param listSample.Ref  a \code{vector} of \code{character} strings
#' corresponding to
#' the sample identifiers that will be used for the PCA.
#'
#' @param listSNP a \code{vector} of \code{character} strings representing
#' the SNV identifiers retained for the PCA.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @return a \code{list} with 3 entries:
#' \itemize{
#' \item{SNP}{ a \code{vector} of \code{character} strings representing the
#' SNV identifiers used in the PCA.}
#' \item{pca.unrel}{ an object of class \code{snpgdsPCAClass} as generated
#' by the
#' \code{\link[SNPRelate:snpgdsPCA]{SNPRelate::snpgdsPCA}} function. }
#' \item{snp.load}{ an object of class \code{snpgdsPCASNPLoading} as generated
#' by the
#' \code{\link[SNPRelate:snpgdsPCASNPLoading]{SNPRelate::snpgdsPCASNPLoading}}
#' function. }
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
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASNPLoading
#' @encoding UTF-8
#' @export
basePCASample <- function(gds, listSample.Ref=NULL, listSNP=NULL, np=1L) {

    listPCA <- list()

    ## Save the SNV list
    listPCA[["SNP"]] <- listSNP

    ## Calculate the PCA and save the results
    listPCA[["pca.unrel"]] <- snpgdsPCA(gds, sample.id=listSample.Ref,
                                            snp.id=listSNP,
                                            num.thread=np,
                                            verbose=TRUE)

    ## Calculate the SNV eigenvectors and save the results
    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                                    gdsobj=gds,
                                                    num.thread=np,
                                                    verbose=TRUE)

    ## Return a list with 3 entries
    return(listPCA)
}


#' @title Extract the specified column from the 1KG GDS 'sample.ref' node
#' for the reference profiles (real ancestry assignation)
#'
#' @description The function extract the specified column for the 'sample.ref'
#' node present in the Reference GDS file. The column must be present in the
#' \code{data.frame} saved in the 'sample.ref' node. Only the information for
#' the reference profiles is returned. The values
#' represent the known ancestry assignation.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Reference GDS file.
#'
#' @param popName a \code{character} string representing the name of the column
#' that will be fetched in the \code{data.frame} present in the Reference GDS
#' "sample.ref" node. The column must be present in the \code{data.frame}.
#'  Default: \code{"superPop"}.
#'
#' @return \code{vector} of \code{character} strings representing the content
#' of the extracted column for the 1KG GDS 'sample.ref' node. The values
#' represent the known ancestry assignation. The profile
#' identifiers are used as names for the \code{vector}.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Open existing demo 1K GDS file with "sample.ref" node
#' nameFileGDS <- file.path(dataDir, "1KG_Demo_with_sampleREF.gds")
#' fileGDS <- snpgdsOpen(nameFileGDS)
#'
#' ## Extract super population information for the 1KG profiles
#' getRef1KGPop(gdsReference=fileGDS, popName="superPop")
#'
#' ## Close 1K GDS file
#' closefn.gds(fileGDS)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats rmultinom
#' @encoding UTF-8
#' @export
getRef1KGPop <- function(gdsReference, popName="superPop") {

    ## The gdsReference must be an object of class "gds.class"
    validateGDSClass(gds=gdsReference, "gdsReference")

    ## The popName is a character string
    if (!is.character(popName)) {
        stop("The \'popName\' parameter must be a single character string.")
    }

    sample.ref <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    dataRef <- read.gdsn(index.gdsn(gdsReference,
                            "sample.annot"))[which(sample.ref == TRUE),]

    if(! popName %in% colnames(dataRef)) {
        stop("The population ", popName, " is not supported ",
                "(not found in the 1KG GDS file).")
    }

    dataRef <- dataRef[, popName]
    names(dataRef) <- read.gdsn(index.gdsn(node=gdsReference,
                            "sample.id"))[which(sample.ref == TRUE)]

    return(dataRef)
}


#' @title Append information associated to blocks, as indexes, into the
#' Population Reference SNV Annotation GDS file
#'
#' @description The function appends the information about the blocks into
#' the Population Reference SNV Annotation GDS file. The information is
#' extracted from the Population Reference GDS file.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Reference GDS file.
#'
#' @param gdsRefAnnotFile a \code{character} string representing the
#' file name corresponding the Reference SNV
#' Annotation GDS file. The function will
#' open it in write mode and close it after. The file must exist.
#'
#' @param winSize a single positive \code{integer} representing the
#' size of the window to use to group the SNVs when the SNVs are in a
#' non-coding region. Default: \code{10000L}.
#'
#' @param ensDb An object with the ensembl genome annotation
#' Default: \code{EnsDb.Hsapiens.v86}.
#'
#' @param suffixBlockName a \code{character} string that identify the source
#' of the block and that will be added to the block description into
#' the Reference SNV Annotation GDS file, as example: Ensembl.Hsapiens.v86.
#'
#' @return The integer \code{OL} when the function is successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt openfn.gds closefn.gds
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
addGeneBlockGDSRefAnnot <- function(gdsReference, gdsRefAnnotFile,
                        winSize=10000, ensDb, suffixBlockName) {

    ## The gdsReference must be an object of class "gds.class"
    if (!inherits(gdsReference, "gds.class")) {
        stop("The \'gdsReference\' must be an object of class \'gds.class\'")
    }

    ## Validate that the Reference Annotation GDS file exists
    if (! file.exists(gdsRefAnnotFile)) {
        stop("The file \'", gdsRefAnnotFile, "\' does not exist.")
    }

    ## The winSize must be a positive single number
    if (!(isSingleNumber(winSize) && (winSize > 0))) {
        stop("The \'winSize\' parameter must be a single numeric value." )
    }

    ## Generate two indexes based on gene annotation for
    ## the Reference GDS Annotation block
    dfGeneBlock <- generateGeneBlock(gdsReference=gdsReference,
                                     winSize=winSize, ensDb=ensDb)

    ## Open GDS Reference Annotation file in writing mode
    gdsRefAnnot <- openfn.gds(gdsRefAnnotFile, readonly=FALSE)

    blockName <- paste0("Gene.", suffixBlockName)
    blockDesc <- paste0("List of blocks including overlapping genes ",
                                suffixBlockName)
    addBlockInGDSAnnot(gdsRefAnnot, dfGeneBlock$Gene, blockName, blockDesc)
    blockName <- paste0("GeneS.", suffixBlockName)
    blockDesc <- paste0("List of blocks of split by genes ", suffixBlockName)
    addBlockInGDSAnnot(gdsRefAnnot, dfGeneBlock$GeneS, blockName, blockDesc)

    ## Close GDS Reference annotation file
    closefn.gds(gdsRefAnnot)

    ## Success
    return(0L)
}


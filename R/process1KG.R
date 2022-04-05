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
#' data.dir <- system.file("extdata", package="aicsPaper")
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
#' @importFrom S4Vectors isSingleInteger
#' @encoding UTF-8
#' @export
prepPed1KG <- function(pedFile, PATHGENO=file.path("data", "sampleGeno"),
                        batch.v=0L) {

    ## Validate that the batch is an integer
    if (! isSingleInteger(batch.v)) {
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
#' data.dir <- system.file("extdata", package="aicsPaper")
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
#' @return None.
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
#' data.dir <- system.file("extdata", package="aicsPaper")
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
                            listSamples=NULL) {

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


    print(paste0("Start ", Sys.time()))

    listSampleGDS <- generateGDSSample(newGDS, ped1KG, listSamples)
    print(paste0("Sample info DONE ", Sys.time()))

    generateGDSSNPinfo(newGDS, fileSNPSel)
    print(paste0("SNP info DONE ", Sys.time()))


    generateGDSgenotype(newGDS, PATHGENO, fileListSNP, listSampleGDS)
    print(paste0("Genotype DONE ", Sys.time()))

    closefn.gds(newGDS)
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
#' @return TODO a \code{vector} of \code{string}
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn readmode.gdsn
#' @encoding  UTF-8
#' @export
generatePhase1KG2GDS <- function(gds, gdsPhase,
                            PATHGENO,
                            fileLSNP){

    sample.id <- read.gdsn(index.gdsn(gds,"sample.id"))
    listSNP <- readRDS(fileLSNP)


    var.phase <- NULL
    for(i in seq_len(length(sample.id))){

        print(paste0("S ", i, " ", Sys.time()))

        file1KG <- file.path(PATHGENO, paste0(sample.id[i],".csv.bz2"))
        matSample <- read.csv2( file1KG,
                                row.names = NULL)[listSNP,, drop=FALSE]
        matSample <- matrix(as.numeric(unlist(strsplit( matSample[,1], "\\|"))),nrow=2)[1,]
        print(paste0("GDS ", i, " ", Sys.time()))
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
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
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
#'
#' @param filePart a \code{character} string representing the path and file
#' name of the RDS file that will be created. The RDS file will contain the
#' information about the 1KG patients that are unrelated.
#' The extension of the file must be '.rds'.
#'
#' @return \code{NULL} invisibly.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
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
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#'
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
#' @param gds an object of class \code{gds} opened
#'
#' @param method a \code{character string} TODO . Default: \code{corr}
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
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @export
pruning1KG.Chr <- function(gds,
                            method="corr",
                            listSamples=NULL,
                            slide.max.bp.v=5e5,
                            ld.threshold.v=sqrt(0.1),
                            np=1,
                            verbose.v=FALSE,
                            chr=NULL,
                            minAF = NULL,
                            outPref = "pruned_1KG",
                            keepObj = FALSE) {

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

            listKeep <- snpID[which(snpCHR == chr & snpAF >= minAF & snpAF <= 1-minAF)]
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




#' @title Create a study at GDS including the reference (first study add)
#'
#' @description TODO
#'
#' @param PATHGENO TODO a PATH to the directory genotype file of 1KG
#' The directory sampleGeno must contain matFreqSNV.txt.bz2
#'
#' @param fileNamePED TODO
#'
#' @param fileNameGDS TODO
#'
#' @param batch TODO
#'
#' @param studyDF TODO
#'
#' @param listSamples A \code{vector} of \code{string} corresponding to
#' the sample.ids. if NULL all samples
#'
#' @param KEEPCOV TODO a \code{logical} if it is keeping the coverage
#'
#' @param PATHSAMPLEGDS TODO a PATH to a directory where a gds specific
#' to the samples with coverage info is keep
#'
#' @return None
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
#'
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds read.gdsn
#' @export
appendStudy2GDS1KG <- function(PATHGENO=file.path("data", "sampleGeno"),
                                fileNamePED,
                                fileNameGDS,
                                batch=1,
                                studyDF,
                                listSamples=NULL,
                                KEEPCOV=TRUE,
                                PATHSAMPLEGDS=NULL){

    # check if file fileGDS
    # It must not exists

    # validate the para

    pedStudy <- readRDS(fileNamePED)


    # list in the file genotype we keep from fileLSNP in generateMapSnvSel




    # Create the GDS file
    gds <- snpgdsOpen(fileNameGDS, readonly=FALSE)

    snpCHR <- index.gdsn(gds, "snp.chromosome")
    snpPOS <- index.gdsn(gds, "snp.position")

    listPos <- data.frame(snp.chromosome=read.gdsn(snpCHR),
                            snp.position=read.gdsn(snpPOS))



    print(paste0("Start ", Sys.time()))

    listSampleGDS <- addStudyGDSSample(gds, pedDF=pedStudy, batch=batch,
                                        listSamples=listSamples,
                                        studyDF=studyDF)


    print(paste0("Sample info DONE ", Sys.time()))

    generateGDS1KGgenotypeFromSNPPileup(gds, PATHGENO,
                                            listSamples=listSamples,
                                            listPos=listPos, offset=-1,
                                            minCov=10, minProb=0.999,
                                            seqError=0.001,
                                            KEEPCOV = TRUE,
                                            PATHGDSSAMPLE=PATHSAMPLEGDS)

    print(paste0("Genotype DONE ", Sys.time()))

    closefn.gds(gds)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened
#'
#' @param listSample.Ref  A \code{vector} of \code{string} corresponding to
#' the sample.ids
#'
#' @param listSNP the list of snp.id keep
#'
#' @param np a \code{integer} reprenting the number of thread
#'
#'
#' @return TODO a \code{list}  with with two objects pca.unrel -> \code{snpgdsPCAClass}
#' and a snp.load -> \code{snpgdsPCASNPLoading}
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
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



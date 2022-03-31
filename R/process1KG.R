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
#' \code{data.frame}. Default: \code{"./data/sampleGeno"}.
#'
#' @param batch.v a\code{integer} that uniquely identifies the source of the
#' pedigree information. The 1KG is usually \code{0L}. Default: \code{0L}.
#'
#' @return a \code{data.frame} containing the needed pedigree information
#' from 1K.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#' @importFrom utils read.delim
#' @importFrom S4Vectors isSingleInteger
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
    if (! file.exists(pedFile)) {
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


#' @title Generate the mapSNVSel file
#'
#' @description TODO
#'
#' @param cutOff a single \code{numeric} value, the cut off
#' for the frequency in at least one super population.
#'
#' @param fileSNV TODO
#'
#' @param fileLSNP TODO
#'
#' @param fileFREQ TODO
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#'
#' @export

generateMapSnvSel <- function(cutOff = 0.01, fileSNV, fileLSNP, fileFREQ){

    # test fileSNV

    mapSNVSel <- read.csv2(fileSNV)
    listSNP <- which( rowSums( mapSNVSel[,c("EAS_AF",
                                            "EUR_AF",
                                            "AFR_AF",
                                            "AMR_AF",
                                            "SAS_AF")] >= (cutOff - 1e-10) &
                                   mapSNVSel[,c("EAS_AF",
                                                "EUR_AF",
                                                "AFR_AF",
                                                "AMR_AF",
                                                "SAS_AF")] <= (1 - cutOff + 1e-10) ) >0)

    mapSNVSel <- mapSNVSel[listSNP,]

    saveRDS(listSNP, fileLSNP)
    saveRDS(mapSNVSel, fileFREQ)

    return(0L)
}


#' @title Generate the base GDS file from 1KG
#'
#' @description TODO
#'
#' @param PATHGENO TODO a PATH to the directory genotype file of 1KG
#' The directory sampleGeno must contain matFreqSNV.txt.bz2
#'
#' @param fileNamePED TODO
#'
#' @param fileListSNP TODO
#'
#' @param fileSNPSel TODO
#'
#' @param fileNameGDS TODO
#'
#' @param listSamples a \code{vector} of \code{string} corresponding to
#' the sample.ids
#' if NULL all the samples
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#'
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds
#'
#' @export
generateGDS1KG <- function(PATHGENO = file.path("data", "sampleGeno"),
                            fileNamePED, fileListSNP,
                            fileSNPSel, fileNameGDS,
                            listSamples=NULL) {

    # check if file fileGDS
    # It must not exists

    # validate the para

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


#' @title Extract a list of id unrelated and a list related
#'
#' @description TODO
#'
#' @param gds TODO
#'
#' @param maf TODO
#'
#' @param thresh TODO
#'
#' @param fileIBD  TODO
#'
#' @param filePart TODO
#'
#' @return \code{NULL} invisibly.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#'
#' @importFrom GENESIS pcairPartition
#'
#' @export
identifyRelative <- function(gds, maf=0.05, thresh=2^(-11/2),
                             fileIBD, filePart) {

    ibd.robust <- runIBDKING(gds=gds, maf=maf)

    matKING <- ibd.robust$kinship
    colnames(matKING) <- ibd.robust$sample.id
    row.names(matKING) <- ibd.robust$sample.id

    part <- pcairPartition(kinobj=matKING, divobj=matKING,
                                kin.thresh=thresh, div.thresh=-1*thresh)

    saveRDS(ibd.robust, fileIBD)
    saveRDS(part, filePart)

}


#' @title This function create the field sample.ref which is 1 when the sample
#' are a reference and 0 otherwise. The sample.ref is fill base of on the file filePart
#'
#' @description This function create the field sample.ref which is 1 when the sample
#' are a reference and 0 otherwise. The sample.ref is fill base of on the file filePart$unrels
#' from  in GENESIS TODO
#'
#' @param fileNameGDS  \code{string} with the path and the name of the gds file
#'
#' @param filePart file save by identifyRelative
#'
#' @return None.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#'
#' @importFrom SNPRelate snpgdsOpen
#'
#' @export
addRef2GDS1KG <- function(fileNameGDS,
                             filePart){



    gds <- snpgdsOpen(fileNameGDS, readonly = FALSE)

    addGDSRef(gds, filePart)

    closefn.gds(gds)
}


#' @title Find the pruned snv in 1KG by chr
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened
#'
#' @param method TODO
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
#' # TODO
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
                            keepObj = FALSE){

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
#' @param PATHGDSSAMPLE TODO a PATH to a directory where a gds specific
#' to the samples with coverage info is keep
#'
#' @return None
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#'
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds
#'
#' @export


appendStudy2GDS1KG <- function(PATHGENO = file.path("data", "sampleGeno"),
                               fileNamePED,
                               fileNameGDS,
                               batch = 1,
                               studyDF,
                               listSamples = NULL,
                               KEEPCOV = TRUE,
                               PATHSAMPLEGDS = NULL){

    # check if file fileGDS
    # It must not exists

    # validate the para

    pedStudy <- readRDS(fileNamePED)


    # list in the file genotype we keep from fileLSNP in generateMapSnvSel




    # Create the GDS file
    gds <- snpgdsOpen(fileNameGDS, readonly=FALSE)

    snpCHR <- index.gdsn(gds, "snp.chromosome")
    snpPOS <- index.gdsn(gds, "snp.position")

    listPos <- data.frame(snp.chromosome = read.gdsn(snpCHR),
                          snp.position = read.gdsn(snpPOS))



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

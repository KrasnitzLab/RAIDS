#' @title Prepare the ped file for 1KG
#'
#' @description TODO
#'
#' @param pedFile TODO a \code{string} path and file name of ped file from 1KG
#'
#' @param PATHGENO TODO
#'
#' @param batch.v TODO
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#' @importFrom utils read.delim
#' @export
prepPed1KG <- function(pedFile, PATHGENO=file.path("data", "sampleGeno"),
                        batch.v=0) {

    # TODO validate

    ped1KG <- read.delim(pedFile)

    pedAll <- data.frame(
            sample.id=c(ped1KG$Individual.ID),
            Name.ID=c(ped1KG$Individual.ID),
            sex=c(ped1KG$Gender),
            pop.group=c(ped1KG$Population),
            superPop=rep(NA, length(c(ped1KG$Population))),
            batch=c(rep(batch.v,nrow(ped1KG))),
            stringsAsFactors=FALSE
    )


    # TODO The population versus super.population is hardcode
    # change to parameters
    listSuperPop1000G <- list()
    listSuperPop1000G[['EAS']] <- c("CHB", "JPT", "CHS", "CDX", "KHV")
    listSuperPop1000G[['EUR']] <- c("CEU", "TSI", "FIN", "GBR", "IBS")
    listSuperPop1000G[['AFR']] <- c("YRI", "LWK", "GWD", "MSL", "ESN", "ASW", "ACB")
    listSuperPop1000G[['AMR']] <- c("MXL", "PUR", "CLM", "PEL")
    listSuperPop1000G[['SAS']] <- c("GIH", "PJL", "BEB", "STU", "ITU")
    listSuperPop <- c("EAS", "EUR", "AFR", "AMR", "SAS")

    for(sp in listSuperPop){
        pedAll[which(pedAll$pop.group %in% listSuperPop1000G[[sp]]), "superPop"] <-
            sp
    }
    row.names(pedAll) <- pedAll$sample.id
    pedAll$sex <- as.character(pedAll$sex)


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
#' @param fileSNV
#'
#' @param fileLSNP
#'
#' @param fileFREQ
#'
#' @return NULL
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
#'
#' @export

generateMapSnvSel <- function(cutOff = 0.01, fileSNV, fileLSNP, fileFREQ){

    # test fileSNV

    mapSNVSel <- read.csv2(fileSNV)
    listSNP <- which(rowSums(mapSNVSel[,c("EAS_AF",
                                          "EUR_AF",
                                          "AFR_AF",
                                          "AMR_AF",
                                          "SAS_AF")] >= cutOff &
                                 mapSNVSel[,c("EAS_AF",
                                              "EUR_AF",
                                              "AFR_AF",
                                              "AMR_AF",
                                              "SAS_AF")] <= 1 - cutOff)>0)

    mapSNVSel <- mapSNVSel[listSNP,]

    saveRDS(listSNP, fileLSNP)
    saveRDS(mapSNVSel, fileFREQ)
}


#' @title Generate the base GDS file from 1KG
#'
#' @description TODO
#'
#' @param PATHGENO TODO a PATH to the directory genotype file of 1KG
#' The directory sampleGeno must contain matFreqSNV.txt.bz2
#'
#' @param fileNamePED
#'
#' @param fileListSNP
#'
#' @param fileSNPSel
#'
#' @param fileNameGDS
#'
#' @param listSamples A \code{vector} of \code{string} corresponding to the sample.ids
#' if NULL all the samples
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
#'
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds
#'
#' @export




generateGDS1KG <- function(PATHGENO = file.path("data", "sampleGeno"),
                           fileNamePED,
                           fileListSNP,
                           fileSNPSel,
                           fileNameGDS,
                           listSamples = NULL){

    # check if file fileGDS
    # It must not exists

    # validate the para

    ped1KG <- readRDS(fileNamePED)


    # list in the file genotype we keep from fileLSNP in generateMapSnvSel
    listKeep <- readRDS(fileListSNP)



    # Create the GDS file
    newGDS <- createfn.gds(fileGDS)
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
#' @param gds
#'
#' @param maf
#'
#' @param thresh
#'
#' @param fileIBD
#'
#' @param filePart
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
#'
#' @importFrom GENESIS pcairPartition
#'
#' @export

identifyRelative <- function(gds,
                             maf = 0.05,
                             thresh = 2^(-11/2),
                             fileIBD,
                             filePart){

    ibd.robust <- runIBDKING(gds=gds, maf=maf)

    matKING <- ibd.robust$kinship
    colnames(matKING) <- ibd.robust$sample.id
    row.names(matKING) <- ibd.robust$sample.id

    part <- pcairPartition(kinobj = matKING, divobj = matKING,
                                kin.thresh = thresh, div.thresh = -1 *thresh)

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
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
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
#' @param method
#'
#' @param listSamples
#'
#' @param slide.max.bp.v
#'
#' @param ld.threshold.v
#'
#' @param np
#'
#' @param verbose.v
#'
#' @param chr
#'
#' @param minAF
#'
#' @param outPref
#'
#' @param keepObj
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @export

pruning1KG.Chr <- function(gds,
                           method="corr",
                           listSamples=NULL,
                           slide.max.bp.v = 5e5,
                           ld.threshold.v=sqrt(0.1),
                           np = 1,
                           verbose.v=FALSE,
                           chr = NULL,
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
#' @param fileNamePED
#'
#' @param fileNameGDS
#'
#' @param batch
#'
#' @param studyDesc
#'
#' @param listSamples A \code{vector} of \code{string} corresponding to the sample.ids
#' if NULL all samples
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
#'
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds
#'
#' @export


appendStudy2GDS1KG <- function(PATHGENO = file.path("data", "sampleGeno"),
                               fileNamePED,
                               fileNameGDS,
                               batch = 1,
                               studyDesc,
                               listSamples = NULL){

    # check if file fileGDS
    # It must not exists

    # validate the para

    pedStudy <- readRDS(fileNamePED)


    # list in the file genotype we keep from fileLSNP in generateMapSnvSel




    # Create the GDS file
    gds <- snpgdsOpen(fileNameGDS, readonly = FALSE)

    snpCHR <- index.gdsn(gds, "snp.chromosome")
    snpPOS <- index.gdsn(gds, "snp.position")

    listPos <- data.frame(snp.chromosome = read.gdsn(snpCHR),
                          snp.position = read.gdsn(snpPOS))



    print(paste0("Start ", Sys.time()))

    listSampleGDS <- addStudyGDSSample(gds, pedDF=pedStudy, batch=batch, listSamples = listSamples, studyDF = studyDF)


    print(paste0("Sample info DONE ", Sys.time()))

    generateGDS1KGgenotypeFromSNPPileup(gds, PATHGENO,
                                        listSamples = listSamples, listPos = listPos, offset = -1,
                                        minCov = 10, minProb = 0.999, seqError = 0.001)
    print(paste0("Genotype DONE ", Sys.time()))

    closefn.gds(gds)
}

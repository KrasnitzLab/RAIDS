#' @title Compute the list of pruned SNVs for a specific profile using the
#' information from the Reference GDS file and a linkage disequilibrium
#' analysis
#'
#' @description This function computes the list of pruned SNVs for a
#' specific profile. When
#' a group of SNVs are in linkage disequilibrium, only one SNV from that group
#' is retained. The linkage disequilibrium is calculated with the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The initial list of
#' SNVs that are passed to the \code{\link[SNPRelate]{snpgdsLDpruning}}()
#' function can be specified by the user.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return The function returns a \code{RAIDSparam} an object with all the RAIDS
#' parameters updated.
#'
#' @examples
#'
#' ## Path to the demo Reference GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG.gds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyDF <- data.frame(study.id = "MYDATA",
#'                         study.desc = "Description",
#'                         study.platform = "PLATFORM",
#'                         stringsAsFactors = FALSE)
#'
#' ## The data.frame containing the information about the samples
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' samplePED <- data.frame(Name.ID = c("ex1", "ex2"),
#'                     Case.ID = c("Patient_h11", "Patient_h12"),
#'                     Diagnosis = rep("Cancer", 2),
#'                     Sample.Type = rep("Primary Tumor", 2),
#'                     Source = rep("Databank B", 2), stringsAsFactors = FALSE)
#' rownames(samplePED) <- samplePED$Name.ID
#'
#' ## Temporary Profile GDS file
#' pathProfileGDS <- tempdir()
#' profileFile <- file.path(pathProfileGDS, "ex1.gds")
#'
#' pRAIDS <- RAIDSparam( studyDF=studyDF,
#'                       profileFile=file.path(dataDir, "ex1.txt.gz"),
#'                       genoSource = "snp-pileup",
#'                       pedStudy=samplePED,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileGDS)
#'
#' ## Copy the Profile GDS file demo that has not been pruned yet
#' file.copy(file.path(dataDir, "ex1_demo.gds"), profileFile)
#'
#'
#' ## Compute the list of pruned SNVs for a specific profile 'ex1'
#' ## and save it in the Profile GDS file 'ex1.gds'
#' pRAIDS <- pruningProfile2(pRAIDS)
#'
#'
#' ## Check content of Profile GDS file
#' ## The 'pruned.study' entry should be present
#' content <- openfn.gds(profileFile)
#' content
#'
#' ## Close the Profile GDS file (important)
#' closefn.gds(content)
#'
#' ## Remove Profile GDS file (created for demo purpose)
#' unlink(profileFile, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn openfn.gds ls.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
pruningProfile2 <- function(pRAIDS) {

    ## Matches a character method against a table of candidate values
    # method <- arg_match(method)
    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))

    ## Profile GDS file name
    fileGDSSample <-  validateProfileGDSExist(pathProfile=pathProfileGDS(pRAIDS),
                                              profile=pedStudy(pRAIDS)$Name.ID[1])

    # filePruned <- file.path(pathPrunedGDS, paste0(outPrefix, ".rds"))
    # fileObj <- file.path(pathPrunedGDS, paste0(outPrefix, ".Obj.rds"))

    snp.id <- read.gdsn(node=index.gdsn(gdsReference, "snp.id"))
    sample.id <- read.gdsn(node=index.gdsn(gdsReference, "sample.id"))
    keepPos <- seq_len(length(snp.id))
    if("snp.KeepDefault" %in% ls.gdsn(gdsReference) ){
        keepPos <- read.gdsn(index.gdsn(gdsReference, "snp.KeepDefault"))
    }
    snp.id <- snp.id[keepPos]
    sample.ref <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    listSamples <- sample.id[which(sample.ref == 1)]
    snpgdsClose(gdsReference)
    ## Open the GDS Sample file
    gdsSample <- openfn.gds(filename=fileGDSSample)

    ## Extract all study information from the GDS Sample file
    study.annot <- read.gdsn(node=index.gdsn(gdsSample, "study.annot"))

    ## Select study information associated to the current profile
    posSample <- which(study.annot$data.id == pedStudy(pRAIDS)$Name.ID[1] &
                           study.annot$study.id == studyDF(pRAIDS)$study.id)

    ## Check that the information is found for the specified profile and study
    if(length(posSample) != 1) {
        closefn.gds(gdsSample)
        stop("In pruningSample the profile \'", pedStudy(pRAIDS)$Name.ID[1],
             "\' doesn't exists for the study \'", studyDF(pRAIDS)$study.id[1], "\'\n")
    }

    ## Get the SNV genotype information for the current profile
    # "geno.ref"
    # TODO use genoType(pRAIDS) when it is ready
    if(genoType(pRAIDS) %in% ls.gdsn(gdsSample)){
        g <- read.gdsn(index.gdsn(gdsSample, genoType(pRAIDS)),
                       start=c(1, posSample), count=c(-1,1))

    }else{
        closefn.gds(gdsSample)
        stop("The genotype type don't exists")
    }

    ## Close the Profile GDS file
    closefn.gds(gdsSample)

    listGeno <- which(g != 3)
    rm(g)
    # This position are the position in keepPos (I suppose the position in genotype match keepPos)
    listKeepPos <- listGeno

    if(!is.null(specificSNV(pRAIDS))){
        # we need to intersect
        if(!("snvKeep" %in% colnames(specificSNV(pRAIDS)))){
            pRAIDS <- specificSNVColKeep(pRAIDS)
        }
        listKeepPos <- intersect(specificSNV(pRAIDS)$snvKeep[specificSNV(pRAIDS)$snvKeep > 0], listKeepPos)
    }



    if (length(listKeepPos) == 0) {
        stop("In pruningSample, the sample ", pedStudy(pRAIDS)$Name.ID[1],
             " doesn't have SNPs after filters\n")
    }
    # Not snp.id match keepPos
    listKeep <- snp.id[listKeepPos]


    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))
    on.exit(snpgdsClose(gdsReference), add=TRUE)
    ## Use a LD analysis to generate a subset of SNPs
    snpset <- runLDPruning(gds=gdsReference, method=pruningMethod(pRAIDS),
                           listSamples=listSamples, listKeep=listKeep,
                           slideWindowMaxBP=slideWindowMaxBP(pRAIDS), thresholdLD=thresholdLD(pRAIDS),
                           np=np(pRAIDS), verbose=verbose(pRAIDS))
    #snpgdsClose(gdsReference)
    pruned <- unlist(snpset, use.names=FALSE)

    gdsProfile <- openfn.gds(filename=fileGDSSample, readonly=FALSE)
    addGDSStudyPruning(gdsProfile=gdsProfile, pruned=pruned)
    closefn.gds(gdsfile=gdsProfile)


    return(pRAIDS)
}

#' @title Run most steps leading to the ancestry inference call based on PCA
#' and KNN on a specific profile
#'
#' @description This function runs most steps leading to the ancestry inference
#' call on a specific profile. First, the function creates the
#' Profile GDS file for the specific profile using the information from a
#' RDS Sample description file and the Population Reference GDS file.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{pcaSample}}{ a \code{list} containing the information related
#' to the eigenvectors. The \code{list} contains those 3 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{eigenvector.ref}}{ a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{\code{eigenvector}}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current profile projected on the PCA from the
#' reference profiles.}
#' }
#' }
#' \item{\code{paraSample}}{ a \code{list} containing the results with
#' different \code{D} and \code{K} values that lead to optimal parameter
#' selection. The \code{list} contains those entries:
#' \describe{
#' \item{\code{dfPCA}}{ a \code{data.frame} containing statistical results
#' on all combined synthetic results done with a fixed value of \code{D} (the
#' number of dimensions). The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{median}}{ a \code{numeric} representing the median of the
#' minimum AUROC obtained (within super populations) for all combination of
#' the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{mad}}{ a \code{numeric} representing the MAD of the minimum
#' AUROC obtained (within super populations) for all combination of the fixed
#' \code{D} value and all tested \code{K} values. }
#' \item{\code{upQuartile}}{ a \code{numeric} representing the upper quartile
#' of the minimum AUROC obtained (within super populations) for all
#' combination of the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{k}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for a fixed \code{D} value. }
#' }
#' }
#' \item{\code{dfPop}}{ a \code{data.frame} containing statistical results on
#' all combined synthetic results done with different values of \code{D} (the
#' number of dimensions) and \code{K} (the number of neighbors).
#' The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{AUROC.min}}{ a \code{numeric} representing the minimum accuracy
#' obtained by grouping all the synthetic results by super-populations, for
#' the specified values of \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing the accuracy obtained
#' by grouping all the synthetic results for the specified values of \code{D}
#' and \code{K}.}
#' \item{\code{Accu.CM}}{ a \code{numeric} representing the value of accuracy
#' of the confusion matrix obtained by grouping all the synthetic results for
#' the specified values of \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{dfAUROC}}{ a \code{data.frame} the summary of the results by
#' super-population. The \code{data.frame} contains
#' those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{Call}}{ a \code{character} string representing the
#' super-population.}
#' \item{\code{L}}{ a \code{numeric} representing the lower value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing  the AUROC obtained for
#' the fixed values of super-population, \code{D} and \code{K}.}
#' \item{\code{H}}{ a \code{numeric} representing the higher value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{D}}{ a \code{numeric} representing the optimal \code{D} value
#' (the number of dimensions) for the specific profile.}
#' \item{\code{K}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for the specific profile.}
#' \item{\code{listD}}{ a \code{numeric} representing the optimal \code{D}
#' values (the number of dimensions) for the specific profile. More than one
#' \code{D} is possible.}
#' }
#' }
#' \item{\code{KNNSample}}{  a \code{data.frame} containing the inferred
#' ancestry for different values of \code{K} and \code{D}. The
#' \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' }
#' }
#' \item{\code{KNNSynthetic}}{  a \code{data.frame} containing the inferred
#' ancestry for each synthetic data for different values of \code{K} and
#' \code{D}.
#' The \code{data.frame}
#' contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current synthetic data.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{infer.superPop}}{ a \code{character} string representing the
#' inferred ancestry for the specified \code{D} and \code{K} values.}
#' \item{\code{ref.superPop}}{ a \code{character} string representing the known
#' ancestry from the reference}
#' }
#' }
#' \item{\code{Ancestry}}{ a \code{data.frame} containing the inferred
#' ancestry for the current profile. The \code{data.frame} contains those
#' columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry.}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry.}
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry.}
#' }
#' }
#' \item{\code{pRAIDS}}{ a \code{list} of class \code{RAIDSparam} a list of
#' the parameters.}
#' }
#'
#' @details
#'
#' The runExomeAncestry() function generates 3 types of files
#' in the OUTPUT directory.
#' \describe{
#' \item{Ancestry Inference}{ The ancestry inference CSV file
#' (".Ancestry.csv" file)}
#' \item{Inference Informaton}{ The inference information RDS file
#' (".infoCall.rds" file)}
#' \item{Synthetic Information}{ The parameter information RDS files
#' from the synthetic inference ("KNN.synt.*.rds" files in a sub-directory)}
#' }
#'
#' In addition, a sub-directory (named using the profile ID) is
#' also created.
#'
#' @references
#'
#' Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
#' Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution
#' of ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
#' doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
#'
#' @examples
#'
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#'
#' #################################################################
#' ## The 1KG GDS file and the 1KG SNV Annotation GDS file
#' ## need to be located in the same directory
#' ## Note that the 1KG GDS file used for this example is a
#' ## simplified version and CANNOT be used for any real analysis
#' #################################################################
#' path1KG <- file.path(dataDir, "tests")
#'
#' fileReferenceGDS <- file.path(path1KG, "ex1_good_small_1KG.gds")
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#'
#' #################################################################
#' ## The Sample SNP pileup files (one per sample) need
#' ## to be located in the same directory.
#' #################################################################
#' demoProfileEx1 <- file.path(dataDir, "example", "snpPileup", "ex1.txt.gz")
#'
#' #################################################################
#' ## The path where the Profile GDS Files (one per sample)
#' ## will be created need to be specified.
#' #################################################################
#' pathProfileGDS <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathProfileGDS)){
#'     dir.create(pathProfileGDS)
#' }
#' ####################################################################
#' ## Fix seed to ensure reproducible results
#' ####################################################################
#' set.seed(3043)
#'
#' pRAIDS <- RAIDS:::RAIDSparam(profileFile=demoProfileEx1,
#'     genoSource="snp-pileup",
#'     pathProfileGDS=pathProfileGDS,
#'     fileReferenceGDS=fileReferenceGDS,
#'     fileReferenceAnnotGDS=fileAnnotGDS
#'     )
#'
#'
#'
#' res <- ancestryInferencePCAKNN(pRAIDS=pRAIDS)
#'
#' res$Ancestry
#'
#' unlink(file.path(pathProfileGDS, paste0(pedStudy(res$pRAIDS)$Name.ID[1], ".gds")) , force=TRUE)
#' if (length(dir(pathProfileGDS)) == 0) {
#' unlink(pathProfileGDS, recursive = TRUE)
#' }
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv
#' @importFrom gdsfmt index.gdsn read.gdsn openfn.gds ls.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
ancestryInferencePCAKNN <- function(pRAIDS=RAIDSparam()) {

    pRAIDS <- profileGenoSynthesis(pRAIDS=pRAIDS)

    ##########################
    ## Section ancestry
    ##########################


    spRef <- getRefSuperPop2(pRAIDS)

    # get a list of each pop.group in pRAIDS$syntheticRefDF
    sampleRM <- splitSelectByPop(syntheticRefDF(pRAIDS))



    if(verbose(pRAIDS)){
        message("SyntheticAncestry start ", Sys.time())
    }


    resSyn <- lapply(seq_len(nrow(sampleRM)),
                     FUN=function(x, sampleRM,
                                  pRAIDS, spRef) {
                         synthKNN <- computePoolSyntheticAncestryGr2(sampleRM=sampleRM[x,],
                                                                     spRef=spRef,
                                                                     pRAIDS = pRAIDS)
                         # If modify to keep the order in parallel
                         # synthKNN$matKNN$pos <- rep(x, nrow(synthKNN$matKNN))
                         # synthKNN$matKNN$spos <- seq_len(nrow(synthKNN$matKNN))
                         return(synthKNN$matKNN)
                     }, sampleRM=sampleRM, pRAIDS=pRAIDS,
                     spRef=spRef)

    if(verbose(pRAIDS)){
        message("SyntheticAncestry end ", Sys.time())
    }
    resSyn <- do.call(rbind, resSyn)
    # If modify to keep the order in parallel
    # resSyn <- resSyn[order(resSyn$pos, resSyn$spos),]
    # resSyn <- resSyn[,-1 * which(colnames(resSyn) %in% c("pos", "spos"))]
    pedSyn <- prepPedSyntheticRef(pRAIDS)

    resCall <- computeAncestryFromSynthetic2(syntheticKNN=resSyn,
                            pedSyn=pedSyn,
                            spRef=spRef,
                            listCatPop=unique(spRef),
                            pRAIDS=pRAIDS)

    if(verbose(pRAIDS)){
        message("Ancestry end ", Sys.time())
    }


    resSyn[[paste0("ref.superPop")]] <- pedSyn[resSyn$sample.id, "superPop"]

    colnames(resSyn) <- c("sample.id", "D", "K", "infer.superPop",
                          "ref.superPop")

    res <- list(pcaSample=resCall$pcaSample, # PCA of the profile + 1KG
                paraSample=resCall$paraSample, # Result of the parameter selection
                KNNSample=resCall$KNNSample$matKNN, # KNN for the profile
                KNNSynthetic=resSyn, # KNN results for synthetic data
                Ancestry=resCall$Ancestry,
                pRAIDS=pRAIDS) # the ancestry call fo the profile
    ## Successful
    return(res)
}

#' @title Create the gds file for the profile with the genotype and synthetic data
#'
#' @description First, the function creates the
#' Profile GDS file for the specific profile using the information from a
#' RDS Sample description file and the Population Reference GDS file.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{OL} when the function is successful.
#'
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#'
#' #################################################################
#' ## The 1KG GDS file and the 1KG SNV Annotation GDS file
#' ## need to be located in the same directory
#' ## Note that the 1KG GDS file used for this example is a
#' ## simplified version and CANNOT be used for any real analysis
#' #################################################################
#' path1KG <- file.path(dataDir, "tests")
#'
#' fileReferenceGDS <- file.path(path1KG, "ex1_good_small_1KG.gds")
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#'
#' #################################################################
#' ## The Sample SNP pileup files (one per sample) need
#' ## to be located in the same directory.
#' #################################################################
#' demoProfileEx1 <- file.path(dataDir, "example", "snpPileup", "ex1.txt.gz")
#'
#' #################################################################
#' ## The path where the Profile GDS Files (one per sample)
#' ## will be created need to be specified.
#' #################################################################
#' pathProfileGDS <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathProfileGDS)){
#'     dir.create(pathProfileGDS)
#' }
#'
#' ####################################################################
#' ## Fix seed to ensure reproducible results
#' ####################################################################
#' set.seed(3043)
#'
#' pRAIDS <- RAIDSparam(profileFile=demoProfileEx1,
#'     genoSource="snp-pileup",
#'     pathProfileGDS=pathProfileGDS,
#'     fileReferenceGDS=fileReferenceGDS,
#'     fileReferenceAnnotGDS=fileAnnotGDS
#'     )
#'
#'
#'
#' pRAIDS <- profileGenoSynthesis(pRAIDS=pRAIDS)
#'
#' unlink(file.path(pathProfileGDS, paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds")) , force=TRUE)
#' if (length(dir(pathProfileGDS)) == 0) {
#' unlink(pathProfileGDS, recursive = TRUE)
#' }
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv
#' @importFrom gdsfmt index.gdsn read.gdsn openfn.gds ls.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
profileGenoSynthesis <- function(pRAIDS) {
    # =RAIDSparam()
    profileBaseName <- basename(profileFile(pRAIDS))
    pathGeno <- dirname(profileFile(pRAIDS))

    if(is.null(studyType(pRAIDS)) || !(studyType(pRAIDS) %in% c("LD", "GeneAware"))){
        stop("The study must be specify as LD or GeneAware in the studyType or paramRAIDS.\n")
    }



    profileName <- gsub("\\.gz$", "", profileBaseName, ignore.case = TRUE)
    for(extCur in c( "\\.vcf$", "\\.txt$", "\\.bam", "\\.tsv", "\\.csv")){
        profileName <- gsub(extCur, "", profileName, ignore.case = TRUE)
    }
    tmp <- pedStudy(pRAIDS)

    tmp$Name.ID[1] <- profileName
    tmp$Case.ID[1] <- profileName
    rownames(tmp) <- tmp$Name.ID

    pedStudy(pRAIDS) <- tmp


    if(file.exists(file.path(pathProfileGDS(pRAIDS), paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds")))){
        stop(paste0("The gds file for ", pedStudy(pRAIDS)$Name.ID[1], " already exist."))
    }

    ## Generate the synthetic dataset if not included
    if(is.null(syntheticRefDF(pRAIDS))){
        syntheticRefDF(pRAIDS) <- selectForSynthetic(pRAIDS = pRAIDS,
                                                    nbProfiles=30L)
    }




    ##########################
    ## Create profile section
    ##########################

    ## parse the read.count if the pRAIDS$listPos is not
    ## it define it and return an update pRAIDS
    ## Create if note exists profileGDS
    ## add in the profile gds the "Ref.count", "Alt.count", and "Total.count"
    pRAIDS <- generateProfileRawGDS2(pRAIDS)

    ## Call the genotype from Ref.count and Alt.count
    generateProfileGenoCall(pRAIDS)

    ## Pruning
    pruningProfile2(pRAIDS)
    ## add genotype structure for the pruning positions
    addGenotypeProfile(pRAIDS)




    ##########################
    ## Section allelic fraction
    ##########################

    studyTypeLeg <- ifelse(studyType(pRAIDS) == "LD", "DNA", "RNA")
    ## Profile GDS file name
    fileProfileGDS <-  validateProfileGDSExist(pathProfile=pathProfileGDS(pRAIDS),
                                               profile=pedStudy(pRAIDS)$Name.ID[1])
    ## Open Profile GDS file
    gdsProfile <- snpgdsOpen(fileProfileGDS, readonly=FALSE)

    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))
    gdsRefAnnot <- openfn.gds(fileReferenceAnnotGDS(pRAIDS))

    estimateAllelicFraction(gdsReference=gdsReference, gdsProfile=gdsProfile,
                            currentProfile=pedStudy(pRAIDS)$Name.ID[1], studyID=studyDF(pRAIDS)$study.id,
                            chrInfo=chrInfo(pRAIDS), studyType=studyTypeLeg, gdsRefAnnot=gdsRefAnnot,
                            blockID=blockTypeId(pRAIDS), verbose=verbose(pRAIDS))
    snpgdsClose(gdsProfile)
    snpgdsClose(gdsReference)
    closefn.gds(gdsRefAnnot)

    ##########################
    ## Section synthetic
    ##########################

    ## Add information related to the synthetic profiles in Profile GDS file
    prepSynthetic(fileProfileGDS=fileProfileGDS,
                  listSampleRef=syntheticRefDF(pRAIDS)$sample.id,  profileID=pedStudy(pRAIDS)$Name.ID[1],
                  studyDF=studyDFSyn(pRAIDS), prefix=prefix(pRAIDS), verbose=verbose(pRAIDS))

    if(verbose(pRAIDS)){
        message("syntheticGeno start ", Sys.time())
    }


    resG <-syntheticGeno2(pRAIDS)


    ## Successful
    return(pRAIDS)
}


#' @title Run a k-nearest neighbors analysis on a subset of the
#' synthetic dataset
#'
#' @description The function runs k-nearest neighbors analysis on a
#' subset of the synthetic data set. The function uses the 'knn' package.
#'
#' @param listEigenvector a \code{list} with 3 entries:
#' 'sample.id', 'eigenvector.ref' and 'eigenvector'. The \code{list} represents
#' the PCA done on the 1KG reference profiles and the synthetic profiles
#' projected onto it.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{vector} of \code{character} strings
#' representing the identifiers of the synthetic profiles analysed.}
#' \item{\code{sample1Kg}}{ a \code{vector} of \code{character} strings
#' representing the identifiers of the 1KG reference profiles used to
#' generate the synthetic profiles.}
#' \item{\code{sp}}{ a \code{vector} of \code{character} strings representing
#' the known super population ancestry of the 1KG reference profiles used
#' to generate the synthetic profiles.}
#' \item{\code{matKNN}}{ a \code{data.frame} containing the super population
#' inference for each synthetic profiles for different values of PCA
#' dimensions \code{D} and k-neighbors values \code{K}. The fourth column title
#' corresponds to the \code{fieldPopInfAnc} parameter.
#' The \code{data.frame} contains 4 columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing
#' the identifier of the synthetic profile analysed.}
#' \item{\code{D}}{ a \code{numeric} strings representing
#' the value of the PCA dimension used to infer the super population.}
#' \item{\code{K}}{ a \code{numeric} strings representing
#' the value of the k-neighbors used to infer the super population.}
#' \item{\code{fieldPopInfAnc} value}{ a \code{character} string representing
#' the inferred ancestry.}
#' }
#' }
#' }
#'
#' @examples
#'
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS)
#' ## Projects synthetic profiles on 1KG PCA
#' results <- computeKNNRefSyn(listEigenvector=demoPCASyntheticProfiles,
#'     pRAIDS=pRAIDS)
#'
#' ## The inferred ancestry for the synthetic profiles for different values
#' ## of D and K
#' head(results$matKNN)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn openfn.gds closefn.gds
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom class knn
#' @encoding UTF-8
#' @export
computeKNNRefSyn <- function(listEigenvector, pRAIDS) {

    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                                    paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(filename=fileGDSProfile)


    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))
    dfRef <- read.gdsn(index.gdsn(gdsReference, "sample.annot"))
    row.names(dfRef) <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
    # If pRAIDS$sampleRef is define we restricted the sample to pRAIDS$sampleRef
    if(! is.null(sampleRef(pRAIDS))){
        dfRef <- dfRef[sampleRef(pRAIDS),]
    }

    # listCatPop <- unique(spRef)
    snpgdsClose(gdsReference)

    ## Get study information from the GDS profile file
    studyAnnotAll <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    snpgdsClose(gdsProfile)

    studyAnnot <- studyAnnotAll[which(studyAnnotAll$study.id ==
                    studyDFSyn(pRAIDS)$study.id & studyAnnotAll$data.id %in%
                    listEigenvector$sample.id), ]
    pcaList <- pcaList(pRAIDS)[pcaList(pRAIDS) <= eigenCountSyn(pRAIDS)]

    listMat <- lapply(pcaList,
                    FUN=function(pcaD, listEigenvector, dfRef, pRAIDS){
                        resMat <- lapply(kList(pRAIDS),
                            FUN=function(kV, pcaD, listEigenvector,
                                    dfRef, pRAIDS){
                                return(computeKNNProfileSubSet(listEigenvector,
                                                kV, pcaD,dfRef, pRAIDS))
                            },
                            pcaD=pcaD,
                            listEigenvector=listEigenvector,
                            dfRef=dfRef,
                            pRAIDS=pRAIDS)
                        resMat <- do.call(rbind, resMat)
                        return(resMat)
                      },
                      listEigenvector=listEigenvector,
                      dfRef=dfRef,
                      pRAIDS=pRAIDS)

    resMat <- do.call(rbind, listMat)

    listKNN <- list(sample.id=listEigenvector$sample.id,
                    sample1Kg=studyAnnot$case.id,
                    sp=dfRef[studyAnnot$case.id,],
                    matKNN=resMat)

    return(listKNN)
}


#' @title Run a PCA analysis and a K-nearest neighbors analysis on a small set
#' of synthetic data using all 1KG profiles except the ones used to generate
#' the synthetic profiles
#'
#' @description The function runs a PCA analysis using 1 synthetic profile
#' from each sub-continental population. The reference profiles used to
#' create those synthetic profiles are first removed from the list
#' of 1KG reference profiles that generates the reference PCA. Then, the
#' retained synthetic
#' profiles are projected on the 1KG PCA space. Finally, a K-nearest neighbors
#' analysis using a range of K and D values is done.
#'
#' @param sampleRM a \code{vector} of \code{character} strings representing
#' the identifiers of the 1KG reference profiles that should not be used to
#' create the reference PCA. There should be one per sub-continental
#' population. Those profiles are
#' removed because those have been used to generate the synthetic profiles
#' that are going to be analysed here. The sub-continental
#' identifiers are used as names for the \code{vector}.
#'
#' @param spRef \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param pRAIDS a \code{RAIDSparam} an object containing all RAIDS
#' parameters.
#'
#' @return a \code{list} containing the following entries:
#' \describe{
#' \item{sample.id}{ a \code{vector} of \code{character} strings representing
#' the identifiers of the synthetic profiles. }
#' \item{sample1Kg}{ a \code{vector} of \code{character} strings representing
#' the identifiers of the reference 1KG profiles used to generate the
#' synthetic profiles. }
#' \item{sp}{ a \code{vector} of \code{character} strings representing the
#' known ancestry for the reference 1KG profiles used to generate the
#' synthetic profiles. }
#' \item{matKNN}{ a \code{data.frame} containing 4 columns. The first column
#' 'sample.id' contains the name of the synthetic profile. The second column
#' 'D' represents the dimension D used to infer the ancestry. The third column
#' 'K' represents the number of neighbors K used to infer the ancestry. The
#' fourth column 'SuperPop' contains the inferred ancestry. }
#' }
#'
#' @references
#'
#' Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
#' Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution
#' of ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
#' doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
#'
#' @examples
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' samplesRM <- c("HG00246", "HG00325", "HG00611", "HG01173", "HG02165",
#'     "HG01112", "HG01615", "HG01968", "HG02658", "HG01850", "HG02013",
#'     "HG02465", "HG02974", "HG03814", "HG03445", "HG03689", "HG03789",
#'     "NA12751", "NA19107", "NA18548", "NA19075", "NA19475", "NA19712",
#'     "NA19731", "NA20528", "NA20908")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS)
#'
#' result <- computePoolSyntheticAncestryGr2(sampleRM=samplesRM,
#'                              spRef=demoKnownSuperPop1KG,
#'                              pRAIDS=pRAIDS)
#' ## The ancestry inference for the synthetic data using
#' ## different K and D values
#' head(result$matKNN)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom rlang arg_match
#' @importFrom gdsfmt openfn.gds closefn.gds
#' @encoding UTF-8
#' @export
computePoolSyntheticAncestryGr2 <- function(sampleRM, spRef,
    pRAIDS) {

    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                                paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(filename=fileGDSProfile)
    ## Validate the input parameters
    ## TODO just validate sampleRM and spRef
    validateComputePoolSyntheticAncestryGr(gdsProfile=gdsProfile,
        sampleRM=sampleRM, spRef=spRef,
        studyIDSyn=studyDFSyn(pRAIDS)$study.id,
        np=np(pRAIDS),
        listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
        pcaList=pcaList(pRAIDS),
        fieldPopInfAnc=fieldPopInfAnc(pRAIDS),
        kList=kList(pRAIDS),
        algorithm=PCAalgorithm(pRAIDS),
        eigenCount=eigenCountSyn(pRAIDS),
        missingRate=PCAmissingRate(pRAIDS),
        verbose=verbose(pRAIDS))

    if(verbose(pRAIDS)){
        message("computePoolSyntheticAncestryGr2 p1 ", Sys.time())
    }
    # Get the list of syn profiles to process
    study.annot <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    study.annot <- study.annot[which(study.annot$study.id == studyDFSyn(pRAIDS)$study.id &
                        study.annot$case.id %in% sampleRM),]
    closefn.gds(gdsProfile)
    ## Set algorithm
    #algorithm <- arg_match(algorithm)

    ## Calculate Principal Component Analysis (PCA) on SNV genotype dataset
    ## excluded the selected profiles used to generate the synthetic profiles
    if(verbose(pRAIDS)){
        message("computePoolSyntheticAncestryGr2 p2 ", Sys.time())
    }
    pcaRef <- computePCARefRMMulti2(listRM=sampleRM, pRAIDS=pRAIDS)

    ## Calculate PCA on the synthetic profiles using 1KG PCA results

    # resPCA <- computePCAMultiSynthetic(gdsProfile=gdsProfile, listPCA=pca1KG,
    #                                    sampleRef=sampleRM, studyIDSyn=studyIDSyn, verbose=verbose)
    if(verbose(pRAIDS)){
        message("computePoolSyntheticAncestryGr2 p3 ", Sys.time())
    }
    resPCA <- computePCAProfile(listPCA=pcaRef, profileId =study.annot$data.id , pRAIDS)
    ## Calculate the k-nearest neighbor analyses on a subset of the
    ## synthetic data set
    # synthKNN <- computeKNNRefSynthetic(gdsProfile=gdsProfile,
    #                                    listEigenvector=resPCA,
    #                                    listCatPop=listCatPop, studyIDSyn=studyIDSyn, spRef=spRef,
    #                                    fieldPopInfAnc=fieldPopInfAnc, kList=kList, pcaList=pcaList)
    if(verbose(pRAIDS)){
        message("computePoolSyntheticAncestryGr2 p4 ", Sys.time())
    }
    synthKNN <- computeKNNRefSyn(listEigenvector=resPCA,
                                 pRAIDS=pRAIDS)
    if(verbose(pRAIDS)){
        message("computePoolSyntheticAncestryGr2 end ", Sys.time())
    }
    return(synthKNN)
}

#' @title Select the optimal K and D parameters using the synthetic data and
#' infer the ancestry of a specific profile
#'
#' @description The function select the optimal K and D parameters for a
#' specific profile. The results on the synthetic data are used for the
#' parameter selection. Once the optimal parameters are selected, the
#' ancestry is inferred for the specific profile.
#'
#' @param syntheticKNN a \code{vector} of \code{character} strings representing
#' the name of files that contain the results of ancestry inference done on
#' the synthetic profiles for multiple values of _D_ and _K_. The files must
#' exist.
#'
#' @param pedSyn a \code{data.frame} containing the columns extracted from the
#' GDS Sample 'study.annot' node with a extra column named as the 'popName'
#' parameter that has been extracted from the 1KG GDS 'sample.annot' node.
#'
#' @param spRef a \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{pcaSample}}{ a \code{list} containing the information related
#' to the eigenvectors. The \code{list} contains those 3 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{eigenvector.ref}}{ a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{\code{eigenvector}}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current profile projected on the PCA from the
#' reference profiles.}
#' }
#' }
#' \item{\code{paraSample}}{ a \code{list} containing the results with
#' different \code{D} and \code{K} values that lead to optimal parameter
#' selection. The \code{list} contains those entries:
#' \describe{
#' \item{\code{dfPCA}}{ a \code{data.frame} containing statistical results
#' on all combined synthetic results done with a fixed value of \code{D} (the
#' number of dimensions). The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{median}}{ a \code{numeric} representing the median of the
#' minimum AUROC obtained (within super populations) for all combination of
#' the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{mad}}{ a \code{numeric} representing the MAD of the minimum
#' AUROC obtained (within super populations) for all combination of the fixed
#' \code{D} value and all tested \code{K} values. }
#' \item{\code{upQuartile}}{ a \code{numeric} representing the upper quartile
#' of the minimum AUROC obtained (within super populations) for all
#' combination of the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{k}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for a fixed \code{D} value. }
#' }
#' }
#' \item{\code{dfPop}}{ a \code{data.frame} containing statistical results on
#' all combined synthetic results done with different values of \code{D} (the
#' number of dimensions) and \code{K} (the number of neighbors).
#' The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{AUROC.min}}{ a \code{numeric} representing the minimum accuracy
#' obtained by grouping all the synthetic results by super-populations, for
#' the specified values of \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing the accuracy obtained
#' by grouping all the synthetic results for the specified values of \code{D}
#' and \code{K}.}
#' \item{\code{Accu.CM}}{ a \code{numeric} representing the value of accuracy
#' of the confusion matrix obtained by grouping all the synthetic results for
#' the specified values of \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{dfAUROC}}{ a \code{data.frame} the summary of the results by
#' super-population. The \code{data.frame} contains
#' those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{Call}}{ a \code{character} string representing the
#' super-population.}
#' \item{\code{L}}{ a \code{numeric} representing the lower value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing  the AUROC obtained for the
#' fixed values of super-population, \code{D} and \code{K}.}
#' \item{\code{H}}{ a \code{numeric} representing the higher value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{D}}{ a \code{numeric} representing the optimal \code{D} value
#' (the number of dimensions) for the specific profile.}
#' \item{\code{K}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for the specific profile.}
#' \item{\code{listD}}{ a \code{numeric} representing the optimal \code{D}
#' values (the number of dimensions) for the specific profile. More than one
#' \code{D} is possible.}
#' }
#' }
#' \item{\code{KNNSample}}{ a \code{list} containing the inferred ancestry
#' using different \code{D} and \code{K} values. The \code{list} contains
#' those entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{matKNN}}{ a \code{data.frame} containing the inferred ancestry
#' for different values of \code{K} and \code{D}. The \code{data.frame}
#' contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' }
#' }
#' }
#' }
#' \item{\code{Ancestry}}{ a \code{data.frame} containing the inferred
#' ancestry for the current profile. The \code{data.frame} contains those
#' columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry.}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry.}
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry.}
#' }
#' }
#' }
#'
#'
#' @references
#'
#' Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
#' Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution
#' of ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
#' doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
#'
#' @examples
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## The Reference GDS file
#' path1KG <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#' pedSyn <- RAIDS:::prepPedSyntheticRef(pRAIDS)
#' spRef <- getRefSuperPop2(pRAIDS)
#' syntheticKNN <- readRDS(file.path(dataDir, "ex1", "KNN.synt.ex1.2.rds"))
#' resCall <- computeAncestryFromSynthetic2(syntheticKNN=syntheticKNN,
#'                       pedSyn=pedSyn,
#'                       spRef=spRef,
#'                       listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
#'                       pRAIDS=pRAIDS)
#'
#' resCall$Ancestry
#'
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
computeAncestryFromSynthetic2 <- function(syntheticKNN,
                                         pedSyn,
                                         spRef,
                                         listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                                         pRAIDS=pRAIDS) {

    ## Validate input parameters
    # validateComputeAncestryFromSynthetic(gdsReference=gdsReference,
    #                                          gdsProfile=gdsProfile, syntheticKNN=syntheticKNN,
    #                                          pedSyn=pedSyn,
    #                                          currentProfile=currentProfile, spRef=spRef, studyIDSyn=studyIDSyn,
    #                                          np=np, listCatPop=listCatPop, fieldPopInRef=fieldPopInRef,
    #                                          fieldPopInfAnc=fieldPopInfAnc, kList=kList, pcaList=pcaList,
    #                                          algorithm=algorithm, eigenCount=eigenCount, missingRate=missingRate,
    #                                          verbose=verbose)

    ## Matches a character method against a table of candidate values
    ## aalready validate
    #algorithm <- arg_match(algorithm)



    ## Compile all the inferred ancestry results for different values of
    ## D and K to select the optimal parameters
    listParaSample <- selParaPCAUpQuartile(matKNN=syntheticKNN,
                                           pedCall=pedSyn, refCall=fieldPopInRef(pRAIDS), predCall=fieldPopInfAnc(pRAIDS),
                                           listCall=listCatPop,kList=kList(pRAIDS),
                                           pcaList=pcaList(pRAIDS)[pcaList(pRAIDS) <= eigenCountSyn(pRAIDS)])
    listPCARef <- computePCARefRMMulti2(listRM = NULL, pRAIDS = pRAIDS)
    listPCAProfile <- computePCAProfile(listPCA=listPCARef, profileId=pedStudy(pRAIDS)$Name.ID[1], pRAIDS=pRAIDS)
    ## Project profile on the PCA created with the reference profiles
    # listPCAProfile <- computePCARefSample(gdsProfile=gdsProfile,
    #                                       currentProfile=currentProfile, studyIDRef="Ref.1KG", np=np,
    #                                       algorithm=algorithm, eigenCount=eigenCount, missingRate=missingRate,
    #                                       verbose=verbose)

    ## Run a k-nearest neighbors analysis on one specific profile
    listKNNSample <- computeKNNRefSample(listEigenvector=listPCAProfile,
                                         listCatPop=listCatPop, spRef=spRef,
                                         fieldPopInfAnc=fieldPopInfAnc(pRAIDS),
                                         kList=kList(pRAIDS), pcaList=pcaList(pRAIDS))

    ## The ancestry call for the current profile
    resCall <- listKNNSample$matKNN[
        which(listKNNSample$matKNN$D == listParaSample$D &
                  listKNNSample$matKNN$K == listParaSample$K ),]
    colnames(listParaSample$dfAUROC) <- c("D", "K", "Call", "L", "AUROC", "H")
    res <- list(pcaSample=listPCAProfile, # PCA of the profile + 1KG
                paraSample=listParaSample, # Result of the parameter selection
                KNNSample=listKNNSample, # KNN for the profile
                Ancestry=resCall) # the ancestry call fo the profile

    return(res)
}

#' @title Extract the sample information from the 1KG GDS file for a list
#' of profiles associated to a specific study in the Profile GDS file
#'
#' @description The function extracts the information for the profiles
#' associated to a specific study in the GDS Sample file. The information is
#' extracted from the 'study.annot' node as a 'data.frame'.
#'
#' Then, the function used the 1KG GDS file to extract specific information
#' about each sample and add it, as an extra column, to the 'data.frame'.
#'
#' As example, this function can extract the synthetic profiles
#' for a GDS Sample and the super-population of the 1KG samples used to
#' generate each synthetic profile would be added
#' as an extra column to the final 'data.frame'.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{data.frame} containing the columns extracted from the
#' GDS Sample 'study.annot' node with a extra column named as the 'popName'
#' parameter that has been extracted from the 1KG GDS 'sample.annot' node.
#' Only the rows corresponding to the specified study ('studyID' parameter)
#' are returned.
#'
#'
#' @details
#'
#' As example, this function can extract the synthetic profiles
#' for a Profile GDS and the super-population of the 1KG samples used to
#' generate each synthetic profile would be added
#' as an extra column to the final 'data.frame'. In that situation, the
#' 'popName' parameter would correspond to the super-population column and the
#' 'studyID' parameter would be the name given to the synthetic dataset.
#'
#'
#' @examples
#'
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## The Reference GDS file
#' path1KG <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#' pedSyn <- RAIDS:::prepPedSyntheticRef(pRAIDS)
#' head(pedSyn)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @export
prepPedSyntheticRef <- function(pRAIDS) {

    studyID <- studyDFSyn(pRAIDS)$study.id
    popName <- fieldPopInRef(pRAIDS)
    fileProfileGDS <-  validateProfileGDSExist(pathProfile=pathProfileGDS(pRAIDS),
                                    profile=pedStudy(pRAIDS)$Name.ID[1])
    ## Open Profile GDS file
    gdsProfile <- snpgdsOpen(fileProfileGDS, readonly=TRUE)

    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))

    ## Extract study information from the Profile GDS file
    studyAnnot <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))

    ## Retain the information associated to the current study
    studyCur <- studyAnnot[which(studyAnnot$study.id == studyID),]
    rm(studyAnnot)
    if(studyCur$case.id[1] == studyCur$data.id[1]){
        tmp <- matrix(unlist(strsplit(x = studyCur$case.id, split = "\\.")),
                      nrow=4)
        studyCur$case.id <- tmp[3,]
    }
    ## Get the information from 1KG GDS file
    dataRef <- read.gdsn(index.gdsn(node=gdsReference, "sample.annot"))

    if(! popName %in% colnames(dataRef)) {
        stop("The population ", popName, " is not supported.")
    }

    ## Assign sample names to the information
    row.names(dataRef) <- read.gdsn(index.gdsn(node=gdsReference, "sample.id"))

    studyCur[[popName]] <- dataRef[studyCur$case.id, popName]
    rownames(studyCur) <- studyCur$data.id

    ## Important to close the files
    closefn.gds(gdsProfile)
    closefn.gds(gdsReference)

    return(studyCur)
}


#' @title compute admixture proportion
#'
#' @description TODO
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 2 entries:
#' \describe{
#' \item{\code{profileQ}}{ a \code{matrix} containing the admixture
#' proportions of the profile }
#' \item{\code{syntheticQ}}{ a \code{matrix} containing for each
#' synthetic profile the admixture proportions of the synthetic profile,
#' the synthetic proportions of the reference using the same position as
#' the profile, and the admixture proportion of the gold standard}
#' }
#'
#'
#' @details TODO
#'
#'
#'
#'
#' @examples
#'
#' \dontrun{
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' pathOut <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathOut)){
#'     dir.create(pathOut)
#' }
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(pathOut) # , "ex1.gds"
#' file.copy(file.path(dataDir, "ex1.gds"), file.path(pathOut))
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#' ## Path to admixture
#' ## pathADMIXTURE <- "PATH/admixture"
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop",
#'                       pathADMIXTURE=pathADMIXTURE)
#' res <- computeAdmixtureProportions(pRAIDS=pRAIDS)
#'
#'
#'
#' dir(file.path(pathProfileGDS))
#' head(res$profileQ)
#' head(res$syntheticQ)
#'
#' unlink(file.path(pathOut), recursive = TRUE)#'
#' }
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @export
computeAdmixtureProportions <- function(pRAIDS) {
    if(is.null(pathADMIXTURE(pRAIDS)) || ! file.exists(pathADMIXTURE(pRAIDS))){
        stop("Fatal, you need to defined pathADMIXTURE in RAIDSparam to use ADMIXTURE")
    }

    fileProfileGDS <-  validateProfileGDSExist(pathProfile=pathProfileGDS(pRAIDS),
                                    profile=pedStudy(pRAIDS)$Name.ID[1])

    pathOut <- file.path(pathProfileGDS(pRAIDS), paste0("tmp.", pedStudy(pRAIDS)$Name.ID[1]))
    if(! dir.exists(pathOut)){
        dir.create(pathOut)
    }
    matPIndex <- computePmatrix(pathOut=pathOut, pRAIDS=pRAIDS)

    generateProfilePlinkFiles(pathOut=pathOut, pRAIDS=pRAIDS)
    profileQ <- runADMIXTURE(pathOut = file.path(pathOut, "profileAdmixture"),
            fileBed = paste0(pedStudy(pRAIDS)$Name.ID[1], ".bed"),
            seedAdmix=seedADMIXTURE(pRAIDS), pRAIDS)

    syntheticQ <- runSynADMIXTURE(pathOut=pathOut, matPIndex=matPIndex, pRAIDS=pRAIDS)
    syntheticQ <- syntheticJSD(syntheticQ=syntheticQ, pRAIDS=pRAIDS)
    if(cleanADMIXTUREFiles(pRAIDS)){
        unlink(pathOut, recursive = TRUE)
    }
    return(list(profileQ=profileQ, syntheticQ=syntheticQ))
}



#' @title Create the GDS Sample file(s) for one or multiple specific samples
#' using the information from a RDS Sample description file and the 1KG
#' GDS file
#'
#' @description The function uses the information for the 1KG GDS file and the
#' RDS Sample Description file to create the GDS Sample file. One GDS Sample
#' file is created per sample. One GDS Sample file will be created for each
#' entry present in the \code{listSamples} parameter.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the VCF output of SNP-pileup for each sample. The
#' SNP-pileup files must be compressed (gz files) and have the name identifiers
#' of the samples. A sample with "Name.ID" identifier would have an
#' associated SNP-pileup file called "Name.ID.txt.gz".
#'
#' @param fileNamePED a \code{character} string representing the path to the
#' RDS file that contains the information about the sample to analyse.
#' The RDS file must
#' include a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings. The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{fileNamePED} or \code{pedStudy}
#' can be defined.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{fileNamePED} or \code{pedStudy}
#' can be defined.
#'
#' @param fileNameGDS a \code{character} string representing the file name of
#' the 1KG GDS file. The file must exist.
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
#' Default: \code{1}.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param listProfiles a \code{vector} of \code{character} string corresponding
#' to the profile identifiers that will have a Profile GDS file created. The
#' profile identifiers must be present in the "Name.ID" column of the Profile
#' RDS file passed to the \code{fileNamePED} parameter.
#' If \code{NULL}, all profiles present in the \code{fileNamePED} are selected.
#' Default: \code{NULL}.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the Profile GDS files will be created.
#' Default: \code{NULL}.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{TRUE}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(data.dir, "ex1_good_small_1KG_GDS.gds")
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
#' ## Create the Sample GDS file for sample in listSamples vector
#' ## (in this case, samples "ex1")
#' ## The Profile GDS file is created in the pathProfileGDS directory
#' result <- createStudy2GDS1KG(pathGeno=data.dir,
#'             pedStudy=samplePED, fileNameGDS=fileGDS,
#'             studyDF=studyDF, listProfiles=c("ex1"),
#'             pathProfileGDS=data.dir, verbose=FALSE)
#'
#' ## The function returns OL when successful
#' result
#'
#' ## The Profile GDS file 'ex1.gds' has been created in the
#' ## specified directory
#' list.files(data.dir)
#'
#' ## Unlink Profile GDS file (created for demo purpose)
#' unlink(file.path(data.dir, "ex1.gds"))
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
createStudy2GDS1KG <- function(pathGeno=file.path("data", "sampleGeno"),
                                fileNamePED=NULL, pedStudy=NULL, fileNameGDS,
                                batch=1, studyDF, listProfiles=NULL,
                                pathProfileGDS=NULL, verbose=TRUE) {

    ## When fileNamePED is defined and pedStudy is null
    if (!(is.null(fileNamePED)) && is.null(pedStudy)) {
        ## The fileNamePED must be a character string and the file must exists
        if (!(is.character(fileNamePED) && (file.exists(fileNamePED)))) {
            stop("The \'fileNamePED\' must be a character string representing",
                    " the RDS Sample information file. The file must exist.")
        }
        ## Open the RDS Sample information file
        pedStudy <- readRDS(file=fileNamePED)
    } else if (!(is.null(fileNamePED) || is.null(pedStudy))) {
        stop("Both \'fileNamePED\' and \'pedStudy\' parameters cannot be ",
                "defined at the same time.")
    } else if (is.null(fileNamePED) && is.null(pedStudy)) {
        stop("One of the parameter \'fineNamePED\' of \'pedStudy\' must ",
                "be defined.")
    }

    ## Validate input parameters
    validateCreateStudy2GDS1KG(pedStudy=pedStudy, fileNameGDS=fileNameGDS,
                batch=batch, studyDF=studyDF, listProfiles=listProfiles,
                pathProfileGDS=pathProfileGDS, verbose=verbose)

    ## Read the 1KG GDS file
    gds <- snpgdsOpen(filename=fileNameGDS)

    ## Extract the chromosome and position information for all SNPs in 1KG GDS
    snpCHR <- index.gdsn(node=gds, "snp.chromosome")
    snpPOS <- index.gdsn(node=gds, "snp.position")

    ## Create a data.frame containing the information
    listPos <- data.frame(snp.chromosome=read.gdsn(snpCHR),
                            snp.position=read.gdsn(snpPOS))

    if(verbose) {
        message("Start ", Sys.time())
        message("Sample info DONE ", Sys.time())
    }

    generateGDS1KGgenotypeFromSNPPileup(pathGeno=pathGeno,
        listSamples=listProfiles, listPos=listPos, offset=-1, minCov=10,
        minProb=0.999, seqError=0.001, pedStudy=pedStudy, batch=batch,
        studyDF=studyDF, PATHGDSSAMPLE=pathProfileGDS, verbose=verbose)

    if(verbose) {
        message("Genotype DONE ", Sys.time())
    }

    ## Close 1KG GDS file
    closefn.gds(gds)

    ## Return successful code
    return(0L)
}

#' @title Create the GDS Sample file(s) for one or multiple specific samples
#' using the information from a RDS Sample description file and the 1KG
#' GDS file
#'
#' @description The function uses the information for the 1KG GDS file and the
#' RDS Sample Description file to create the GDS Sample file. One GDS Sample
#' file is created per sample. One GDS Sample file will be created for each
#' entry present in the \code{listSamples} parameter.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the output of SNP-pileup, a VCF Sample file, for
#' each sample. The
#' SNP-pileup files must be compressed (gz files) and have the name identifiers
#' of the samples. A sample with "Name.ID" identifier would have an
#' associated SNP-pileup file called "Name.ID.txt.gz".
#'
#' @param fileNamePED a \code{character} string representing the path to the
#' RDS file that contains the information about the sample to analyse.
#' The RDS file must
#' include a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings. The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter.
#'
#' @param fileNameGDS a \code{character} string representing the file name of
#' the 1KG GDS file. The file must exist.
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
#' Default: \code{1}.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings.
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to the sample identifiers that will have a GDS Sample file created. The
#' sample identifiers must be present in the "Name.ID" column of the RDS file
#' passed to the \code{fileNamePED} parameter.
#' If \code{NULL}, all samples in the \code{fileNamePED} are selected.
#' Default: \code{NULL}.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Sample files will be created.
#' Default: \code{NULL}.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{TRUE}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyInfo <- data.frame(study.id="Pancreatic.WES",
#'                 study.desc="Pancreatic study",
#'                 study.platform="WES",
#'                 stringsAsFactors=FALSE)
#'
#' ## TODO
#' fileNamePED <- "TODO"
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds read.gdsn
#' @encoding UTF-8
#' @export
appendStudy2GDS1KG <- function(pathGeno=file.path("data", "sampleGeno"),
                                fileNamePED, fileNameGDS, batch=1,
                                studyDF, listSamples=NULL,
                                pathProfileGDS=NULL, verbose=TRUE) {

    ## Validate inputs
    validateAppendStudy2GDS1KG(pathGeno=pathGeno, fileNamePED=fileNamePED,
        fileNameGDS=fileNameGDS, batch=batch, studyDF=studyDF,
        listSamples=listSamples, pathProfileGDS=pathProfileGDS, verbose=verbose)

    ## Open the RDS Sample information file
    pedStudy <- readRDS(file=fileNamePED)

    ## Read the 1KG GDS file
    gds <- snpgdsOpen(filename=fileNameGDS)

    ## Extract the chromosome and position information for all SNPs in 1KG GDS
    ## Create a data.frame containing the information
    snpCHR <- index.gdsn(node=gds, "snp.chromosome")
    snpPOS <- index.gdsn(node=gds, "snp.position")

    listPos <- data.frame(snp.chromosome=read.gdsn(snpCHR),
                            snp.position=read.gdsn(snpPOS))

    if (verbose) {
        message("Start ", Sys.time())
        message("Sample info DONE ", Sys.time())
    }

    generateGDS1KGgenotypeFromSNPPileup(pathGeno=pathGeno,
        listSamples=listSamples, listPos=listPos, offset=-1,
        minCov=10, minProb=0.999, seqError=0.001, pedStudy=pedStudy,
        batch=batch, studyDF=studyDF, PATHGDSSAMPLE=pathProfileGDS)

    if (verbose) {
        message("Genotype DONE ", Sys.time())
    }

    ## Close 1KG GDS file
    closefn.gds(gds)

    ## Return successful code
    return(0L)
}


#' @title Compute the list of pruned SNVs for a specific profile using the
#' information from the 1KG GDS file and a linkage disequilibrium analysis
#'
#' @description This function computes the list of pruned SNVs for a
#' specific profile. When
#' a group of SNVs are in linkage disequilibrium, only one SNV from that group
#' is retained. The linkage disequilibrium is calculated with the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The initial list of
#' SNVs that are passed to the \code{\link[SNPRelate]{snpgdsLDpruning}}()
#' function can be specified by the user.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1 KG GDS file (reference data set).
#'
#' @param method a \code{character} string that represents the method that will
#' be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime" and "composite". Default: \code{"corr"}.
#'
#' @param currentProfile  a \code{character} string
#' corresponding to the profile identifier used in LD pruning done by the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. A Profile GDS file
#' corresponding to the profile identifier must exist and be located in the
#' \code{pathProfileGDS} directory.
#'
#' @param studyID a \code{character} string corresponding to the study
#' identifier used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#' The study identifier must be present in the Profile GDS file.
#'
#' @param listSNP a \code{vector} of SNVs identifiers specifying selected to
#' be passed the the pruning function;
#' if \code{NULL}, all SNVs are used in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}} function. Default: \code{NULL}.
#'
#' @param slideWindowMaxBP a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done in the \code{\link[SNPRelate]{snpgdsLDpruning}}
#' function.
#' Default: \code{500000L}.
#'
#' @param thresholdLD a single \code{numeric} value that represents the LD
#' threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#' Default: \code{sqrt(0.1)}.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param verbose a \code{logicial} indicating if information is shown
#' during the process in the \code{\link[SNPRelate]{snpgdsLDpruning}}
#' function.  Default: \code{FALSE}.
#'
#' @param chr a \code{character} string representing the chromosome where the
#' selected SNVs should belong. Only one chromosome can be handled. If
#' \code{NULL}, the chromosome is not used as a filtering criterion.
#' Default: \code{NULL}.
#'
#' @param superPopMinAF a single positive \code{numeric} representing the
#' minimum allelic frequency used to select the SNVs. If \code{NULL}, the
#' allelic frequency is not used as a filtering criterion. Default: \code{NULL}.
#'
#' @param keepPrunedGDS a \code{logicial} indicating if the information about
#' the pruned SNVs should be added to the GDS Sample file.
#' Default: \code{TRUE}.
#'
#' @param pathProfileGDS a \code{character} string representing the directory
#' where the Profile GDS files will be created. The directory must exist.
#'
#' @param keepFile a \code{logical} indicating if RDS files containing the
#' information about the pruned SNVs must be
#' created. Default: \code{FALSE}.
#'
#' @param pathPrunedGDS a \code{character} string representing an existing
#' directory. The directory must exist. Default: \code{"."}.
#'
#' @param outPrefix a \code{character} string that represents the prefix of the
#' RDS files that will be generated. The RDS files are only generated when
#' the parameter \code{keepFile}=\code{TRUE}. Default: \code{"pruned"}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(data.dir, "ex1_good_small_1KG_GDS.gds")
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
#' ## Copy the Profile GDS file demo that has not been pruned yet
#' ## into a test directory (deleted after the example has been run)
#' data.dir.pruning <- file.path(system.file("extdata", package="RAIDS"),
#'                  "demoPruning")
#' dir.create(data.dir.pruning, showWarnings=FALSE,
#'                  recursive=FALSE, mode="0777")
#' file.copy(file.path(data.dir, "ex1_demo.gds"),
#'                  file.path(data.dir.pruning, "ex1.gds"))
#'
#' ## Open 1KG file
#' gds1KG <- snpgdsOpen(fileGDS)
#'
#' ## Compute the list of pruned SNVs for a specific profile 'ex1'
#' ## and save it in the Profile GDS file 'ex1.gds'
#' pruningSample(gds=gds1KG, currentProfile=c("ex1"),
#'               studyID = studyDF$study.id, pathProfileGDS=data.dir.pruning)
#'
#' ## Close the 1KG GDS file (it is important to always close the GDS files)
#' closefn.gds(gds1KG)
#'
#' ## Check content of Profile GDS file
#' ## The 'pruned.study' entry should be present
#' content <- openfn.gds(file.path(data.dir.pruning, "ex1.gds"))
#' content
#'
#' ## Close the Profile GDS file (it is important to always close the GDS files)
#' closefn.gds(content)
#'
#' ## Unlink Profile GDS file (created for demo purpose)
#' unlink(file.path(data.dir.pruning, "ex1.gds"))
#' unlink(data.dir.pruning)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
pruningSample <- function(gds, method=c("corr", "r", "dprime", "composite"),
                            currentProfile,
                            studyID,
                            listSNP=NULL,
                            slideWindowMaxBP=500000L,
                            thresholdLD=sqrt(0.1),
                            np=1L,
                            verbose=FALSE,
                            chr=NULL,
                            superPopMinAF=NULL,
                            keepPrunedGDS=TRUE,
                            pathProfileGDS=NULL,
                            keepFile=FALSE,
                            pathPrunedGDS=".",
                            outPrefix="pruned") {

    ## Validate input parameters
    validatePruningSample(gds=gds, method=method, currentProfile=currentProfile,
            studyID=studyID, listSNP=listSNP, slideWindowMaxBP=slideWindowMaxBP,
            thresholdLD=thresholdLD, np=np, verbose=verbose, chr=chr,
            superPopMinAF=superPopMinAF, keepPrunedGDS=keepPrunedGDS,
            pathProfileGDS=pathProfileGDS, keepFile=keepFile,
            pathPrunedGDS=pathPrunedGDS, outPrefix=outPrefix)

    ## Matches a character method against a table of candidate values
    method <- match.arg(method, several.ok=FALSE)

    ## Profile GDS file name
    fileGDSSample <- file.path(pathProfileGDS, paste0(currentProfile, ".gds"))

    ## The Profile GDS file must exists
    if (!(file.exists(fileGDSSample))) {
        stop("The Profile GDS file \'", fileGDSSample, " does not exist.")
    }

    filePruned <- file.path(pathPrunedGDS, paste0(outPrefix, ".rds"))
    fileObj <- file.path(pathPrunedGDS, paste0(outPrefix, ".Obj.rds"))

    snp.id <- read.gdsn(node=index.gdsn(gds, "snp.id"))

    sample.id <- read.gdsn(node=index.gdsn(gds, "sample.id"))

    ## Open the GDS Sample file
    gdsSample <- openfn.gds(filename=fileGDSSample)

    ## Extract all study information from the GDS Sample file
    study.annot <- read.gdsn(node=index.gdsn(gdsSample, "study.annot"))

    ## Select study information associated to the current profile
    posSample <- which(study.annot$data.id == currentProfile &
                            study.annot$study.id == studyID)

    ## Check that the information is found for the specified profile and study
    if(length(posSample) != 1) {
        closefn.gds(gdsSample)
        stop("In pruningSample the profile \'", currentProfile,
                "\' doesn't exists for the study \'", studyID, "\'\n")
    }

    ## Get the SNV genotype information for the current profile
    g <- read.gdsn(index.gdsn(gdsSample, "geno.ref"),
                    start=c(1, posSample), count=c(-1,1))

    ## Close the Profile GDS file
    closefn.gds(gdsSample)

    listGeno <- which(g != 3)
    rm(g)

    listKeepPos <- listGeno

    ## Select SNVs based on the chromosome
    if(!is.null(chr)) {
        snpCHR <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
        listKeepPos <- intersect(which(snpCHR == chr), listKeepPos)
    }

    ## Select SNVs based on the minimum allele frequency in the populations
    if(!is.null(superPopMinAF)) {
        listTMP <- NULL
        for(sp in c("EAS", "EUR", "AFR", "AMR", "SAS")) {
            snpAF <- read.gdsn(index.gdsn(gds, paste0("snp.", sp, "_AF")))
            listTMP <- union(listTMP,
                which(snpAF >= superPopMinAF & snpAF <= 1 - superPopMinAF))
        }
        listKeepPos <- intersect(listTMP, listKeepPos)
    }

    if(length(listKeepPos) == 0) {
        stop("In pruningSample, the sample ", currentProfile,
                " doesn't have SNPs after filters\n")
    }
    listKeep <- snp.id[listKeepPos]

    sample.ref <- read.gdsn(index.gdsn(gds, "sample.ref"))
    listSamples <- sample.id[which(sample.ref == 1)]

    ## Use a LD analysis to generate a subset of SNPs
    snpset <- runLDPruning(gds=gds, method=method, listSamples=listSamples,
                listKeep=listKeep, slideWindowMaxBP=slideWindowMaxBP,
                thresholdLD=thresholdLD, np=np, verbose=verbose)

    pruned <- unlist(snpset, use.names=FALSE)

    ## When TRUE, generate 2 RDS file with the pruned SNVs information
    if(keepFile) {
        saveRDS(pruned, filePruned)
        saveRDS(snpset, fileObj)
    }

    ## When TRUE, add the pruned SNvs information to the Profile GDS file
    if(keepPrunedGDS) {
        gdsSample <- openfn.gds(filename=fileGDSSample, readonly=FALSE)
        addGDSStudyPruning(gds=gdsSample, pruned=pruned)
        closefn.gds(gdsfile=gdsSample)
    }

    return(0L)
}


#' @title Add the genotype information for the list of pruned SNVs
#' into the Profile GDS file
#'
#' @description The function extracts the information about the pruned SNVs
#' from the 1KG GDS file and adds entries related to the pruned SNVs in
#' the Profile GDS file.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param fileProfileGDS a \code{character} string representing the path and
#' file name of the Profile GDS file. The Profile GDS file must exist.
#'
#' @param currentProfile a \code{character} string corresponding to the sample
#' identifier associated to the current list of pruned SNVs.
#'
#' @param studyID a \code{character} string corresponding to the study
#' identifier associated to the current list of pruned SNVs.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(data.dir, "ex1_good_small_1KG_GDS.gds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "studyID", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyDF <- data.frame(study.id="MYDATA",
#'                         study.desc="Description",
#'                         study.platform="PLATFORM",
#'                         stringsAsFactors=FALSE)
#'
#' ## Copy the Profile GDS file demo that has been pruned
#' ## into a test directory (deleted after the example has been run)
#' data.dir.genotype <- file.path(system.file("extdata", package="RAIDS"),
#'                  "demoAddGenotype")
#' dir.create(data.dir.genotype, showWarnings=FALSE,
#'                  recursive=FALSE, mode="0777")
#' file.copy(file.path(data.dir, "ex1_demo_with_pruning.gds"),
#'                  file.path(data.dir.genotype, "ex1.gds"))
#'
#' ## Open 1KG file
#' gds1KG <- snpgdsOpen(fileGDS)
#'
#' ## Compute the list of pruned SNVs for a specific profile 'ex1'
#' ## and save it in the Profile GDS file 'ex1.gds'
#' add1KG2SampleGDS(gds=gds1KG,
#'          fileProfileGDS=file.path(data.dir.genotype, "ex1.gds"),
#'          currentProfile=c("ex1"),
#'          studyID=studyDF$study.id)
#'
#' ## Close the 1KG GDS file (it is important to always close the GDS files)
#' closefn.gds(gds1KG)
#'
#' ## Check content of Profile GDS file
#' ## The 'pruned.study' entry should be present
#' content <- openfn.gds(file.path(data.dir.genotype, "ex1.gds"))
#' content
#'
#' ## Close the Profile GDS file (it is important to always close the GDS files)
#' closefn.gds(content)
#'
#' ## Unlink Profile GDS file (created for demo purpose)
#' unlink(file.path(data.dir.genotype, "ex1.gds"))
#' unlink(data.dir.genotype)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn objdesp.gdsn
#' @encoding UTF-8
#' @export
add1KG2SampleGDS <- function(gds, fileProfileGDS, currentProfile,
                                studyID) {

    ## Validate inputs
    validateAdd1KG2SampleGDS(gds=gds, gdsProfileFile=fileProfileGDS,
            currentProfile=currentProfile, studyID=studyID)

    ## Open Profile GDS file
    gdsSample <- openfn.gds(fileProfileGDS, readonly=FALSE)

    ## Extract needed information from 1KG GDS file
    snp.id <- read.gdsn(index.gdsn(gds,"snp.id"))

    ## Extract list of pruned SNVs from the GDS Sample file
    pruned <- read.gdsn(index.gdsn(gdsSample, "pruned.study"))

    listSNP <- which(snp.id %in% pruned)
    listRef <- which(read.gdsn(index.gdsn(gds, "sample.ref")) == 1)
    sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))

    snp.chromosome <- read.gdsn(index.gdsn(gds,"snp.chromosome"))[listSNP]
    snp.position <-  read.gdsn(index.gdsn(gds,"snp.position"))[listSNP]

    add.gdsn(gdsSample, "sample.id", c(sample.id[listRef], currentProfile))

    add.gdsn(gdsSample, "snp.id", snp.id[listSNP])
    add.gdsn(gdsSample, "snp.chromosome", snp.chromosome)
    add.gdsn(gdsSample, "snp.position", snp.position)
    # snp.index is the index of the snp pruned in snp.id fro 1KG gds
    add.gdsn(gdsSample, "snp.index", listSNP)

    var.geno <- NULL

    j <- 1
    for(i in listRef) {
        g <- read.gdsn(index.gdsn(gds, "genotype"), start=c(1,i),
                            count = c(-1,1))[listSNP]

        if(! ("genotype" %in% ls.gdsn(gdsSample))){
            var.geno <- add.gdsn(gdsSample, "genotype",
                                    valdim=c(length(listSNP), 1),
                                    g, storage="bit2")

        }else {
            if(is.null(var.geno)) {
                var.geno <- index.gdsn(gdsSample, "genotype")
            }
            append.gdsn(var.geno, g)
        }
        if(j %% 5 == 0) {
            sync.gds(gdsSample)
        }
        j <- j + 1
    }

    # add.gdsn(gdsSample, "SamplePos", objdesp.gdsn(index.gdsn(gdsSample,
    #  "genotype"))$dim[2] + 1,
    #          storage="int32")
    study.annot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    posCur <- which(study.annot$data.id == currentProfile &
                            study.annot$study.id == studyID)

    g <- read.gdsn(index.gdsn(gdsSample, "geno.ref"), start=c(1, posCur),
                        count=c(-1, 1))[listSNP]
    append.gdsn(var.geno, g)

    add.gdsn(gdsSample, "lap",
        rep(0.5, objdesp.gdsn(index.gdsn(gdsSample, "genotype"))$dim[1]),
        storage="packedreal8")

    ## Close the GDS Sample file
    closefn.gds(gdsSample)

    return(0L)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory that contains the Profile GDS files. The directory must
#' exist.
#'
#' @param pathGeno a \code{character} string representing the path to
#' the directory that contains TODO
#'
#' @param fileLSNP TODO
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{FALSE}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", 'RAIDS')
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
addPhase1KG2SampleGDSFromFile <- function(gds, pathProfileGDS, pathGeno,
                                            fileLSNP, verbose=FALSE) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gds, name="gds")

    ## Verbose must be a logical
    if (!is.logical(verbose)) {
        stop("The \'verbose\' parameter must be a logical (TRUE or FALSE).")
    }

    listGDSSample <- dir(pathProfileGDS, pattern = ".+.gds")

    indexAll <- NULL
    for(fileProfileGDS in listGDSSample) {
        gdsSample <- openfn.gds(filename=file.path(pathProfileGDS,
                                                        fileProfileGDS))

        snp.index <- read.gdsn(node=index.gdsn(node=gdsSample, "snp.index"))

        indexAll <- union(indexAll, snp.index)
        closefn.gds(gdsfile=gdsSample)
    }

    gdsSample <- createfn.gds(file.path(pathProfileGDS, "phase1KG.gds"))
    indexAll <- indexAll[order(indexAll)]
    snp.id <- read.gdsn(index.gdsn(gds,"snp.id"))[indexAll]
    add.gdsn(gdsSample, "snp.id", snp.id)
    add.gdsn(gdsSample, "snp.index", indexAll)
    listRef <- which(read.gdsn(index.gdsn(gds, "sample.ref"))==1)
    listSample <- read.gdsn(index.gdsn(gds, "sample.id"))[listRef]
    listSNP <- readRDS(file=fileLSNP)
    i<-1
    for(sample1KG in listSample){
        if(verbose) { message("P ", i, " ", Sys.time()) }
        i <- i + 1
        file1KG <- file.path(pathGeno, paste0(sample1KG,".csv.bz2"))
        matSample <- read.csv2(file=file1KG, row.names=NULL)
        matSample <- matSample[listSNP[indexAll],, drop=FALSE]
        matSample <- matrix(as.numeric(unlist(strsplit(matSample[,1],
                                                        "\\|"))), nrow=2)[1,]
        var.phase <- NULL
        if (!("phase" %in% ls.gdsn(gdsSample))) {
            var.phase <- add.gdsn(gdsSample, "phase",
                                    valdim=c(length(indexAll), 1),
                                    matSample, storage="bit2")
        } else {
            if (is.null(var.phase)) {
                var.phase <- index.gdsn(node=gdsSample, "phase")
            }
            append.gdsn(node=var.phase, val=matSample)
        }
    }

    closefn.gds(gdsfile=gdsSample)

    ## Success
    return(0L)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param gdsPhase TODO
#'
#' @param pathProfileGDS the path of an object of class \code{gds} related to
#' the sample
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{TRUE}.
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", "RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
addPhase1KG2SampleGDSFromGDS <- function(gds, gdsPhase, pathProfileGDS,
                                            verbose=FALSE) {

    listGDSSample <- dir(pathProfileGDS, pattern = ".+.gds")

    indexAll <- NULL
    for(fileProfileGDS in listGDSSample){
        gdsSample <- openfn.gds(file.path(pathProfileGDS, fileProfileGDS))

        snp.index <- read.gdsn(index.gdsn(gdsSample,"snp.index"))

        indexAll <- union(indexAll, snp.index)
        closefn.gds(gdsSample)
    }

    gdsSample <- createfn.gds(file.path(pathProfileGDS, "phase1KG.gds"))
    indexAll <- indexAll[order(indexAll)]
    snp.id <- read.gdsn(index.gdsn(gds,"snp.id"))[indexAll]
    add.gdsn(gdsSample, "snp.id", snp.id)
    add.gdsn(gdsSample, "snp.index", indexAll)
    listRef <- which(read.gdsn(index.gdsn(gds, "sample.ref"))==1)
    listSample <- read.gdsn(index.gdsn(gds, "sample.id"))[listRef]
    #listSNP <- readRDS(fileLSNP)
    i<-1
    for(sample1KG in listSample){
        if(verbose) { message("P ", i, " ", Sys.time()) }

        matSample <- read.gdsn(index.gdsn(gdsPhase, "phase"),
                            start=c(1, listRef[i]), count=c(-1,1))[indexAll]
        i<-i+1

        var.phase <- NULL
        if (! ("phase" %in% ls.gdsn(gdsSample))) {
            var.phase <- add.gdsn(gdsSample, "phase",
                                    valdim=c(length(indexAll), 1),
                                    matSample, storage="bit2")

        } else {
            if (is.null(var.phase)) {
                var.phase <- index.gdsn(gdsSample, "phase")
            }
            append.gdsn(var.phase, matSample)
        }
    }

    closefn.gds(gdsSample)

    ## Successful
    return(0L)
}


#' @title Compute principal component axes (PCA) on pruned SNV with the
#' reference samples
#'
#' @description This function compute the PCA on pruned SNV with the
#' reference samples
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param listRef a \code{vector} of string representing the
#' identifier of the profiles in the reference (unrelated).
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @param verbose a \code{logical} indicating if the PCA functions should be
#' verbose. Default: \code{FALSE}.
#'
#' @return listPCA  a \code{list} containing two objects
#' pca.unrel -> \code{snpgdsPCAClass}
#' and a snp.load -> \code{snpgdsPCASNPLoading}
#'
#' @details
#'
#' More information about the method used to calculate the patient eigenvectors
#' can be found at the Bioconductor SNPRelate website:
#' https://bioconductor.org/packages/SNPRelate/
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", "RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASNPLoading
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
computePrunedPCARef <- function(gds, listRef, np=1L, verbose=FALSE) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gds, name="gds")

    ## Validate that np is a single positive integer
    if(! (isSingleNumber(np) && np > 0)) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    if(! is.logical(verbose)) {
        stop("The \'verbose\' parameter must be logical (TRUE or FALSE).")
    }

    listPCA <- list()

    listPruned <- read.gdsn(index.gdsn(gds, "pruned.study"))

    ## Calculate the eigenvectors using the specified SNP loadings for
    ## the reference profiles
    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsobj=gds,
                                            sample.id=listRef,
                                            snp.id=listPruned,
                                            num.thread=np,
                                            verbose=verbose)

    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(pcaobj=listPCA[["pca.unrel"]],
                                                    gdsobj=gds,
                                                    num.thread=np,
                                                    verbose=verbose)
    return(listPCA)
}



#' @title Project profile onto existing principal component axes (PCA)
#'
#' @description This function calculates the profile eigenvectors using
#' the specified SNP loadings.
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, an
#' opened Profile GDS file.
#'
#' @param listPCA  a \code{list} containing two entries:
#' \itemize{
#' \item{pca.unrel} {\code{snpgdsPCAClass} object}
#' \item{snp.load} {\code{snpgdsPCASNPLoading} object}
#' }
#'
#' @param currentProfile a \code{character} string representing the
#' identifiant of the profile to be projected in the PCA.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @param verbose a \code{logical} passed to the PCA function.
#' Default: \code{FALSE}.
#'
#' @return a \code{snpgdsPCAClass} object, a \code{list} that contains:
#' \itemize{
#'    \item{sample.id} {the sample ids used in the analysis}
#'    \item{snp.id} {the SNP ids used in the analysis}
#'    \item{eigenvalues} {eigenvalues}
#'    \item{eigenvect} {eigenvactors, “# of samples” x “eigen.cnt”}
#'    \item{TraceXTX} {the trace of the genetic covariance matrix}
#'    \item{Bayesian} {whether use bayerisan normalization}
#'}
#'
#' @details
#'
#' More information about the method used to calculate the patient eigenvectors
#' can be found at the Bioconductor SNPRelate website:
#' https://bioconductor.org/packages/SNPRelate/
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", "RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom SNPRelate snpgdsPCASampLoading
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
projectSample2PCA <- function(gds, listPCA, currentProfile, np=1L,
                                verbose=FALSE) {


    ## Validate that currentProfile is a character string
    if(! is.character(currentProfile)) {
        stop("The \'currentProfile\' parameter must be a character string.")
    }

    ## Validate that np is a single positive integer
    if(! (isSingleNumber(np) && np > 0)) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    if(! is.logical(verbose)) {
        stop("The \'verbose\' parameter must be logical (TRUE or FALSE).")
    }

    ## Calculate the sample eigenvectors using the specified SNP loadings
    samplePCA <- snpgdsPCASampLoading(listPCA[["snp.load"]],
                                gdsobj=gds, sample.id=currentProfile,
                                num.thread=1, verbose=verbose)

    return(samplePCA)
}


#' @title Estimate the allelic fraction of the pruned SNVs for a specific
#' sample
#'
#' @description The function estimates the allelic fraction of the pruned
#' SNVs for a specific sample and add the information to the associated
#' GDS Sample file. The allelic fraction estimation method is adapted to
#' the type of study (DNA or RNA).
#'
#' @param gds an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the GDS Sample file.
#'
#' @param currentProfile a \code{character} string corresponding to
#' the sample identifier as used in \code{\link{pruningSample}} function.
#'
#' @param studyID a \code{character} string corresponding to the name of
#' the study as
#' used in \code{\link{pruningSample}} function.
#'
#' @param chrInfo a \code{vector} of \code{integer} values representing
#' the length of the chromosomes. See 'details' section.
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "DNA" and "RNA". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"DNA"}.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#'
#' @param minProb a single \code{numeric} between 0 and 1 representing TODO.
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
#' the window used to compute an empty box. Default: \code{9}.
#'
#' @param cutOffAR a single \code{numeric} representing the cutoff, in
#' log score, that the SNVs in a gene are allelic fraction different 0.5
#' Default: \code{3}.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file.
#'  **This parameter is RNA specific.**
#' Default: \code{NULL}.
#'
#' @param block.id a \code{character} string corresponding to the block
#' identifier in \code{gdsRefAnnot}. **This parameter is RNA specific.**
#' Default: \code{NULL}
#'
#' @return The integer \code{0L} when successful.
#'
#' @details
#'
#' The `chrInfo` parameter contains the length of the chromosomes. The
#' length of the chromosomes can be obtain through the
#' \code{\link[BSgenome]{BSgenome-class}}
#' library.
#'
#' As example:
#'
#' ```
#'
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' chrInfo <- integer(25L)
#'
#' for(i in seq_len(22L)){ chrInfo[i] <- length(Hsapiens[[paste0("chr", i)]])}
#'
#' chrInfo[23] <- length(Hsapiens[["chrX"]])
#'
#' chrInfo[24] <- length(Hsapiens[["chrY"]])
#'
#' chrInfo[25] <- length(Hsapiens[["chrM"]])
#'
#' ```
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", "RAIDS")
#'
#' ## Chromosome length information
#' ## chr23 is chrX, chr24 is chrY and chrM is 25
#' chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'     181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
#'     135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
#'     83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
#'     156040895L, 57227415L,  16569L)
#'
#' ## A formal way to get the chormosome length information
#' ## library(BSgenome.Hsapiens.UCSC.hg38)
#' ## chrInfo <- integer(25L)
#' ## for(i in seq_len(22L)){ chrInfo[i] <-
#' ##                          length(Hsapiens[[paste0("chr", i)]])}
#' ## chrInfo[23] <- length(Hsapiens[["chrX"]])
#' ## chrInfo[24] <- length(Hsapiens[["chrY"]])
#' ## chrInfo[25] <- length(Hsapiens[["chrM"]])
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
estimateAllelicFraction <- function(gds, gdsSample, currentProfile, studyID,
    chrInfo, studyType=c("DNA", "RNA"), minCov=10L, minProb=0.999,
    eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=9, cutOffAR=3,
    gdsRefAnnot=NULL, block.id=NULL) {

    ## Validate input parameters
    validateEstimateAllelicFraction(gds=gds, gdsSample=gdsSample,
        currentProfile=currentProfile, studyID=studyID, chrInfo=chrInfo,
        studyType=studyType, minCov=minCov, minProb=minProb, eProb=eProb,
        cutOffLOH=cutOffLOH, cutOffHomoScore=cutOffHomoScore, wAR=wAR,
        cutOffAR=cutOffAR, gdsRefAnnot=gdsRefAnnot, block.id=block.id)

    ## Set study type
    studyType <- match.arg(studyType)

    snp.pos <- NULL

    ## The type of study affects the allelic fraction estimation
    if(studyType == "DNA") {
        snp.pos <- computeAllelicFractionDNA(gds=gds, gdsSample=gdsSample,
                        currentProfile=currentProfile, studyID=studyID,
                        chrInfo=chrInfo, minCov=minCov, minProb=minProb,
                        eProb=eProb, cutOffLOH=cutOffLOH,
                        cutOffHomoScore=cutOffHomoScore, wAR=wAR)

    } else if(studyType == "RNA") {
        snp.pos <- computeAllelicFractionRNA(gds=gds, gdsSample=gdsSample,
                        gdsRefAnnot=gdsRefAnnot, currentProfile=currentProfile,
                        studyID=studyID, block.id=block.id, chrInfo=chrInfo,
                        minCov=minCov, minProb=minProb, eProb=eProb,
                        cutOffLOH=cutOffLOH, cutOffAR=cutOffAR)
    }

    snp.pos$seg <- rep(0, nrow(snp.pos))
    k <- 1
    # Find segment with same lap
    for(chr in seq_len(22)) {
        snpChr <- snp.pos[snp.pos$snp.chr == chr, ]
        tmp <- c(0, abs(snpChr[2:nrow(snpChr), "lap"] -
                            snpChr[seq_len(nrow(snpChr)- 1),  "lap"]) > 1e-3)
        snp.pos$seg[snp.pos$snp.chr == chr] <- cumsum(tmp) + k
        k <- max(snp.pos$seg[snp.pos$snp.chr == chr]) + 1
    }

    ## Save information into the "lap" node in the GDS Sample file
    ## Save information into the "segment" node in the GDS Sample file
    ## Suppose we keep only the pruned SNVs
    addUpdateLap(gdsSample, snp.pos$lap[which(snp.pos$pruned == TRUE)])
    addUpdateSegment(gdsSample, snp.pos$seg[which(snp.pos$pruned == TRUE)])

    # Successful
    return(0L)
}


#' @title Append information about the 1KG samples into
#' the GDS Sample file
#'
#' @description The information about the samples present in the 1KG GDS file
#' is added into the GDS Sample file. Only the information about the
#' unrelated samples
#' from the 1OOO Genome Study are copied into the GDS Sample file. The
#' information is only added to the GDS Sample file when the 1KG Study is not
#' already present in the GDS Sample file. The sample information for all
#' selected samples is appended to the GDS Sample file "study.annot" node.
#' The study information is appended to the GDS Sample file "study.list" node.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param fileProfileGDS a \code{character} string representing the path and
#' file name of the GDS Sample file. The GDS Sample file must exist.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @export
addStudy1Kg <- function(gds, fileProfileGDS) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds, "gds")

    ## The fileProfileGDS must be a character string and the file must exists
    if(!(is.character(fileProfileGDS) && (file.exists(fileProfileGDS)))) {
        stop("The \'fileProfileGDS\' must be a character string representing ",
                "the GDS Sample file. The file must exist.")
    }

    ## Open GDS Sample file
    gdsSample <- openfn.gds(filename=fileProfileGDS, readonly=FALSE)

    ## Extract study information from GDS Sample file
    snp.study <- read.gdsn(index.gdsn(node=gdsSample, "study.list"))

    ## When the 1KG Study is not already present in the GDS Sample file
    if(length(which(snp.study$study.id == "Ref.1KG")) == 0) {

        ## Extract information about all samples from 1KG that are unrelated
        ## and can be used in the ancestry analysis
        sample.ref <- read.gdsn(index.gdsn(node=gds, "sample.ref"))
        sample.id <- read.gdsn(index.gdsn(node=gds,
                                        "sample.id"))[which(sample.ref == 1)]

        ## Create study information for the 1KG Study
        study.list <- data.frame(study.id="Ref.1KG",
                        study.desc="Unrelated samples from 1000 Genomes",
                        study.platform="GRCh38 1000 genotypes",
                        stringsAsFactors=FALSE)

        ## Create the pedigree information  for the 1KG samples
        ped1KG <- data.frame(Name.ID=sample.id, Case.ID=sample.id,
                            Sample.Type=rep("Reference", length(sample.id)),
                            Diagnosis=rep("Reference", length(sample.id)),
                            Source=rep("IGSR", length(sample.id)),
                            stringsAsFactors=FALSE)

        ## Row names must be the sample identifiers
        rownames(ped1KG) <- ped1KG$Name.ID

        ## Add the information about the 1KG samples into the
        ## GDS sample
        addStudyGDSSample(gds=gdsSample, pedDF=ped1KG, batch=1,
                            listSamples=NULL, studyDF=study.list)

        sync.gds(gdsSample)
    }

    ## Close GDS Sample file (important)
    closefn.gds(gdsSample)

    ## Return success
    return(0L)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param gdsSample an object of class \code{gds} opened related to
#' the sample
#'
#' @param listPCA TODO
#'
#' @param sampleRef a \code{vector} of sample.id from 1KG with one entry for
#' each synthetic to project.
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier.
#' The study identifier must be present in the GDS Sample file.
#'
#' @return a \code{list} containing 3 entries:
#' \itemize{
#' \item{sample.id} { TODO }
#' \item{eigenvector.ref} { TODO }
#' \item{eigenvecto} { TODO }
#' }
#'
#' @examples
#'
#' # TODO
#' gds <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCASNPLoading snpgdsPCASampLoading
#' @encoding UTF-8
#' @export
computePCAMultiSynthetic <- function(gdsSample, listPCA,
                                        sampleRef, studyIDSyn) {

    if(length(sampleRef) < 1) {
        stop("Number of sample in study.annot not equal to 1\n")
    }

    study.annot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))
    study.annot <- study.annot[which(study.annot$study.id == studyIDSyn &
                                        study.annot$case.id %in% sampleRef),]


    ## SNP loading in principal component analysis
    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                                    gdsobj=gdsSample,
                                                    num.thread=1, verbose=TRUE)

    ## Project samples onto existing principal component axes
    listPCA[["samp.load"]] <- snpgdsPCASampLoading(listPCA[["snp.load"]],
                                                gdsobj=gdsSample,
                                                sample.id=study.annot$data.id,
                                                num.thread=1L, verbose=TRUE)

    rownames(listPCA[["pca.unrel"]]$eigenvect) <-
                                            listPCA[["pca.unrel"]]$sample.id

    rownames(listPCA[["samp.load"]]$eigenvect) <-
                                            listPCA[["samp.load"]]$sample.id

    listRes <- list(sample.id=listPCA[["samp.load"]]$sample.id,
                        eigenvector.ref=listPCA[["pca.unrel"]]$eigenvect,
                        eigenvector=listPCA[["samp.load"]]$eigenvect)

    return(listRes)
}


#' @title Compute Principal Component Analysis (PCA) on SNV genotype for
#' data from a sample
#'
#' @description TODO
#'
#' @param gdsSample an object of class \link[gdsfmt]{gds.class},
#' a GDS Sample file.
#'
#' @param name.id a single \code{character} string representing the sample
#' identifier.
#'
#' @param studyIDRef a single \code{character} string representing the
#' study identifier.
#'
#' @param np a single positive \code{integer} representing the number of CPU
#' that will be used. Default: \code{1L}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA}
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
#' Default: \code{NaN}.
#'
#' @references
#'
#' Galinsky KJ, Bhatia G, Loh PR, Georgiev S, Mukherjee S, Patterson NJ,
#' Price AL. Fast Principal-Component Analysis Reveals Convergent Evolution
#' of ADH1B in Europe and East Asia. Am J Hum Genet. 2016 Mar 3;98(3):456-72.
#' doi: 10.1016/j.ajhg.2015.12.022. Epub 2016 Feb 25.
#'
#' @return a \code{list} containing 3 entries:
#' \itemize{
#' \item{sample.id} { TODO }
#' \item{eigenvector.ref} { TODO }
#' \item{eigenvector} { TODO }
#' }
#'
#' @examples
#'
#' # TODO
#' gds <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASNPLoading snpgdsPCASampLoading
#' @encoding UTF-8
#' @export
computePCARefSample <- function(gdsSample, name.id, studyIDRef="Ref.1KG",
                            np=1L, algorithm=c("exact","randomized"),
                            eigen.cnt=32L, missing.rate=NaN) {

    ## Validate parameters
    validateComputePCARefSample(gdsSample=gdsSample, name.id=name.id,
        studyIDRef=studyIDRef, np=np, algorithm=algorithm,
        eigen.cnt=eigen.cnt, missing.rate=missing.rate)

    ## Set algorithm
    algorithm <- match.arg(algorithm)

    sample.id <- read.gdsn(index.gdsn(gdsSample, "sample.id"))

    sample.pos <- which(sample.id == name.id)

    study.annot.all <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    sample.Unrel <- study.annot.all[which(study.annot.all$study.id ==
                                                studyIDRef), "data.id"]

    listPCA <- list()

    listPCA[["pruned"]] <- read.gdsn(index.gdsn(gdsSample, "pruned.study"))

    ## Calculate the eigenvectors and eigenvalues for PCA
    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsSample, sample.id=sample.Unrel,
                                snp.id=listPCA[["pruned"]], num.thread=np,
                                algorithm=algorithm, eigen.cnt=eigen.cnt,
                                missing.rate=missing.rate, verbose=TRUE)

    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                gdsobj=gdsSample, num.thread=np, verbose=TRUE)

    listPCA[["samp.load"]] <- snpgdsPCASampLoading(listPCA[["snp.load"]],
                                gdsobj=gdsSample,
                                sample.id=sample.id[sample.pos],
                                num.thread=np, verbose=TRUE)

    rownames(listPCA[["pca.unrel"]]$eigenvect) <-
                                        listPCA[["pca.unrel"]]$sample.id
    rownames(listPCA[["samp.load"]]$eigenvect) <-
                                        listPCA[["samp.load"]]$sample.id

    listRes <- list(sample.id=sample.id[sample.pos],
                        eigenvector.ref=listPCA[["pca.unrel"]]$eigenvect,
                        eigenvector=listPCA[["samp.load"]]$eigenvect)

    return(listRes)
}


#' @title Run a k-nearest neighbors analysis on a subset of the
#' synthetic dataset
#'
#' @description TODO
#'
#' @param gdsSample an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, the
#' GDS Sample file.
#'
#' @param listEigenvector TODO see return of computePCAsynthetic
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param studyIDSyn a the studyID of the synthetic data
#'
#' @param spRef TODO
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset. Default: \code{"SuperPop"}.
#'
#' @param kList  a \code{vector} of \code{integer} representing  the list of
#' values tested for the  K parameter. The K parameter represents the
#' number of neighbors used in the K-nearest neighbors analysis. If
#' \code{NULL}, the value \code{seq(2, 15, 1)} is assigned.
#' Default: \code{seq(2, 15, 1)}.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  D parameter. The D parameter represents the
#' number of dimensions used in the PCA analysis.  If \code{NULL},
#' the value \code{seq(2, 15, 1)} is assigned.
#' Default: \code{seq(2, 15, 1)}.
#'
#' @return a \code{list} TODO with the sample.id and eigenvectors
#' and a table with KNN call for different K and PCA dimensions.
#' \itemize{
#' \item{sample.id} {TODO}
#' \item{sample1Kg} {TODO}
#' \item{sp} {TODO}
#' \item{matKNN} {TODO}
#' }
#'
#' @examples
#'
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom class knn
#' @encoding UTF-8
#' @export
computeKNNRefSynthetic <- function(gdsSample, listEigenvector,
                                    listCatPop, studyIDSyn,
                                    spRef, fieldPopInfAnc="SuperPop",
                                    kList=seq(2, 15, 1),
                                    pcaList=seq(2, 15, 1)) {

    ## The number of rows in study.annot must be one.
    # if(nrow(study.annot) < 1) {
    #     stop("Number of samples in study.annot not equal to 1\n")
    # }

    ## Assign default value if kList is NULL
    if(is.null(kList)) {
        kList <- seq(2, 15, 1) #c(seq_len(14), seq(15,100, by=5))
    }

    ## Assign default value if pcaList is NULL
    if(is.null(pcaList)) {
        pcaList <- seq(2, 15, 1)
    }

    ## Get study information from the GDS Sample file
    study.annot.all <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    study.annot <- study.annot.all[which(study.annot.all$study.id ==
                                                studyIDSyn &
                                study.annot.all$data.id %in%
                                                listEigenvector$sample.id), ]

    listMat <- list()
    for(i in seq_len(length(listEigenvector$sample.id))){
        resMat <- data.frame(sample.id=rep(listEigenvector$sample.id[i],
                                            length(pcaList) * length(kList)),
                                D=rep(0,length(pcaList) * length(kList)),
                                K=rep(0,length(pcaList) * length(kList)),
                        # SuperPop=character(length(pcaList) * length(kList)),
                                stringsAsFactors=FALSE)
        resMat[[fieldPopInfAnc]] <- character(length(pcaList) * length(kList))

        eigenvect <- rbind(listEigenvector$eigenvector.ref,
                                listEigenvector$eigenvector[i,,drop=FALSE])

        totR <- 1
        for (pcaD in pcaList) {
            for(kV in  seq_len(length(kList))) {
                dCur <- paste0("d", pcaD)
                kCur <- paste0("k", kList[kV])
                resMat[totR,c("D", "K")] <- c(pcaD, kList[kV])

                pcaND <- eigenvect[ ,seq_len(pcaD)]
                y_pred <-
                    knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                    test=pcaND[rownames(eigenvect)[nrow(eigenvect)],,
                                                                drop=FALSE],
                    cl=factor(spRef[rownames(eigenvect)[-1*nrow(eigenvect)]],
                                        levels=listCatPop, labels=listCatPop),
                    k=kList[kV],
                    prob=FALSE)

                resMat[totR, fieldPopInfAnc] <- listCatPop[as.integer(y_pred)]

                totR <- totR + 1
            } # end k
        } # end pca Dim
        listMat[[i]] <- resMat
    }
    resMat <- do.call(rbind, listMat)

    listKNN <- list(sample.id=listEigenvector$sample.id,
                    sample1Kg=study.annot$case.id,
                    sp=spRef[study.annot$case.id], matKNN=resMat)

    return(listKNN)
}



#' @title Run a k-nearest neighbors analysis on a referenc sample
#'
#' @description TODO
#'
#' @param listEigenvector TODO see return of computePCARefSample
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param spRef TODO
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param pcaList TODO array of the pca dimension possible values
#'
#' @return A \code{list} TODO with the sample.id and eigenvectors
#' and a table with KNN callfor different K and pca dimension.
#'
#' @examples
#'
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom class knn
#' @encoding UTF-8
#' @export
computeKNNRefSample <- function(listEigenvector,
                            listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                            spRef, fieldPopInfAnc="SuperPop",
                            kList=seq(2, 15, 1), pcaList=seq(2, 15, 1)) {

    if(is.null(kList)){
        kList <- seq_len(15)#c(seq_len(14), seq(15,100, by=5))
    }
    if(is.null(pcaList)){
        pcaList <- 2:15
    }
    if(length(listEigenvector$sample.id) != 1) {
        stop("Number of sample in study.annot not equal to 1\n")
    }

    resMat <- data.frame(sample.id=rep(listEigenvector$sample.id,
                                        length(pcaList) * length(kList)),
                            D=rep(0,length(pcaList) * length(kList)),
                            K=rep(0,length(pcaList) * length(kList)),
                        # SuperPop=character(length(pcaList) * length(kList)),
                            stringsAsFactors=FALSE)
    resMat[[fieldPopInfAnc]] <- character(length(pcaList) * length(kList))

    #curPCA <- listPCA.Samples[[sample.id[sample.pos]]]
    eigenvect <- rbind(listEigenvector$eigenvector.ref,
                            listEigenvector$eigenvector)

    totR <- 1
    for(pcaD in pcaList) {
        for(kV in  seq_len(length(kList))) {
            dCur <- paste0("d", pcaD)
            kCur <- paste0("k", kList[kV])
            resMat[totR,c("D", "K")] <- c(pcaD, kList[kV])

            pcaND <- eigenvect[ ,seq_len(pcaD)]
            y_pred <- knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                    test=pcaND[rownames(eigenvect)[nrow(eigenvect)],,
                                drop=FALSE],
                    cl=factor(spRef[rownames(eigenvect)[-1*nrow(eigenvect)]],
                                    levels=listCatPop, labels=listCatPop),
                    k=kList[kV],
                    prob=FALSE)

            resMat[totR, fieldPopInfAnc] <- listCatPop[as.integer(y_pred)]

            totR <- totR + 1
        } # end k
    } # end pca Dim

    listKNN <- list(sample.id=listEigenvector$sample.id, matKNN=resMat)

    return(listKNN)
}

#' @title Run a PCA analysis and a K nearest neighbor analysis on a small set
#' of synthetic data
#'
#' @description The function runs a PCA analysis using 1 synthetic profile
#' from each sub-continental population. The reference samples used to
#' create those synthetic profiles are first removed from the 1KG list
#' of samples that generates the reference PCA. Then, the retained synthetic
#' profiles are projected on the 1KG PCA space. Finally, a K nearest neighbor
#' analysis using a range of K and D values is done.
#'
#' @param gdsSample an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, the
#' GDS Sample file.
#'
#' @param sampleRM a \code{vector} of \code{character} strings representing
#' the identifiers of the 1KG reference samples that should not be used to
#' create the reference PCA.
#'
#' @param spRef TODO
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier.
#' The study identifier must be present in the GDS Sample file.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param fieldPopIn1KG TODO. Default: \code{"superPop"}.
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.  If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA} function
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
#' Default: \code{0.025}.
#'
#' @return a \code{list} TODO with the sample.id and eigenvectors
#' and a table with KNN call for different K and pca dimension.
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
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
computePoolSyntheticAncestryGr <- function(gdsSample, sampleRM, spRef,
                            studyIDSyn, np=1L,
                            listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                            fieldPopIn1KG="superPop",
                            fieldPopInfAnc="SuperPop",
                            kList=seq(2,15,1),
                            pcaList=seq(2,15,1),
                            algorithm=c("exact", "randomized"),
                            eigen.cnt=32L,
                            missing.rate=0.025) {

    ## Assign default value is kList is NULL
    if(is.null(kList)) {
        kList <- seq(2,15,1)
    }

    ## Assign default value is pcaList is NULL
    if(is.null(pcaList)) {
        pcaList <- seq(2,15,1)
    }

    ## Validate the input parameters
    validateComputePoolSyntheticAncestryGr(gdsSample=gdsSample,
        sampleRM=sampleRM, spRef=spRef, studyIDSyn=studyIDSyn,
        np=np, listCatPop=listCatPop, fieldPopIn1KG=fieldPopIn1KG,
        fieldPopInfAnc=fieldPopInfAnc, kList=kList, pcaList=pcaList,
        algorithm=algorithm, eigen.cnt=eigen.cnt, missing.rate=missing.rate)

    ## Set algorithm
    algorithm <- match.arg(algorithm)

    ## Calculate Principal Component Analysis (PCA) on SNV genotype dataset
    pca1KG <- computePCARefRMMulti(gdsSample=gdsSample,
                        sample.ref=names(spRef), listRM=sampleRM, np=np,
                        algorithm=algorithm, eigen.cnt=eigen.cnt,
                        missing.rate=missing.rate)

    resPCA <- computePCAMultiSynthetic(gdsSample=gdsSample, listPCA=pca1KG,
                        sampleRef=sampleRM, studyIDSyn=studyIDSyn)

    ## Calculate the k-nearest neighbor analyses on a subset of the
    ## synthetic dataset
    KNN.synt <- computeKNNRefSynthetic(gdsSample=gdsSample,
                        listEigenvector=resPCA, listCatPop=listCatPop,
                        studyIDSyn=studyIDSyn, spRef=spRef,
                        fieldPopInfAnc=fieldPopInfAnc, kList=kList,
                        pcaList=pcaList)

    return(KNN.synt)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1KG GDS file.
#'
#' @param gdsSample an object of class \code{gds} opened related to
#' the sample
#'
#' @param sample.ana.id TODO
#'
#' @param dataRef a \code{data.frame} TODO
#'
#' @param spRef TODO
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier. The study identifier must be present in the GDS Sample file.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param fieldPopIn1KG TODO
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbors analysis. If
#' \code{NULL}, the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.  If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA}
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
#' Default: \code{0.025}.
#'
#'
#' @return A \code{list} TODO with the sample.id and eigenvectors
#' and a table with KNN callfor different K and pca dimension.
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
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
computePoolSyntheticAncestry <- function(gds, gdsSample,
                                            sample.ana.id,
                                            dataRef, spRef,
                                            studyIDSyn,
                                            np = 1L,
                                            listCatPop = c("EAS", "EUR",
                                                        "AFR", "AMR", "SAS"),
                                            fieldPopIn1KG = "superPop",
                                            fieldPopInfAnc = "SuperPop",
                                            kList = seq(2, 15, 1),
                                            pcaList = seq(2, 15, 1),
                                            algorithm="exact",
                                            eigen.cnt=32L,
                                            missing.rate=0.025) {

    ## TODO Add parameter validation (not all done)

    ## The gds and gdsSample must be objects of class "gds.class"
    validateGDSClass(gds=gds, "gds")
    validateGDSClass(gds=gdsSample, "gdsSample")

    ## The dataRef must be an data.frame object
    if (!is.data.frame(dataRef)) {
        stop("The \'dataRef\' must be a data.frame object.")
    }

    ## The studyID must be a character string
    if (!(is.character(studyIDSyn) && length(studyIDSyn) == 1)) {
        stop("The \'studyIDSyn\' parameter must be a character string.")
    }

    ## The population name in 1KG must be a character string
    if (!(is.character(fieldPopIn1KG) && length(fieldPopIn1KG) == 1)) {
        stop("The \'fieldPopIn1KG\' parameter must be a character string.")
    }



    sampleRM <- splitSelectByPop(dataRef)

    KNN.list <- list()
    for(j in seq_len(nrow(sampleRM))) {
        ## Run a PCA analysis using 1 synthetic profile from each
        ##  sub-continental ancestry
        ## The synthetic profiles are projected on the 1KG PCA space
        ##  (the reference samples used to generate the synthetic profiles are
        ##  removed from this PCA)
        KNN.list[[j]] <- computePoolSyntheticAncestryGr(gdsSample=gdsSample,
                            sampleRM=sampleRM[j,],
                            spRef=spRef, studyIDSyn=studyIDSyn,
                            np=np, listCatPop=listCatPop,
                            fieldPopIn1KG=fieldPopIn1KG,
                            fieldPopInfAnc=fieldPopInfAnc, kList=kList,
                            pcaList=pcaList, algorithm=algorithm,
                            eigen.cnt=eigen.cnt, missing.rate=missing.rate)
    }

    KNN.sample.syn <- do.call(rbind, KNN.list)

    pedSyn <- prepPedSynthetic1KG(gds=gds, gdsSample=gdsSample,
                            studyID=studyIDSyn, popName=fieldPopIn1KG)

    listParaSample <- selParaPCAUpQuartile(KNN.sample.syn, pedSyn,
                                            fieldPopIn1KG, fieldPopInfAnc,
                                            listCatPop)

    listPCASample <- computePCARefSample(gdsSample=gdsSample,
                            name.id=sample.ana.id, studyIDRef="Ref.1KG",
                            np=np, algorithm=algorithm,
                            eigen.cnt=eigen.cnt, missing.rate=missing.rate)

    listKNNSample <- computeKNNSuperPopSample(gdsSample=gdsSample,
                                                sample.ana.id,
                                                spRef)

    return(listKNNSample)
}

#' @title Select the optimal K and D parameters using the synthetic data and
#' infer the ancestry of a specific sample
#'
#' @description TODO
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1KG GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the GDS Sample file.
#'
#' @param listFiles TODO.
#'
#' @param sample.ana.id TODO
#'
#' @param spRef TODO
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier. The study identifier must be present in the GDS Sample file.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param fieldPopIn1KG a \code{character} string representing the name of TODO
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.  If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA}
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
#' Default: \code{NaN}.
#'
#' @return a \code{list} TODO with the sample.id and eigenvectors
#' and a table with KNN callfor different K and pca dimension.
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
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
computeAncestryFromSyntheticFile <- function(gds, gdsSample,
                            listFiles,
                            sample.ana.id,
                            spRef,
                            studyIDSyn,
                            np=1L,
                            listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                            fieldPopIn1KG="superPop",
                            fieldPopInfAnc="SuperPop",
                            kList=seq(2, 15, 1),
                            pcaList=seq(2, 15, 1),
                            algorithm=c("exact", "randomized"),
                            eigen.cnt=32L,
                            missing.rate=NaN) {

    if(is.null(pcaList)) {
        pcaList <- seq(2, 15, 1)
    }

    if (is.null(kList)) {
        kList <- seq(2, 15, 1)
    }

    ## Validate input parameters
    validateComputeAncestryFromSyntheticFile(gds=gds, gdsSample=gdsSample,
        listFiles=listFiles, sample.ana.id=sample.ana.id, spRef=spRef,
        studyIDSyn=studyIDSyn, np=np, listCatPop=listCatPop,
        fieldPopIn1KG=fieldPopIn1KG, fieldPopInfAnc=fieldPopInfAnc, kList=kList,
        pcaList=pcaList, algorithm=algorithm, eigen.cnt=eigen.cnt,
        missing.rate=missing.rate)

    ## Matches a character method against a table of candidate values
    algorithm <- match.arg(algorithm, several.ok=FALSE)

    ## Merge results from PCA run on synthetic data present in RDS files
    KNN.list <- list()
    for(j in seq_len(length(listFiles))) {
        # We have to test if the file exist and format is OK
        KNN.list[[j]] <- readRDS(listFiles[j])
    }
    KNN.sample.syn <- do.call(rbind, KNN.list)

    ## Extract the sample super-population information from the 1KG GDS file
    ## for profiles associated to the specified study in the GDS Sample file
    pedSyn <- prepPedSynthetic1KG(gds=gds, gdsSample=gdsSample,
        studyID=studyIDSyn, popName=fieldPopIn1KG)

    listParaSample <- selParaPCAUpQuartile(matKNN.All=KNN.sample.syn,
        pedCall=pedSyn, refCall=fieldPopIn1KG, predCall=fieldPopInfAnc,
        listCall=listCatPop)

    listPCASample <- computePCARefSample(gdsSample=gdsSample,
        name.id=sample.ana.id, studyIDRef="Ref.1KG", np=np,
        algorithm=algorithm, eigen.cnt=eigen.cnt, missing.rate=missing.rate)

    listKNNSample <- computeKNNRefSample(listEigenvector=listPCASample,
        listCatPop=listCatPop, spRef=spRef,fieldPopInfAnc=fieldPopInfAnc,
        kList=kList, pcaList=pcaList)

    resCall <- listKNNSample$matKNN[
        which(listKNNSample$matKNN$D == listParaSample$D &
                        listKNNSample$matKNN$K == listParaSample$K ) ,]

    res <- list(pcaSample=listPCASample, # PCA of the sample + 1KG
                paraSample=listParaSample, # Result of the parameter selection
                KNNSample=listKNNSample, # KNN for the sample
                Ancestry=resCall # the ancestry call fo rthe sample
                )

    return(res)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{fileNamePED} or \code{pedStudy}
#' can be defined.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "studyID", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Profile files will be created.
#' Default: \code{NULL}.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the VCF output of SNP-pileup for each sample. The
#' SNP-pileup files must be compressed (gz files) and have the name identifiers
#' of the samples. A sample with "Name.ID" identifier would have an
#' associated SNP-pileup file called "Name.ID.txt.gz".
#'
#' @param pathOut a \code{character} string representing the path to
#' the directory where the output files are created.
#'
#' @param fileReferenceGDS  a \code{character} string representing the file
#' name of the 1KG GDS file. The file must exist.
#'
#' @param fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the 1KG GDS annotation file. The file must exist.
#'
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section.
#'
#' @param dataRefSyn a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The \code{data.frame} must have those columns:
#' \itemize{
#' \item{sample.id} { a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group} { a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop} { a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @return The integer \code{0L} when successful. See details section for
#' more information about the generated output files.
#'
#' @details
#'
#' The runExomeAncestry() function generates 3 types of files
#' in the OUTPUT directory.
#' \itemize{
#' \item{Ancestry Inference}{The ancestry inference CSV file
#' (".Ancestry.csv" file)}
#' \item{Inference Informaton}{The inference information RDS file
#' (".infoCall.rds" file)}
#' \item{Synthetic Information}{The parameter information RDS files
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
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' #################################################################
#' ## The path and file name for the PED RDS file
#' ## will the information about the analyzed samples
#' #################################################################
#' filePED <- file.path(dataDir, "example", "pedEx.rds")
#' ped <- readRDS(filePED)
#' head(ped)
#'
#' #################################################################
#' ## The 1KG GDS file and the 1KG SNV Annotation GDS file
#' ## need to be located in the same directory
#' ## Note that the 1KG GDS file used for this example is a
#' ## simplified version and CANNOT be used for any real analysis
#' #################################################################
#' path1KG <- file.path(dataDir, "example", "gdsRef")
#'
#' fileReferenceGDS  <- file.path(path1KG, "ex1kg.gds")
#' fileAnnotGDS <- file.path(path1KG, "exAnnot1kg.gds")
#'
#' #################################################################
#' ## The Sample SNP pileup files (one per sample) need
#' ## to be located in the same directory.
#' #################################################################
#' pathGeno <- file.path(dataDir, "example", "snpPileup")
#'
#' #################################################################
#' ## The path where the Sample GDS files (one per sample)
#' ## will be created need to be specified.
#' #################################################################
#' pathProfileGDS <- file.path(dataDir, "example", "out.tmp")
#'
#' pathOut <- file.path(dataDir, "example", "res.out")
#'
#' #################################################################
#' ## A data frame containing general information about the study
#' ## is also required. The data frame must have
#' ## those 3 columns: "studyID", "study.desc", "study.platform"
#' #################################################################
#' studyDF <- data.frame(study.id="MYDATA",
#'                         study.desc="Description",
#'                         study.platform="PLATFORM",
#'                         stringsAsFactors=FALSE)
#'
#' ####################################################################
#' ## Fix seed to ensure reproducible results
#' ####################################################################
#' set.seed(3043)
#'
#' gds1KG <- snpgdsOpen(fileReferenceGDS)
#' dataRef <- select1KGPop(gds1KG, nbSamples=2L)
#' closefn.gds(gds1KG)
#'
#' ## Chromosome length information
#' ## chr23 is chrX, chr24 is chrY and chrM is 25
#' chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'     181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
#'     135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
#'     83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
#'     156040895L, 57227415L,  16569L)
#'
#' ## A formal way to get the chormosome length information
#' ## library(BSgenome.Hsapiens.UCSC.hg38)
#' ## chrInfo <- integer(25L)
#' ## for(i in seq_len(22L)){ chrInfo[i] <-
#' ##                          length(Hsapiens[[paste0("chr", i)]])}
#' ## chrInfo[23] <- length(Hsapiens[["chrX"]])
#' ## chrInfo[24] <- length(Hsapiens[["chrY"]])
#' ## chrInfo[25] <- length(Hsapiens[["chrM"]])
#'
#' \dontrun{
#' runExomeAncestry(pedStudy=ped, studyDF=studyDF,
#'                     pathProfileGDS=pathProfileGDS,
#'                     pathGeno=pathGeno,
#'                     pathOut=pathOut,
#'                     fileReferenceGDS=fileReferenceGDS,
#'                     fileReferenceAnnotGDS=fileAnnotGDS,
#'                     chrInfo=chrInfo,
#'                     dataRefSyn=dataRef)
#'
#' unlink(pathProfileGDS, recursive=TRUE, force=TRUE)
#' unlink(pathOut, recursive=TRUE, force=TRUE)
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv
#' @encoding UTF-8
#' @export
runExomeAncestry <- function(pedStudy, studyDF, pathProfileGDS,
                    pathGeno, pathOut, fileReferenceGDS, fileReferenceAnnotGDS,
                        chrInfo, dataRefSyn) {

    ## Validate parameters
    validateRunExomeAncestry(pedStudy, studyDF, pathProfileGDS,
        pathGeno=pathGeno, pathOut=pathOut, fileReferenceGDS=fileReferenceGDS,
        fileReferenceAnnotGDS=fileReferenceAnnotGDS, chrInfo=chrInfo,
        dataRefSyn=dataRefSyn)

    listProfiles <- pedStudy[, "Name.ID"]

    createStudy2GDS1KG(pathGeno=pathGeno, pedStudy=pedStudy,
                fileNameGDS=fileReferenceGDS, listProfiles=listProfiles,
                studyDF=studyDF, pathProfileGDS=pathProfileGDS)

    ## Open the 1KG GDS file (demo version)
    gds1KG <- snpgdsOpen(fileReferenceGDS)
    ## Open the 1KG GDS file and 1KG SNV Annotation file
    gdsAnnot1KG <- openfn.gds(fileReferenceAnnotGDS)

    listProfileRef <- dataRefSyn$sample.id
    studyDF.syn <- data.frame(study.id=paste0(studyDF$study.id, ".Synthetic"),
            study.desc=paste0(studyDF$study.id, " synthetic data"),
            study.platform=studyDF$study.platform, stringsAsFactors=FALSE)

    for(i in seq_len(length(listProfiles))) {
        pruningSample(gds=gds1KG, currentProfile=listProfiles[i],
                    studyID=studyDF$study.id, pathProfileGDS=pathProfileGDS)
        file.GDSProfile <- file.path(pathProfileGDS,
                                        paste0(listProfiles[i], ".gds"))
        add1KG2SampleGDS(gds=gds1KG, fileProfileGDS=file.GDSProfile,
                            currentProfile=listProfiles[i],
                            studyID=studyDF$study.id)
        addStudy1Kg(gds1KG, file.GDSProfile)

        gdsProfile <- openfn.gds(file.GDSProfile, readonly=FALSE)

        estimateAllelicFraction(gds=gds1KG, gdsSample=gdsProfile,
                                    currentProfile=listProfiles[i],
                                    studyID=studyDF$study.id, chrInfo=chrInfo)
        closefn.gds(gdsProfile)

        prepSynthetic(fileProfileGDS=file.GDSProfile,
                        listSampleRef=listProfileRef,
                        data.id.profile=listProfiles[i],
                        studyDF=studyDF.syn, prefId="1")

        resG <- syntheticGeno(gds=gds1KG, gdsRefAnnot=gdsAnnot1KG,
                                fileProfileGDS=file.GDSProfile,
                                data.id.profile=listProfiles[i],
                                listSampleRef=listProfileRef, prefId="1")

        if(! file.exists(pathOut)) {
            dir.create(pathOut)
        }
        spRef <- getRef1KGPop(gds1KG, "superPop")
        sampleRM <- splitSelectByPop(dataRefSyn)

        pathOutProfile <- file.path(pathOut, listProfiles[i])
        if(! file.exists(pathOutProfile)) {
            dir.create(pathOutProfile)
        }

        ## Open the Profile GDS file
        gdsProfile <- snpgdsOpen(file.GDSProfile)

        ## This variable will contain the results from the PCA analyses
        KNN.list <- list()
        ## For each row of the sampleRM matrix
        for(j in seq_len(nrow(sampleRM))) {
            ## Run a PCA analysis using 1 synthetic profile from each
            ##  sub-continental ancestry
            ## The synthetic profiles are projected on the 1KG PCA space
            ##  (the reference samples used to generate the synthetic profiles
            ##  are removed from this PCA)
            ## The K-nearest neighbor analysis is done using
            ##  a range of K and D values
            KNN.synt <- computePoolSyntheticAncestryGr(gdsSample=gdsProfile,
                sampleRM=sampleRM[j,], studyIDSyn=studyDF.syn$study.id,
                np=1L, spRef=spRef, eigen.cnt=15L)

            ## Results are saved
            saveRDS(KNN.synt$matKNN, file.path(pathOutProfile,
                        paste0("KNN.synt.", listProfiles[i], ".", j, ".rds")))
        }

        ## Directory where the KNN results have been saved
        PATHKNN <- file.path(pathOut, listProfiles[i])
        listFilesName <- dir(file.path(PATHKNN), ".rds")
        ## List of the KNN result files from PCA on synthetic data
        listFiles <- file.path(file.path(PATHKNN) , listFilesName)

        resCall <- computeAncestryFromSyntheticFile(gds=gds1KG,
                        gdsSample=gdsProfile, listFiles=listFiles,
                        sample.ana.id=listProfiles[i], spRef=spRef,
                        studyIDSyn=studyDF.syn$study.id, np=1L)

        saveRDS(resCall, file.path(pathOut,
                            paste0(listProfiles[i], ".infoCall", ".rds")))

        write.csv(resCall$Ancestry,
            file.path(pathOut, paste0(listProfiles[i], ".Ancestry",".csv")),
            quote=FALSE, row.names=FALSE)

        ## Close Profile GDS file (important)
        closefn.gds(gdsProfile)
    }
    closefn.gds(gds1KG)
    closefn.gds(gdsAnnot1KG)
    # Successful
    return(0L)
}



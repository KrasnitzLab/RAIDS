#' @title Create the GDS Sample file(s) for one or multiple specific samples
#' using the information from a RDS Sample description file and the 1KG
#' GDS file
#'
#' @description The function uses the information for the 1KG GDS file and the
#' RDS Sample Description file to create the GDS Sample file. One GDS Sample
#' file is created per sample. One GDS Sample file will be created for each
#' entry present in the \code{listSamples} parameter.
#'
#' @param PATHGENO a \code{character} string representing the path to the
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
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to the sample identifiers that will have a GDS Sample file created. The
#' sample identifiers must be present in the "Name.ID" column of the RDS file
#' passed to the \code{fileNamePED} parameter.
#' If \code{NULL}, all samples in the \code{fileNamePED} are selected.
#' Default: \code{NULL}.
#'
#' @param PATHSAMPLEGDS a \code{character} string representing the path to
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
#' ## Path to the demo 1KG GDS file is located in this package
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFile <- file.path(data.dir, "ex1_good_small_1KG_GDS.gds")
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
#' ## The Sample GDS file is created in the PATHSAMPLEGDS directory
#' result <- createStudy2GDS1KG(PATHGENO=data.dir,
#'             pedStudy=samplePED, fileNameGDS=gdsFile,
#'             studyDF=studyDF, listSamples=c("ex1"),
#'             PATHSAMPLEGDS=data.dir, verbose=FALSE)
#'
#' ## The function returns OL when successful
#' result
#'
#' ## The Sample GDS file 'ex1.gds' has been created in the
#' ## specified directory
#' list.files(data.dir)
#'
#' ## Unlink Sample GDS file (created for demo purpose)
#' unlink(file.path(data.dir, "ex1.gds"))
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
createStudy2GDS1KG <- function(PATHGENO=file.path("data", "sampleGeno"),
                                fileNamePED=NULL, pedStudy=NULL, fileNameGDS,
                                batch=1, studyDF, listSamples=NULL,
                                PATHSAMPLEGDS=NULL, verbose=TRUE) {

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
                batch=batch, studyDF=studyDF, listSamples=listSamples,
                PATHSAMPLEGDS=PATHSAMPLEGDS, verbose=verbose)

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

    generateGDS1KGgenotypeFromSNPPileup(PATHGENO=PATHGENO,
        listSamples=listSamples, listPos=listPos, offset=-1, minCov=10,
        minProb=0.999, seqError=0.001, pedStudy=pedStudy, batch=batch,
        studyDF=studyDF, PATHGDSSAMPLE=PATHSAMPLEGDS, verbose=verbose)

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
#' @param PATHGENO a \code{character} string representing the path to the
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
#' @param PATHSAMPLEGDS a \code{character} string representing the path to
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
appendStudy2GDS1KG <- function(PATHGENO=file.path("data", "sampleGeno"),
                                fileNamePED, fileNameGDS, batch=1,
                                studyDF, listSamples=NULL,
                                PATHSAMPLEGDS=NULL, verbose=TRUE) {

    ## Validate inputs
    validateAppendStudy2GDS1KG(PATHGENO=PATHGENO, fileNamePED=fileNamePED,
        fileNameGDS=fileNameGDS, batch=batch, studyDF=studyDF,
        listSamples=listSamples, PATHSAMPLEGDS=PATHSAMPLEGDS, verbose=verbose)

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

    generateGDS1KGgenotypeFromSNPPileup(PATHGENO=PATHGENO,
        listSamples=listSamples, listPos=listPos, offset=-1,
        minCov=10, minProb=0.999, seqError=0.001, pedStudy=pedStudy,
        batch=batch, studyDF=studyDF, PATHGDSSAMPLE=PATHSAMPLEGDS)

    if (verbose) {
        message("Genotype DONE ", Sys.time())
    }

    ## Close 1KG GDS file
    closefn.gds(gds)

    ## Return successful code
    return(0L)
}


#' @title Compute the list of pruned SNVs for a specific sample using the
#' information from the 1KG GDS file and a linkage disequilibrium analysis
#'
#' @description This function computes the list of pruned SNVs for a
#' specific sample. When
#' a group of SNVs are in linkage disequilibrium, only one SNV from that group
#' is retained. The linkage disequilibrium is calculated with the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The initial list of
#' SNVs that are passed to the \code{\link[SNPRelate]{snpgdsLDpruning}}()
#' function can be specified by the user.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1 KG GDS file.
#'
#' @param method a \code{character} string that represents the method that will
#' be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime" and "composite". Default: \code{"corr"}.
#'
#' @param sampleCurrent  a \code{character} string
#' corresponding to the sample identifier used in LD pruning done by the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. A GDS Sample file
#' corresponding to the sample identifier must exist and be located in the
#' \code{PATHSAMPLEGDS} directory.
#'
#' @param study.id a \code{character} string corresponding to the study
#' identifier used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#' The study identifier must be present in the GDS Sample file.
#'
#' @param listSNP a \code{vector} of SNVs identifiers specifying selected to
#' be passed the the pruning function;
#' if \code{NULL}, all SNVs are used in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}} function. Default: \code{NULL}.
#'
#' @param slide.max.bp.v a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done in the \code{\link[SNPRelate]{snpgdsLDpruning}}
#' function.
#' Default: \code{500000L}.
#'
#' @param ld.threshold.v a single \code{numeric} value that represents the LD
#' threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#' Default: \code{sqrt(0.1)}.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param verbose.v a \code{logicial} indicating if information is shown
#' during the process in the \code{\link[SNPRelate]{snpgdsLDpruning}}
#' function.  Default: \code{FALSE}.
#'
#' @param chr a \code{character} string representing the chromosome where the
#' selected SNVs should belong. Only one chromosome can be handled. If
#' \code{NULL}, the chromosome is not used as a filtering criterion.
#' Default: \code{NULL}.
#'
#' @param minAF.SuperPop a single positive \code{numeric} representing the
#' minimum allelic frequency used to select the SNVs. If \code{NULL}, the
#' allelic frequency is not used as a filtering criterion. Default: \code{NULL}.
#'
#' @param keepGDSpruned a \code{logicial} indicating if the information about
#' the pruned SNVs should be added to the GDS Sample file.
#' Default: \code{TRUE}.
#'
#' @param PATHSAMPLEGDS a \code{character} string representing the directory
#' where the GDS Sample file will be created. The directory must exist.
#'
#' @param keepFile a \code{logical} indicating if RDS files containing the
#' information about the pruned SNVs must be
#' created. Default: \code{FALSE}.
#'
#' @param PATHPRUNED a \code{character} string representing an existing
#' directory. The directory must exist. Default: \code{"."}.
#'
#' @param outPref a \code{character} string that represents the prefix of the
#' RDS files that will be generated. The RDS files are only generated when
#' the parameter \code{keepFile}=\code{TRUE}. Default: \code{"pruned"}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Open the 1KG GDS file (demo version)
#' gdsFile <- file.path(data.dir, "ex1_good_small_1KG_GDS.gds")
#' gds_1KG <- snpgdsOpen(gdsFile)

#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyDF <- data.frame(study.id = "MYDATA",
#'                         study.desc = "Description",
#'                         study.platform = "PLATFORM",
#'                         stringsAsFactors = FALSE)
#'
#'
#' ## TODO
#' fileNamePED <- "TODO"
#'
#'
#' ## Close the 1KG GDS file (it is important to always close the GDS files)
#' closefn.gds(gds_1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
pruningSample <- function(gds, method=c("corr", "r", "dprime", "composite"),
                            sampleCurrent,
                            study.id,
                            listSNP=NULL,
                            slide.max.bp.v=500000L,
                            ld.threshold.v=sqrt(0.1),
                            np=1L,
                            verbose.v=FALSE,
                            chr=NULL,
                            minAF.SuperPop=NULL,
                            keepGDSpruned=TRUE,
                            PATHSAMPLEGDS=NULL,
                            keepFile=FALSE,
                            PATHPRUNED=".",
                            outPref="pruned") {

    ## Validate input parameters
    validatePruningSample(gds=gds, method=method, sampleCurrent=sampleCurrent,
            study.id=study.id, listSNP=listSNP, slide.max.bp.v=slide.max.bp.v,
            ld.threshold.v=ld.threshold.v, np=np, verbose.v=verbose.v, chr=chr,
            minAF.SuperPop=minAF.SuperPop, keepGDSpruned=keepGDSpruned,
            PATHSAMPLEGDS=PATHSAMPLEGDS, keepFile=keepFile,
            PATHPRUNED=PATHPRUNED, outPref=outPref)

    ## Matches a character method against a table of candidate values
    method <- match.arg(method, several.ok=FALSE)

    ## GDS sample file name
    fileGDSSample <- file.path(PATHSAMPLEGDS, paste0(sampleCurrent, ".gds"))

    ## The GDS Sample file must exists
    if (!(file.exists(fileGDSSample))) {
        stop("The GDS Sample file \'", fileGDSSample, " does not exist.")
    }

    filePruned <- file.path(PATHPRUNED, paste0(outPref, ".rds"))
    fileObj <- file.path(PATHPRUNED, paste0(outPref, ".Obj.rds"))

    snp.id <- read.gdsn(node=index.gdsn(gds, "snp.id"))

    sample.id <- read.gdsn(node=index.gdsn(gds, "sample.id"))

    ## Open the GDS Sample file
    gdsSample <- openfn.gds(filename=fileGDSSample)

    ## Extract all study information from the GDS Sample file
    study.annot <- read.gdsn(node=index.gdsn(gdsSample, "study.annot"))

    ## Select study information associated to the current sample
    posSample <- which(study.annot$data.id == sampleCurrent &
                            study.annot$study.id == study.id)

    if(length(posSample) != 1) {
        stop("In pruningSample the sample ",
                sampleCurrent, " doesn't exists\n")
    }

    ## Get the SNV genotype information for the current sample
    g <- read.gdsn(index.gdsn(gdsSample, "geno.ref"),
                    start=c(1, posSample), count=c(-1,1))

    ## Close the GDS Sample file
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
    if(!is.null(minAF.SuperPop)) {
        listTMP <- NULL
        for(sp in c("EAS", "EUR", "AFR", "AMR", "SAS")) {
            snpAF <- read.gdsn(index.gdsn(gds, paste0("snp.", sp, "_AF")))
            listTMP <- union(listTMP,
                which(snpAF >= minAF.SuperPop & snpAF <= 1 - minAF.SuperPop))
        }
        listKeepPos <- intersect(listTMP, listKeepPos)
    }

    if(length(listKeepPos) == 0) {
        stop("In pruningSample, the sample ", sampleCurrent,
                " doesn't have SNPs after filters\n")
    }
    listKeep <- snp.id[listKeepPos]

    sample.ref <- read.gdsn(index.gdsn(gds, "sample.ref"))
    listSamples <- sample.id[which(sample.ref == 1)]

    ## Use a LD analysis to generate a subset of SNPs
    snpset <- runLDPruning(gds=gds, method=method, listSamples=listSamples,
                listKeep=listKeep, slide.max.bp.v=slide.max.bp.v,
                ld.threshold.v=ld.threshold.v, np=np, verbose.v=verbose.v)

    pruned <- unlist(snpset, use.names=FALSE)

    ## When TRUE, generate 2 RDS file with the pruned SNVs information
    if(keepFile) {
        saveRDS(pruned, filePruned)
        saveRDS(snpset, fileObj)
    }

    ## When TRUE, add the pruned SNvs information to the GDS Sample file
    if(keepGDSpruned) {
        gdsSample <- openfn.gds(filename=fileGDSSample, readonly=FALSE)
        addGDSStudyPruning(gds=gdsSample, pruned=pruned)
        closefn.gds(gdsfile=gdsSample)
    }

    return(0L)
}


#' @title Add the information about the pruned SNVs into the GDS Sample file
#'
#' @description The function extracts the information about the pruned SNVs
#' from the 1KG GDS file and adds entries related to the pruned SNVs in
#' the GDS Sample file.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param gdsSampleFile a \code{character} string representing the path and
#' file name of the GDS Sample file. The GDS Sample file must exist.
#'
#' @param sampleCurrent a \code{character} string corresponding to the sample
#' identifier associated to the current list of pruned SNVs.
#'
#' @param study.id a \code{character} string corresponding to the study
#' identifier associated to the current list of pruned SNVs.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn objdesp.gdsn
#' @encoding UTF-8
#' @export
add1KG2SampleGDS <- function(gds, gdsSampleFile, sampleCurrent,
                                study.id) {

    ## Validate inputs
    validateAdd1KG2SampleGDS(gds=gds, gdsSampleFile=gdsSampleFile,
                sampleCurrent=sampleCurrent, study.id=study.id)

    ## Open GDS Sample file
    gdsSample <- openfn.gds(gdsSampleFile, readonly=FALSE)

    ## Extract needed information from 1KG GDS file
    snp.id <- read.gdsn(index.gdsn(gds,"snp.id"))

    ## Extract list of pruned SNVs from the GDS Sample file
    pruned <- read.gdsn(index.gdsn(gdsSample, "pruned.study"))

    listSNP <- which(snp.id %in% pruned)
    listRef <- which(read.gdsn(index.gdsn(gds, "sample.ref")) == 1)
    sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))

    snp.chromosome <- read.gdsn(index.gdsn(gds,"snp.chromosome"))[listSNP]
    snp.position <-  read.gdsn(index.gdsn(gds,"snp.position"))[listSNP]

    add.gdsn(gdsSample, "sample.id", c(sample.id[listRef], sampleCurrent))

    add.gdsn(gdsSample, "snp.id", snp.id[listSNP])
    add.gdsn(gdsSample, "snp.chromosome", snp.chromosome)
    add.gdsn(gdsSample, "snp.position", snp.position)
    # snp.index is the index of the snp pruned in snp.id fro 1KG gds
    add.gdsn(gdsSample, "snp.index", listSNP)


    var.geno <- NULL

    j <- 1
    for(i in listRef){
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

    posCur <- which(study.annot$data.id == sampleCurrent &
                            study.annot$study.id == study.id)

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
#' @param PATHSAMPLEGDS a \code{character} string representing the path to
#' the directory that contains the GDS Sample files. The directory must
#' exist.
#'
#' @param PATHGENO TODO
#'
#' @param fileLSNP TODO
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{TRUE}.
#'
#' @return The integer \code{0} when successful.
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
addPhase1KG2SampleGDSFromFile <- function(gds, PATHSAMPLEGDS, PATHGENO,
                                            fileLSNP, verbose=FALSE) {

    listGDSSample <- dir(PATHSAMPLEGDS, pattern = ".+.gds")


    indexAll <- NULL
    for(gdsSampleFile in listGDSSample) {
        gdsSample <- openfn.gds(filename=file.path(PATHSAMPLEGDS,
                                                        gdsSampleFile))

        snp.index <- read.gdsn(node=index.gdsn(node=gdsSample, "snp.index"))

        indexAll <- union(indexAll, snp.index)
        closefn.gds(gdsfile=gdsSample)
    }

    gdsSample <- createfn.gds(file.path(PATHSAMPLEGDS, "phase1KG.gds"))
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
        file1KG <- file.path(PATHGENO, paste0(sample1KG,".csv.bz2"))
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
#' @param PATHSAMPLEGDS the path of an object of class \code{gds} related to
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
addPhase1KG2SampleGDSFromGDS <- function(gds, gdsPhase, PATHSAMPLEGDS,
                                            verbose=FALSE) {

    listGDSSample <- dir(PATHSAMPLEGDS, pattern = ".+.gds")


    indexAll <- NULL
    for(gdsSampleFile in listGDSSample){
        gdsSample <- openfn.gds(file.path(PATHSAMPLEGDS, gdsSampleFile))

        snp.index <- read.gdsn(index.gdsn(gdsSample,"snp.index"))

        indexAll <- union(indexAll, snp.index)
        closefn.gds(gdsSample)
    }

    gdsSample <- createfn.gds(file.path(PATHSAMPLEGDS, "phase1KG.gds"))
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

        #file1KG <- file.path(PATHGENO, paste0(sample1KG,".csv.bz2"))
        #matSample <- read.csv2( file1KG,
        #                        row.names = NULL)
        #matSample <- matSample[listSNP[indexAll],, drop=FALSE]
        #matSample <- matrix(as.numeric(unlist(strsplit( matSample[,1],
        #                        "\\|"))),nrow=2)[1,]
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
#' identifier of the samples in the reference (unrelated).
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
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
computePrunedPCARef <- function(gds, listRef, np=1L) {

    ## Validate that np is a single positive integer
    if(! (isSingleNumber(np) && np > 0)) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    listPCA <- list()

    listPruned <- read.gdsn(index.gdsn(gds, "pruned.study"))


    ## Calculate the eigenvectors using the specified SNP loadings for
    ## the reference samples
    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsobj=gds,
                                            sample.id=listRef,
                                            snp.id=listPruned,
                                            num.thread=np,
                                            verbose=TRUE)

    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(pcaobj=listPCA[["pca.unrel"]],
                                                    gdsobj=gds,
                                                    num.thread=np,
                                                    verbose=TRUE)
    return(listPCA)
}



#' @title Project patients onto existing principal component axes (PCA)
#'
#' @description This function calculates the patient eigenvectors using
#' the specified SNP loadings.
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param listPCA  a \code{list} containing two objects
#' pca.unrel -> \code{snpgdsPCAClass}
#' and a snp.load -> \code{snpgdsPCASNPLoading}
#'
#' @param sample.current a \code{character} string representing the
#' identifiant of the sample to be projected in the PCA.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
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
projectSample2PCA <- function(gds, listPCA, sample.current, np=1L) {


    ## Validate that sample.current is a character string
    if(! is.character(sample.current)) {
        stop("The \'sample.current\' parameter must be a character string.")
    }

    ## Validate that np is a single positive integer
    if(! (isSingleNumber(np) && np > 0)) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## Calculate the sample eigenvectors using the specified SNP loadings
    samplePCA <- snpgdsPCASampLoading(listPCA[["snp.load"]],
                                gdsobj=gds, sample.id=sample.current,
                                num.thread=1, verbose=TRUE)

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
#' @param sampleCurrent a \code{character} string corresponding to
#' the sample identifier as used in \code{\link{pruningSample}} function.
#'
#' @param study.id a \code{character} string corresponding to the name of
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
estimateAllelicFraction <- function(gds, gdsSample, sampleCurrent, study.id,
    chrInfo, studyType=c("DNA", "RNA"), minCov=10L, minProb=0.999,
    eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=9, cutOffAR=3,
    gdsRefAnnot=NULL, block.id=NULL) {

    ## Validate input parameters
    validateEstimateAllelicFraction(gds=gds, gdsSample=gdsSample,
        sampleCurrent=sampleCurrent, study.id=study.id, chrInfo=chrInfo,
        studyType=studyType, minCov=minCov, minProb=minProb, eProb=eProb,
        cutOffLOH=cutOffLOH, cutOffHomoScore=cutOffHomoScore, wAR=wAR,
        cutOffAR=cutOffAR, gdsRefAnnot=gdsRefAnnot, block.id=block.id)

    ## Set study type
    studyType <- match.arg(studyType)

    snp.pos <- NULL

    ## The type of study affects the allelic fraction estimation
    if(studyType == "DNA") {
        snp.pos <- computeAllelicFractionDNA(gds=gds, gdsSample=gdsSample,
                        sampleCurrent=sampleCurrent, study.id=study.id,
                        chrInfo=chrInfo, minCov=minCov, minProb=minProb,
                        eProb=eProb, cutOffLOH=cutOffLOH,
                        cutOffHomoScore=cutOffHomoScore, wAR=wAR)

    } else if(studyType == "RNA") {
        snp.pos <- computeAllelicFractionRNA(gds=gds, gdsSample=gdsSample,
                        gdsRefAnnot=gdsRefAnnot, sampleCurrent=sampleCurrent,
                        study.id=study.id, block.id=block.id, chrInfo=chrInfo,
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
#' @param gdsSampleFile a \code{character} string representing the path and
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
addStudy1Kg <- function(gds, gdsSampleFile) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds, "gds")

    ## The gdsSampleFile must be a character string and the file must exists
    if(!(is.character(gdsSampleFile) && (file.exists(gdsSampleFile)))) {
        stop("The \'gdsSampleFile\' must be a character string representing ",
                "the GDS Sample file. The file must exist.")
    }

    ## Open GDS Sample file
    gdsSample <- openfn.gds(filename=gdsSampleFile, readonly=FALSE)

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
#' @param study.id.syn a \code{character} string corresponding to the study
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
                                        sampleRef, study.id.syn) {

    if(length(sampleRef) < 1) {
        stop("Number of sample in study.annot not equal to 1\n")
    }

    study.annot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))
    study.annot <- study.annot[which(study.annot$study.id == study.id.syn &
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
#' @param study.id.ref a single \code{character} string representing the
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
computePCARefSample <- function(gdsSample, name.id, study.id.ref="Ref.1KG",
                            np=1L, algorithm=c("exact","randomized"),
                            eigen.cnt=32L, missing.rate=NaN) {

    ## Validate parameters
    validateComputePCARefSample(gdsSample=gdsSample, name.id=name.id,
        study.id.ref=study.id.ref, np=np, algorithm=algorithm,
        eigen.cnt=eigen.cnt, missing.rate=missing.rate)

    ## Set algorithm
    algorithm <- match.arg(algorithm)

    sample.id <- read.gdsn(index.gdsn(gdsSample, "sample.id"))

    sample.pos <- which(sample.id == name.id)

    study.annot.all <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    sample.Unrel <- study.annot.all[which(study.annot.all$study.id ==
                                                study.id.ref), "data.id"]

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
#' @param study.id.syn a the study.id of the synthetic data
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
                                    listCatPop, study.id.syn,
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
                                                study.id.syn &
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
#' @param study.id.syn a \code{character} string corresponding to the study
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
                            study.id.syn, np=1L,
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
        sampleRM=sampleRM, spRef=spRef, study.id.syn=study.id.syn,
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
                        sampleRef=sampleRM, study.id.syn=study.id.syn)

    ## Calculate the k-nearest neighbor analyses on a subset of the
    ## synthetic dataset
    KNN.synt <- computeKNNRefSynthetic(gdsSample=gdsSample,
                        listEigenvector=resPCA, listCatPop=listCatPop,
                        study.id.syn=study.id.syn, spRef=spRef,
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
#' @param study.id.syn a \code{character} string corresponding to the study
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
                                            study.id.syn,
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

    ## The study.id must be a character string
    if (!(is.character(study.id.syn) && length(study.id.syn) == 1)) {
        stop("The \'study.id.syn\' parameter must be a character string.")
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
                            spRef=spRef, study.id.syn=study.id.syn,
                            np=np, listCatPop=listCatPop,
                            fieldPopIn1KG=fieldPopIn1KG,
                            fieldPopInfAnc=fieldPopInfAnc, kList=kList,
                            pcaList=pcaList, algorithm=algorithm,
                            eigen.cnt=eigen.cnt, missing.rate=missing.rate)
    }

    KNN.sample.syn <- do.call(rbind, KNN.list)

    pedSyn <- prepPedSynthetic1KG(gds=gds, gdsSample=gdsSample,
                            study.id=study.id.syn, popName=fieldPopIn1KG)

    listParaSample <- selParaPCAUpQuartile(KNN.sample.syn, pedSyn,
                                            fieldPopIn1KG, fieldPopInfAnc,
                                            listCatPop)

    listPCASample <- computePCARefSample(gdsSample=gdsSample,
                            name.id=sample.ana.id, study.id.ref="Ref.1KG",
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
#' @param study.id.syn a \code{character} string corresponding to the study
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
                            study.id.syn,
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
        study.id.syn=study.id.syn, np=np, listCatPop=listCatPop,
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
        study.id=study.id.syn, popName=fieldPopIn1KG)

    listParaSample <- selParaPCAUpQuartile(matKNN.All=KNN.sample.syn,
        pedCall=pedSyn, refCall=fieldPopIn1KG, predCall=fieldPopInfAnc,
        listCall=listCatPop)

    listPCASample <- computePCARefSample(gdsSample=gdsSample,
        name.id=sample.ana.id, study.id.ref="Ref.1KG", np=np,
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





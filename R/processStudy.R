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
    if(!(is.null(fileNamePED)) && is.null(pedStudy)) {
        ## The fileNamePED must be a character string and the file must exists
        if (!(is.character(fileNamePED) && (file.exists(fileNamePED)))) {
            stop("The \'fileNamePED\' must be a character string ",
                    "representing the RDS Sample information file. ",
                    "The file must exist.")
        }
        ## Open the RDS Sample information file
        pedStudy <- readRDS(file=fileNamePED)
    } else if(!(is.null(fileNamePED) || is.null(pedStudy))) {
        stop("Both \'fileNamePED\' and \'pedStudy\' parameters cannot be ",
                "defined at the same time.")
    } else if(is.null(fileNamePED) && is.null(pedStudy)) {
        stop("One of the parameter \'fineNamePED\' of \'pedStudy\' must ",
                "be defined.")
    }

    ## The PED study must have the mandatory columns
    if(!(all(c("Name.ID", "Case.ID", "Sample.Type", "Diagnosis", "Source")
             %in% colnames(pedStudy)))) {
        stop(paste0("The PED study data frame is incomplete. ",
                        "One or more mandatory columns are missing."))
    }

    ## The fileNameGDS must be a character string and the file must exists
    if(!(is.character(fileNameGDS) && (file.exists(fileNameGDS)))) {
        stop("The \'fileNameGDS\' must be a character string representing ",
             "the GDS 1KG file. The file must exist.")
    }

    ## The batch must be a single numeric
    if(!(isSingleNumber(batch))) {
        stop("The \'batch\' must be a single integer.")
    }

    ## The listSamples must be a vector of character string
    if(!(is.character(listSamples) || is.null(listSamples))) {
        stop("The \'listSamples\' must be a vector ",
             "of character strings (1 entry or more) or NULL.")
    }

    ## The verbose parameter must be a logical
    if(!(is.logical(verbose))) {
        stop("The \'verbose\' parameter must be a logical (TRUE or FALSE).")
    }

    ## Read the 1KG GDS file
    gds <- snpgdsOpen(filename=fileNameGDS)

    ## Extract the chromosome and position information for all SNPs in 1KG GDS
    ## Create a data.frame containing the information
    snpCHR <- index.gdsn(node=gds, "snp.chromosome")
    snpPOS <- index.gdsn(node=gds, "snp.position")

    listPos <- data.frame(snp.chromosome=read.gdsn(snpCHR),
                            snp.position=read.gdsn(snpPOS))

    if(verbose) {
        message("Start ", Sys.time())
        message("Sample info DONE ", Sys.time())
    }

    generateGDS1KGgenotypeFromSNPPileup(PATHGENO=PATHGENO,
        listSamples=listSamples, listPos=listPos, offset=-1,
        minCov=10, minProb=0.999, seqError=0.001, pedStudy=pedStudy,
        batch=batch, studyDF=studyDF, PATHGDSSAMPLE=PATHSAMPLEGDS)

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
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
appendStudy2GDS1KG <- function(PATHGENO=file.path("data", "sampleGeno"),
                                fileNamePED, fileNameGDS, batch=1,
                                studyDF, listSamples=NULL,
                                PATHSAMPLEGDS=NULL, verbose=TRUE) {

    ## The fileNamePED must be a character string and the file must exists
    if (!(is.character(fileNamePED) && (file.exists(fileNamePED)))) {
        stop("The \'fileNamePED\' must be a character string representing ",
             "the RDS Sample information file. The file must exist.")
    }

    ## The fileNameGDS must be a character string and the file must exists
    if (!(is.character(fileNameGDS) && (file.exists(fileNameGDS)))) {
        stop("The \'fileNameGDS\' must be a character string representing ",
                "the GDS 1KG file. The file must exist.")
    }

    ## The batch must be a single numeric
    if (!(isSingleNumber(batch))) {
        stop("The \'batch\' must be a single integer.")
    }

    ## The listSamples must be a vector of character string
    if (!(is.character(listSamples) || is.null(listSamples))) {
        stop("The \'listSamples\' must be a vector ",
                "of character strings (1 entry or more) or NULL.")
    }

    ## The verbose parameter must be a logical
    if (!(is.logical(verbose))) {
        stop("The \'verbose\' parameter must be a logical (TRUE or FALSE).")
    }

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
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
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

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'.")
    }

    ## The parameter sampleCurrent must be a character string
    if(!(is.character(sampleCurrent))) {
        stop("The \'sampleCurrent\' parameter must be a character string.")
    }

    ## The parameter method must be a character string
    if(!(is.character(method))) {
        stop("The \'method\' parameter must be a character string.")
    }

    ## Matches a character method against a table of candidate values
    method <- match.arg(method, several.ok=FALSE)

    ## The parameter ld.threshold.v must be a single positive integer
    if(!(isSingleNumber(ld.threshold.v) && (ld.threshold.v >= 0.0))) {
        stop("The \'ld.threshold.v\' parameter must be a single positive ",
                "numeric value.")
    }

    ## The parameter slide.max.bp.v must be a single positive integer
    if(!(isSingleNumber(slide.max.bp.v) && (slide.max.bp.v >= 0.0))) {
        stop("The \'slide.max.bp.v\' parameter must be a single positive ",
                "numeric value.")
    }

    ## The parameter np must be a single positive integer
    if(!(isSingleNumber(np) && (np >= 0.0))) {
        stop("The \'np\' parameter must be a single positive numeric value.")
    }

    ## The parameter keepGDSpruned must be a logical
    if(!is.logical(keepGDSpruned)) {
        stop("The \'keepGDSpruned\' parameter must be a logical ",
                "(TRUE or FALSE).")
    }

    ## The parameter PATHSAMPLEGDS must be a character string representing an
    ## existing path
    if(!(is.character(PATHSAMPLEGDS) && dir.exists(PATHSAMPLEGDS))) {
        stop("The \'PATHSAMPLEGDS\' parameter must be a character string ",
             "representing an existing directory.")
    }

    ## The parameter keepFile must be a logical
    if(!is.logical(keepFile)) {
        stop("The \'keepFile\' parameter must be a logical (TRUE or FALSE).")
    }

    ## The parameter PATHPRUNED must be a character string representing an
    ## existing path
    if(!(is.character(PATHPRUNED) && dir.exists(PATHPRUNED))) {
        stop("The \'PATHPRUNED\' parameter must be a character string ",
                "representing an existing directory.")
    }

    fileGDSSample <- file.path(PATHSAMPLEGDS, paste0(sampleCurrent, ".gds"))

    ## The GDS Sample file must exists
    if(!(file.exists(fileGDSSample))) {
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
            snpAF <- read.gdsn(index.gdsn(gds, paste0("snp.", sp, "_AF") ))
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
    snpset <- runLDPruning(gds=gds, method=method,
                            listSamples=listSamples,
                            listKeep=listKeep,
                            slide.max.bp.v=slide.max.bp.v,
                            ld.threshold.v=ld.threshold.v, np=np,
                            verbose.v=verbose.v)

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

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'.")
    }

    ## The gdsSampleFile must be a character string and the file must exists
    if(!(is.character(gdsSampleFile) && (file.exists(gdsSampleFile)))) {
        stop("The \'gdsSampleFile\' must be a character string representing ",
                "the GDS Sample file. The file must exist.")
    }

    ## The sampleCurrent must be a character string
    if(!(is.character(sampleCurrent))) {
        stop("The \'sampleCurrent\' must be a character string.")
    }

    ## The study.id must be a character string
    if(!(is.character(study.id))) {
        stop("The \'study.id\' must be a character string.")
    }

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
#' @param gds an object of class \code{gds} opened
#'
#' @param PATHSAMPLEGDS the path of an object of class \code{gds} related to
#' the sample
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
        if(!("phase" %in% ls.gdsn(gdsSample))) {
            var.phase <- add.gdsn(gdsSample, "phase",
                                 valdim=c(length(indexAll),
                                          1),
                                 matSample, storage="bit2")
        }else {
            if(is.null(var.phase)) {
                var.phase <- index.gdsn(node=gdsSample, "phase")
            }
            append.gdsn(node=var.phase, val=matSample)
        }
    }

    closefn.gds(gdsfile=gdsSample)

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
        if(! ("phase" %in% ls.gdsn(gdsSample))){
            var.phase <- add.gdsn(gdsSample, "phase",
                                  valdim=c(length(indexAll),
                                           1),
                                  matSample, storage="bit2")

        }else{
            if(is.null(var.phase)){
                var.phase <- index.gdsn(gdsSample, "phase")
            }
            append.gdsn(var.phase, matSample)
        }
    }

    closefn.gds(gdsSample)

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
#' identifiant of the samples in the reference (unrelated).
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


#' @title Deprecated - Project patients onto existing principal component
#' axes (PCA)
#'
#' @description This function calculates the patient eigenvectors using
#' the specified SNP loadings. Deprecated
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param PATHSAMPLEGDS the path of an object of class \code{gds} related to
#' the sample
#'
#' @param listSamples a \code{vector} of string representing the samples for
#' which compute the PCA.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @return The integer \code{0} when successful.
#'
#' @details
#'
#' More information about the method used to calculate the patient eigenvectors
#' can be found at the Bioconductor SNPRelate website:
#' https://bioconductor.org/packages/SNPRelate/
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom SNPRelate snpgdsPCASampLoading
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computePCAForSamples <- function(gds, PATHSAMPLEGDS, listSamples, np=1L) {

    ## Validate that np is a single positive integer
    if(! (isSingleNumber(np) && np > 0)) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    #sample.ref <- read.gdsn(index.gdsn(gds, "sample.ref"))

    #listRef <- read.gdsn(index.gdsn(gds, "sample.id"))[which(sample.ref == 1)]

    for(i in seq_len(length(listSamples)) ){

        gdsSample <- openfn.gds(file.path(PATHSAMPLEGDS,
                                            paste0(listSamples[i], ".gds")))
        study.annot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

        if(length(which(study.annot$study.id == "Ref.1KG")) == 0) {
            stop("The study Ref.1KG is not define you must run the ",
                                    "function addStudy1Kg \n")
        }

        sample.Unrel.All <- study.annot$data.id[study.annot$study.id ==
                                                                    "Ref.1KG"]
        #sample.ref <- sample.Unrel.All$data.id
        listPCA <- computePrunedPCARef(gdsSample, sample.Unrel.All, np)

        listPCA[["samp.load"]] <- projectSample2PCA(gdsSample, listPCA,
                                                            listSamples[i], np)
        closefn.gds(gdsSample)

        saveRDS(listPCA, file.path(PATHSAMPLEGDS, paste0(listSamples[i],
                                                        ".pca.pruned.rds")))

    }

    return(0L)
}


#' @title Estimate the allelic fraction of the pruned SNVs for a specific
#' sample
#'
#' @description The function estimates the allelic fraction of the pruned
#' SNVs for a specific sample and add the information to the associated
#' GDS Sample file.
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
#' @param study.id a \code{character} string corresponding to the study
#' identifier.
#'
#' @param chrInfo a \code{vector}
#'  chrInfo[i] = length(Hsapiens[[paste0("chr", i)]])
#'         Hsapiens library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' @param studyType a \code{character} string representing the type of study.
#' It can be values such as "DNA", "whole-exome", "RNA", etc..
#' Default: \code{"DNA"}.
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
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file. RNA specific
#' Default: \code{NULL}.
#'
#' @param block.id a \code{character} string corresponding to the block
#' identifier in \code{gdsRefAnnot}.  RNA specific
#' Default: \code{NULL}
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
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
estimateAllelicFraction <- function(gds, gdsSample, sampleCurrent, study.id,
                                    chrInfo, studyType="DNA",
                                    minCov=10L, minProb=0.999, eProb=0.001,
                                    cutOffLOH=-5, cutOffHomoScore=-3,
                                    wAR=9, gdsRefAnnot=NULL, block.id = NULL) {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'.")
    }

    ## The gdsSample must be an object of class "gds.class"
    if (!inherits(gdsSample, "gds.class")) {
        stop("The \'gdsSample\' must be an object of class \'gds.class\'.")
    }

    ## The sampleCurrent must be a character string
    if (!is.character(sampleCurrent)) {
        stop("The \'sampleCurrent\' must be a character string.")
    }

    ## The study.id must be a character string
    if (!is.character(study.id)) {
        stop("The \'study.id\' must be a character string.")
    }

    ## The studyType must be a character string
    if (!is.character(studyType)) {
        stop("The \'studyType\' must be a character string.")
    }

    ## The minCov parameter must be a single positive integer
    if (!(isSingleNumber(minCov) && (minCov >= 0.0))) {
        stop("The \'minCov\' must be a single numeric positive value")
    }

    ## The minProb parameter must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(minProb) && (minProb >= 0.0) && (minProb <= 1.0))) {
        stop("The \'minProb\' must be a single numeric positive ",
             "value between 0 and 1.")
    }

    ## The eProb parameter must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(eProb) && (eProb >= 0.0) && (eProb <= 1.0))) {
        stop("The \'eProb\' must be a single numeric positive ",
             "value between 0 and 1.")
    }

    ## The wAR parameter must be a single positive numeric superior to 1
    if (!(isSingleNumber(wAR) && (wAR >= 1))) {
        stop("The \'wAR\' parameter must be a single numeric positive value.")
    }

    ## The cutOffLOH parameter must be a single numeric
    if (!(isSingleNumber(cutOffLOH))) {
        stop("The \'cutOffLOH\' parameter must be a single numeric value.")
    }

    ## The cutOffHomoScore parameter must be a single numeric
    if (!(isSingleNumber(cutOffHomoScore))) {
        stop("The \'cutOffHomoScore\' parameter must be a single ",
                        "numeric value.")
    }

    snp.pos <- NULL
    if(studyType == "DNA") {
        snp.pos <- computeAllelicFractionDNA(gds=gds, gdsSample=gdsSample,
                                            sampleCurrent=sampleCurrent,
                                            study.id, chrInfo,
                                            minCov=minCov, minProb=minProb,
                                            eProb=eProb,
                                            cutOffLOH=cutOffLOH,
                                            cutOffHomoScore=cutOffHomoScore,
                                            wAR=wAR)

        snp.pos$seg <- rep(0, nrow(snp.pos))
        k <- 1
        # Find segment with same lap
        for(chr in seq_len(22)) {
            snpChr <- snp.pos[snp.pos$snp.chr == chr, ]
            tmp <- c(0,
                     abs(snpChr[2:nrow(snpChr), "lap"] -
                             snpChr[seq_len(nrow(snpChr)- 1),  "lap"]) > 1e-3)
            snp.pos$seg[snp.pos$snp.chr == chr] <- cumsum(tmp) + k
            k <- max(snp.pos$seg[snp.pos$snp.chr == chr]) + 1
        }
    } else if(studyType == "RNA"){

    }

    ## Save information into the "lap" node in the GDS Sample file
    ## Suppose we keep only the pruned SNV
    addUpdateLap(gdsSample, snp.pos$lap[which(snp.pos$pruned == TRUE)])

    ## Save information into the "segment" node in the GDS Sample file
    ## Suppose we keep only the pruned SNV
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
#' @return The integer \code{0} when successful.
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
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'.")
    }

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
        ped1KG <- data.frame(Name.ID=sample.id,
                            Case.ID=sample.id,
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


#' @title Deprecated Function
#'
#' @description TODO Deprecated
#'
#' @param gdsSample an object of class \code{gds} opened related to
#' the sample
#'
#' @param pruned TODO
#'
#' @param sample.id TODO
#'
#' @param sample.ref TODO
#'
#' @param study.annot a  \code{data.frame} with one entry from study.annot in
#' the gds
#'
#' @param algorithm algorithm of the PCA "exact", "randomized"
#'
#' @param eigen.cnt number of eigenvectors in PCA
#'
#' @return A \code{list} TODO with the sample.id and eigenvectors.
#'
#' @examples
#'
#' # TODO
#' gds <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASampLoading snpgdsPCASampLoading
#' @encoding UTF-8
#' @keywords internal
computePCAsynthetic <- function(gdsSample, pruned, sample.id,
                                    sample.ref, study.annot,
                                    algorithm="exact", eigen.cnt=32L) {

    if(nrow(study.annot) != 1) {
        stop("Number of sample in study.annot not equal to 1\n")
    }

    sample.pos <- which(sample.id == study.annot$data.id[1])
    sample.Unrel <- sample.ref[which(sample.ref != study.annot$case.id[1])]

    g <- read.gdsn(index.gdsn(gdsSample, "genotype"),
                    start=c(1, sample.pos), count=c(-1, 1))

    listPCA <- list()

    listPCA[["pruned"]] <- pruned[which(g != 3)]
    rm(g)

    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsSample,
                                    sample.id=sample.Unrel,
                                    snp.id=listPCA[["pruned"]],
                                    num.thread=1,
                                    algorithm=algorithm,
                                    eigen.cnt=eigen.cnt,
                                    verbose=TRUE)

    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                        gdsobj=gdsSample,
                                        num.thread=1,
                                        verbose=TRUE)

    listPCA[["samp.load"]] <- snpgdsPCASampLoading(listPCA[["snp.load"]],
                                        gdsobj=gdsSample,
                                        sample.id=sample.id[sample.pos],
                                        num.thread=1, verbose=TRUE)

    listRes <- list(sample.id=sample.id[sample.pos],
                        eigenvector.ref=listPCA[["pca.unrel"]]$eigenvect,
                        eigenvector=listPCA[["samp.load"]]$eigenvect)
    return(listRes)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param gdsSample an object of class \code{gds} opened related to
#' the sample
#'
#' @param sample.ref TODO
#'
#' @param listRM a  \code{} list of sample from the Ref to remove
#' before the PCA
#'
#' @param np TODO
#'
#' @param algorithm algorithm of the PCA "exact", "randomized"
#' (para snpgdsPCA)
#'
#' @param eigen.cnt number of eigenvectors in PCA
#' (para snpgdsPCA)
#'
#' @param missing.rate number of eigenvectors in PCA
#' (para snpgdsPCA)
#'
#' @return A \code{list} TODO with the sample.id and eigenvectors.
#'
#' @examples
#'
#' # TODO
#' gds <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA
#' @encoding UTF-8
#' @export
computePCARefRMMulti <- function(gdsSample,
                                sample.ref, listRM, np=1L,
                                algorithm="exact", eigen.cnt=32L,
                                missing.rate=0.025) {

    if(length(listRM) < 1) {
        stop("Number of sample in study.annot not equal 0\n")
    }


    sample.Unrel <- sample.ref[which(!(sample.ref %in% listRM) )]


    listPCA <- list()

    listPCA[["pruned"]] <- read.gdsn(index.gdsn(gdsSample, "pruned.study"))

    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsSample,
                                        sample.id=sample.Unrel,
                                        snp.id=listPCA[["pruned"]],
                                        num.thread=np,
                                        missing.rate=missing.rate,
                                        algorithm=algorithm,
                                        eigen.cnt=eigen.cnt,
                                        verbose=TRUE)


    return(listPCA)
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
#' @param study.id.syn a the study.id of the synthetic data
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


    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                                 gdsobj=gdsSample,
                                                 num.thread=1,
                                                 verbose=TRUE)

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


#' @title TODO
#'
#' @description TODO
#'
#' @param gdsSample an object of class \link[gdsfmt]{gds.class},
#' a GDS Sample file.
#'
#' @param name.id a \code{character} string representing the sample
#' identifier.
#'
#' @param study.id.ref id of the reference in study.annot
#'
#' @param np a single positive \code{integer} representing the number of CPU
#' that will be used. Default: \code{1L}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigen.cnt TODO. Default: \code{32L}.
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
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
computePCARefSample <- function(gdsSample, name.id, study.id.ref="Ref.1KG",
                            np=1L, algorithm=c("exact","randomized"),
                            eigen.cnt=32L) {

    ## Validate that name.id is a string
    if(!(is.character(name.id))) {
        stop("The \'name.id\' parameter must be a character string.")
    }

    ## Validate that np is a positive number
    if(!(isSingleNumber(np) && np > 0)) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    ## Validate that algorithm is valid choice
    algorithm <- match.arg(algorithm)

    if(length(name.id) != 1) {
        stop("Number of sample in study.annot not equal to 1\n")
    }

    sample.id <- read.gdsn(index.gdsn(gdsSample, "sample.id"))

    sample.pos <- which(sample.id == name.id)

    study.annot.all <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    sample.Unrel <- study.annot.all[which(study.annot.all$study.id ==
                                                study.id.ref), "data.id"]


    listPCA <- list()

    listPCA[["pruned"]] <- read.gdsn(index.gdsn(gdsSample, "pruned.study"))

    ## Calculate the eigenvectors and eigenvalues for PCA
    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsSample,
                                            sample.id=sample.Unrel,
                                            snp.id=listPCA[["pruned"]],
                                            num.thread=np,
                                            algorithm=algorithm,
                                            eigen.cnt=eigen.cnt,
                                            verbose=TRUE)

    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                                    gdsobj=gdsSample,
                                                    num.thread=np,
                                                    verbose=TRUE)

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


#' @title Deprecated
#'
#' @description Deprecated
#'
#' @param listEigenvector TODO see return of computePCAsynthetic
#'
#' @param sample.ref TODO
#'
#' @param study.annot a  \code{data.frame} with one entry from study.annot in
#' the gds
#'
#' @param spRef TODO
#'
#' @param kList TODO array of the k possible values
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
#' @importFrom gdsfmt add.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASampLoading snpgdsPCASampLoading
#' @importFrom class knn
#' @encoding UTF-8
#' @keywords internal
computeKNNSuperPoprSynthetic <- function(listEigenvector, sample.ref,
                                         study.annot, spRef,
                                         kList = seq_len(15), pcaList = 2:15) {

    ## The number of rows in study.annot must be one.
    if(nrow(study.annot) != 1) {
        stop("Number of samples in study.annot not equal to 1\n")
    }

    if(is.null(kList)){
        kList <- seq_len(15)#c(seq_len(14), seq(15,100, by=5))
    }
    if(is.null(pcaList)){
        pcaList <- 2:15
    }

    resMat <- data.frame(sample.id=rep(listEigenvector$sample.id,
                                            length(pcaList) * length(kList)),
                            D=rep(0,length(pcaList) * length(kList)),
                            K=rep(0,length(pcaList) * length(kList)),
                            SuperPop=character(length(pcaList) * length(kList)),
                            stringsAsFactors=FALSE)

    listSuperPop <- c("EAS", "EUR", "AFR", "AMR", "SAS")

    #curPCA <- listPCA.Samples[[sample.id[sample.pos]]]
    eigenvect <- rbind(listEigenvector$eigenvector.ref,
                            listEigenvector$eigenvector)

    rownames(eigenvect) <- c(sample.ref[which(sample.ref !=
                                                    study.annot$case.id[1])],
                                listEigenvector$sample.id)

    totR <- 1
    for(pcaD in pcaList) {
        for(kV in  seq_len(length(kList))) {
            dCur <- paste0("d", pcaD)
            kCur <- paste0("k", kList[kV])
            resMat[totR,c("D", "K")] <- c(pcaD, kList[kV])

            pcaND <- eigenvect[ ,seq_len(pcaD)]
            y_pred <- knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                test=pcaND[rownames(eigenvect)[nrow(eigenvect)],, drop=FALSE],
                cl=factor(spRef[rownames(eigenvect)[-1*nrow(eigenvect)]],
                            levels=listSuperPop, labels=listSuperPop),
                k=kList[kV],
                prob=FALSE)

            resMat[totR, paste0("SuperPop")] <- listSuperPop[as.integer(y_pred)]

            totR <- totR + 1
        } # end k
    } # end pca Dim
    listKNN <- list(sample.id=listEigenvector$sample.id,
                        sample1Kg=study.annot$case.id[1],
                        sp=spRef[study.annot$case.id[1]],
                        matKNN=resMat)

    return(listKNN)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param gdsSample an object of class \code{gds} opened related to
#' the sample
#'
#' @param listEigenvector TODO see return of computePCAsynthetic
#'
#' @param listCatPop TODO
#'
#' @param study.id.syn a the study.id of the synthetic data
#'
#' @param spRef TODO
#'
#' @param fieldPopInfAnc TODO
#'
#' @param kList TODO array of the k possible values
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
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom class knn
#' @encoding UTF-8
#' @export
computeKNNRefSynthetic <- function(gdsSample, listEigenvector,
                                   listCatPop, study.id.syn,
                                   spRef, fieldPopInfAnc="SuperPop",
                                   kList = seq_len(15),
                                   pcaList = 2:15) {

    ## The number of rows in study.annot must be one.
    # if(nrow(study.annot) < 1) {
    #     stop("Number of samples in study.annot not equal to 1\n")
    # }

    if(is.null(kList)){
        kList <- seq_len(15)#c(seq_len(14), seq(15,100, by=5))
    }
    if(is.null(pcaList)){
        pcaList <- 2:15
    }

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
        for(pcaD in pcaList) {
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
                    sp=spRef[study.annot$case.id],
                    matKNN=resMat)

    return(listKNN)
}


#' @title Deprecated
#'
#' @description Deprecated
#'
#' @param listEigenvector TODO see return of computePCARefSample
#'
#' @param sample.ref TODO
#'
#' @param study.annot a  \code{data.frame} with one entry from study.annot in
#' the gds
#'
#' @param spRef TODO
#'
#' @param kList TODO array of the k possible values
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
#' @importFrom gdsfmt add.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASampLoading snpgdsPCASampLoading
#' @importFrom class knn
#' @encoding UTF-8
#' @keywords internal
computeKNNSuperPopSample <- function(gdsSample, listEigenvector, name.id,
                                    spRef, study.id.ref="Ref.1KG",
                                    kList=seq_len(15), pcaList=2:15) {

    if(is.null(kList)) {
        kList <- seq_len(15)#c(seq_len(14), seq(15,100, by=5))
    }
    if(is.null(pcaList)) {
        pcaList <- 2:15
    }
    if(length(name.id) != 1) {
        stop("Number of sample in study.annot not equal to 1\n")
    }

    #sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))


    study.annot.all <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    sample.ref <- study.annot.all[which(study.annot.all$study.id ==
                                                study.id.ref), "data.id"]

    resMat <- data.frame(sample.id=rep(listEigenvector$sample.id,
                                       length(pcaList) * length(kList)),
                         D=rep(0,length(pcaList) * length(kList)),
                         K=rep(0,length(pcaList) * length(kList)),
                         SuperPop=character(length(pcaList) * length(kList)),
                         stringsAsFactors=FALSE)

    listSuperPop <- c("EAS", "EUR", "AFR", "AMR", "SAS")

    #curPCA <- listPCA.Samples[[sample.id[sample.pos]]]
    eigenvect <- rbind(listEigenvector$eigenvector.ref,
                        listEigenvector$eigenvector)

    rownames(eigenvect) <- c(sample.ref,
                                listEigenvector$sample.id)

    totR <- 1
    for(pcaD in pcaList) {
        for(kV in  seq_len(length(kList))) {
            dCur <- paste0("d", pcaD)
            kCur <- paste0("k", kList[kV])
            resMat[totR,c("D", "K")] <- c(pcaD, kList[kV])

            pcaND <- eigenvect[ ,seq_len(pcaD)]
            y_pred <- knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                test=pcaND[rownames(eigenvect)[nrow(eigenvect)],, drop=FALSE],
                cl=factor(spRef[rownames(eigenvect)[-1*nrow(eigenvect)]],
                    levels=listSuperPop, labels=listSuperPop),
                k=kList[kV],
                prob=FALSE)

            resMat[totR, paste0("SuperPop")] <- listSuperPop[as.integer(y_pred)]

            totR <- totR + 1
        } # end k
    } # end pca Dim
    listKNN <- list(sample.id=listEigenvector$sample.id,
                        matKNN=resMat)

    return(listKNN)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param listEigenvector TODO see return of computePCARefSample
#'
#' @param listCatPop TODO
#'
#' @param spRef TODO
#'
#' @param fieldPopInfAnc TODO
#'
#' @param kList TODO array of the k possible values
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
computeKNNRefSample <- function(listEigenvector, listCatPop,
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
#                         SuperPop=character(length(pcaList) * length(kList)),
                         stringsAsFactors=FALSE)
    resMat[[fieldPopInfAnc]] <- character(length(pcaList) * length(kList))

    listSuperPop <- c("EAS", "EUR", "AFR", "AMR", "SAS")

    #curPCA <- listPCA.Samples[[sample.id[sample.pos]]]
    eigenvect <- rbind(listEigenvector$eigenvector.ref,
                       listEigenvector$eigenvector)

    # rownames(eigenvect) <- c(sample.ref,
    #                          listEigenvector$sample.id)

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

            resMat[totR, fieldPopInfAnc] <- listSuperPop[as.integer(y_pred)]

            totR <- totR + 1
        } # end k
    } # end pca Dim

    listKNN <- list(sample.id=listEigenvector$sample.id, matKNN=resMat)

    return(listKNN)
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
#' @param sampleRM TODO
#'
#' @param spRef TODO
#'
#' @param study.id.syn TODO
#'
#' @param np TODO
#'
#' @param listCatPop TODO
#'
#' @param fieldPopIn1KG TODO
#'
#' @param fieldPopInfAnc TODO
#'
#' @param kList TODO array of the k possible values
#'
#' @param pcaList TODO array of the pca dimension possible values
#'
#' @param algorithm algorithm of the PCA "exact", "randomized"
#' (para snpgdsPCA)
#'
#' @param eigen.cnt number of eigenvectors in PCA
#' (para snpgdsPCA)
#'
#' @param missing.rate number of eigenvectors in PCA
#' (para snpgdsPCA)
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
#' @encoding UTF-8
#' @export
computePoolSyntheticAncestryGr <- function(gds, gdsSample,
                            sampleRM, spRef,
                            study.id.syn,
                            np = 1L,
                            listCatPop = c("EAS", "EUR", "AFR", "AMR", "SAS"),
                            fieldPopIn1KG = "superPop",
                            fieldPopInfAnc = "SuperPop",
                            kList = seq(2,15,1),
                            pcaList = 2:15,
                            algorithm="exact",
                            eigen.cnt=32L,
                            missing.rate=0.025) {


    pca1KG <- computePCARefRMMulti(gdsSample, names(spRef),
                                    sampleRM, np=np,
                                    algorithm=algorithm,
                                    eigen.cnt=eigen.cnt,
                                    missing.rate=missing.rate)

    resPCA <- computePCAMultiSynthetic(gdsSample, pca1KG,
                                        sampleRM, study.id.syn)

    KNN.synt <- computeKNNRefSynthetic(gdsSample, resPCA,
                                        listCatPop,
                                        study.id.syn, spRef,
                                        fieldPopInfAnc=fieldPopInfAnc)

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
#' @param dataRef TODO
#'
#' @param spRef TODO
#'
#' @param study.id.syn TODO
#'
#' @param np TODO
#'
#' @param listCatPop TODO
#'
#' @param fieldPopIn1KG TODO
#'
#' @param fieldPopInfAnc TODO
#'
#' @param kList TODO array of the k possible values
#'
#' @param pcaList TODO array of the pca dimension possible values
#'
#' @param algorithm algorithm of the PCA "exact", "randomized"
#' (para snpgdsPCA)
#'
#' @param eigen.cnt number of eigenvectors in PCA
#' (para snpgdsPCA)
#'
#' @param missing.rate number of eigenvectors in PCA
#' (para snpgdsPCA)
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
                                         kList = seq(2,15,1),
                                         pcaList = 2:15,
                                         algorithm="exact",
                                         eigen.cnt=32L,
                                         missing.rate=0.025) {

    sampleRM <- splitSelectByPop(dataRef)
    KNN.list <- list()
    for(j in seq_len(nrow(sampleRM))){
        KNN.list[[j]] <- computePoolSyntheticAncestryGr(gds, gdsSample,
                                                        sampleRM[j,],
                                                        spRef, study.id.syn,
                                                        np, listCatPop,
                                                        fieldPopIn1KG,
                                                        fieldPopInfAnc,
                                                        kList,
                                                        pcaList,
                                                        algorithm,
                                                        eigen.cnt,
                                                        missing.rate)
    }

    KNN.sample.syn <- do.call(rbind, KNN.list)
    pedSyn <- prepPedSynthetic1KG(gds, gdsSample,
                                    study.id.syn, fieldPopIn1KG)



    listParaSample <- selParaPCAUpQuartile(KNN.sample.syn, pedSyn,
                                            fieldPopIn1KG, fieldPopInfAnc,
                                            listCatPop)

    listPCASample <- computePCARefSample(gdsSample=gdsSample,
                            name.id=sample.ana.id, study.id.ref="Ref.1KG",
                            np=np, algorithm=algorithm,
                            eigen.cnt=eigen.cnt)


    listKNNSample <- computeKNNSuperPopSample(gdsSample=gdsSample,
                                                sample.ana.id,
                                                spRef)


    return(listKNNSample)
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
#' @param listFiles TODO.
#'
#' @param sample.ana.id TODO
#'
#' @param spRef TODO
#'
#' @param study.id.syn TODO
#'
#' @param np a single \code{integer} representing the number of CPU to use.
#' Default: \code{1L}.
#'
#' @param listCatPop TODO,
#' Default: \code{c("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param fieldPopIn1KG TODO
#'
#' @param fieldPopInfAnc TODO
#'
#' @param kList TODO array of the k possible values
#'
#' @param pcaList TODO array of the pca dimension possible values
#'
#' @param algorithm algorithm of the PCA "exact", "randomized"
#' (para snpgdsPCA)
#'
#' @param eigen.cnt number of eigenvectors in PCA
#' (para snpgdsPCA)
#'
#' @param missing.rate number of eigenvectors in PCA
#' (para snpgdsPCA)
#'
#' @return a \code{list} TODO with the sample.id and eigenvectors
#' and a table with KNN callfor different K and pca dimension.
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
                            kList=seq(2,15,1),
                            pcaList=2:15,
                            algorithm="exact",
                            eigen.cnt=32L,
                            missing.rate=0.025) {

    KNN.list <- list()

    for(j in seq_len(length(listFiles))) {
        # We have to test if the file exist and format is OK
        KNN.list[[j]] <- readRDS(listFiles[j])
    }

    KNN.sample.syn <- do.call(rbind, KNN.list)

    pedSyn <- prepPedSynthetic1KG(gds, gdsSample,
                                    study.id.syn, fieldPopIn1KG)



    listParaSample <- selParaPCAUpQuartile(KNN.sample.syn, pedSyn,
                                            fieldPopIn1KG, fieldPopInfAnc,
                                            listCatPop)

    listPCASample <- computePCARefSample(gdsSample, sample.ana.id,
                                            study.id.ref="Ref.1KG", np=np,
                                            algorithm=algorithm,
                                            eigen.cnt=eigen.cnt)


    listKNNSample <- computeKNNRefSample(listPCASample,
                                            listCatPop, spRef,
                                            kList = seq(2, 15, 1),
                                            pcaList = seq(2, 15, 1))
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
#' @param matKNN.All TODO see it is rbind matKNN of the
#' computeKNNSuperPoprSynthetic return from group synthetic data
#'
#' @param pedCall TODO see return of prepPedSynthetic1KG
#'
#' @param refCall TODO column name in pedCall with the call
#'
#' @param predCall TODO column name in matKNN with the call
#'
#' @param listCall TODO array of the possible call
#'
#' @param kList TODO array of the k possible values
#'
#' @param pcaList TODO array of the pca dimension possible values
#'
#' @return a \code{list} containing 5 entries:
#' \itemize{
#' \item{dfPCA} {TODO}
#' \item{dfPop} {TODO}
#' \item{D} {TODO}
#' \item{K} {TODO}
#' \item{listD} {TODO}
#' }
#'
#' @examples
#'
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom stats mad median quantile
#' @encoding UTF-8
#' @export
selParaPCAUpQuartile <- function(matKNN.All, pedCall, refCall,
                                    predCall, listCall,
                                    kList = 3:15, pcaList = 2:15) {
    if(min(kList) < 3) {
        warning("A K smaller than 3 could not give robust results.\n")
    }
    tableSyn <- list()
    tableCall <- list()
    tableAUROC <- list()
    i <- 1

    for(D in pcaList) {
        matKNNCurD <- matKNN.All[which(matKNN.All$D == D ), ]
        listTMP <- list()
        listTMP.AUROC <- list()
        j <- 1
        for(K in kList) {
            matKNNCur <- matKNNCurD[which(matKNNCurD$K == K), ]
            res <- computeSyntheticConfMat(matKNNCur, pedCall, refCall,
                                            predCall, listCall)
            resROC <- computeSyntheticROC(matKNNCur, pedCall, refCall,
                                            predCall, listCall)

            df <- data.frame(D=D,
                                K=K,
                                AUROC.min=min(resROC$matAUROC.Call$AUC),
                                AUROC=resROC$matAUROC.All$ROC.AUC,
                                Accu.CM=res$matAccuracy$Accu.CM)

            listTMP[[j]] <- df
            listTMP.AUROC[[j]] <- resROC$matAUROC.Call
            j <- j + 1
        }
        df <- do.call(rbind, listTMP)

        tableCall[[i]] <- df
        tableAUROC[[i]] <- do.call(rbind, listTMP.AUROC)
        maxAUROC <- max(df[df$K %in% kList, "AUROC.min"])
        kMax <- df[df$K %in% kList & abs(df$AUROC.min-maxAUROC) < 1e-3, "K"]
        kV <- kMax[(length(kMax) + length(kMax)%%2)/2]
        dfPCA <- data.frame(D = D,
                           median = median(df[df$K %in% kList, "AUROC.min"]),
                           mad = mad(df[df$K %in% kList, "AUROC.min"]),
                           upQuartile = quantile(df[df$K %in% kList,
                                                        "AUROC.min"], 0.75),
                           K = kV)
        tableSyn[[i]] <- dfPCA
        i <- i + 1
    }

    dfPCA <- do.call(rbind, tableSyn)
    dfCall <- do.call(rbind, tableCall)
    dfAUROC <- do.call(rbind, tableAUROC)
    selD <- dfPCA$D[which.max(dfPCA$upQuartile)]
    selK <- dfPCA$K[which.max(dfPCA$upQuartile)]
    tmp <- max(dfPCA$upQuartile)
    listD <- dfPCA$D[which(abs(dfPCA$upQuartile - tmp) < 1e-3)]

    res <- list(dfPCA=dfPCA,
                dfPop=dfCall,
                dfAUROC=dfAUROC,
                D = selD,
                K = selK,
                listD = listD)
    return(res)
}



#' @title Validate input parameters for pruningSample() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{pruningSample}} function.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1 KG GDS file.
#'
#' @param method a \code{character} string that represents the method that will
#' be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function.
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
#' \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#'
#' @param slide.max.bp.v a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done in the \code{\link[SNPRelate]{snpgdsLDpruning}}
#' function.
#'
#' @param ld.threshold.v a single \code{numeric} value that represents the LD
#' threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used.
#'
#' @param verbose.v a \code{logicial} indicating if information is shown
#' during the process in the \code{\link[SNPRelate]{snpgdsLDpruning}}
#' function.
#'
#' @param chr a \code{character} string representing the chromosome where the
#' selected SNVs should belong. Only one chromosome can be handled. If
#' \code{NULL}, the chromosome is not used as a filtering criterion.
#'
#' @param minAF.SuperPop a single positive \code{numeric} representing the
#' minimum allelic frequency used to select the SNVs. If \code{NULL}, the
#' allelic frequency is not used as a filtering criterion.
#'
#' @param keepGDSpruned a \code{logicial} indicating if the information about
#' the pruned SNVs should be added to the GDS Sample file.
#'
#' @param PATHSAMPLEGDS a \code{character} string representing the directory
#' where the GDS Sample file will be created. The directory must exist.
#'
#' @param keepFile a \code{logical} indicating if RDS files containing the
#' information about the pruned SNVs must be created.
#'
#' @param PATHPRUNED a \code{character} string representing an existing
#' directory. The directory must exist.
#'
#' @param outPref a \code{character} string that represents the prefix of the
#' RDS files that will be generated. The RDS files are only generated when
#' the parameter \code{keepFile}=\code{TRUE}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validatePruningSample(gds=gds1KG, method="corr",
#'      sampleCurrent="TGCA_01", study.id="TCGA",
#'      listSNP=c("sr10103", "sr10202"), slide.max.bp.v=1000L,
#'      ld.threshold.v=0.008, np=1L, verbose.v=TRUE, chr=1,
#'      minAF.SuperPop=0.002, keepGDSpruned=TRUE, PATHSAMPLEGDS=data.dir,
#'      keepFile=FALSE, PATHPRUNED=data.dir, outPref="test")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validatePruningSample <- function(gds, method, sampleCurrent, study.id,
                        listSNP, slide.max.bp.v, ld.threshold.v, np, verbose.v,
                        chr, minAF.SuperPop, keepGDSpruned, PATHSAMPLEGDS,
                        keepFile, PATHPRUNED, outPref) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gds, name="gds")

    ## The parameter sampleCurrent must be a character string
    if (!(is.character(sampleCurrent))) {
        stop("The \'sampleCurrent\' parameter must be a character string.")
    }

    ## The parameter method must be a character string
    if(!(is.character(method))) {
        stop("The \'method\' parameter must be a character string.")
    }

    ## The parameter ld.threshold.v must be a single positive integer
    if (!(isSingleNumber(ld.threshold.v) && (ld.threshold.v >= 0.0))) {
        stop("The \'ld.threshold.v\' parameter must be a single positive ",
                "numeric value.")
    }

    ## The parameter slide.max.bp.v must be a single positive integer
    if (!(isSingleNumber(slide.max.bp.v) && (slide.max.bp.v >= 0.0))) {
        stop("The \'slide.max.bp.v\' parameter must be a single positive ",
                "numeric value.")
    }

    ## The parameter np must be a single positive integer
    if (!(isSingleNumber(np) && (np >= 0.0))) {
        stop("The \'np\' parameter must be a single positive numeric value.")
    }

    ## The parameter keepGDSpruned must be a logical
    validateLogical(logical=keepGDSpruned, "keepGDSpruned")

    ## The parameter PATHSAMPLEGDS must be a character string for existing path
    if (!(is.character(PATHSAMPLEGDS) && dir.exists(PATHSAMPLEGDS))) {
        stop("The \'PATHSAMPLEGDS\' parameter must be a character string ",
                "representing an existing directory.")
    }

    ## The parameter keepFile must be a logical
    validateLogical(logical=keepFile, "keepFile")

    ## The parameter PATHPRUNED must be a character string for existing path
    if (!(is.character(PATHPRUNED) && dir.exists(PATHPRUNED))) {
        stop("The \'PATHPRUNED\' parameter must be a character string ",
                "representing an existing directory.")
    }

    return(0L)
}


#' @title Validate input parameters for pruningSample() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{pruningSample}} function.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1KG GDS file.
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
#' threads.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations.
#'
#' @param fieldPopIn1KG TODO. Default: \code{"superPop"}.
#'
#' @param fieldPopInfAnc TODO. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.  If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016).
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA} only
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
#'
#' @return The function returns \code{0L} when successful.
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
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validateComputePoolSyntheticAncestryGr(gds=gds1KG,
#'      gdsSample=gdsSample, sampleRM="TGCA_01", spRef="TCGA",
#'      study.id.syn="TCGA", np=1L, listCatPop=c("AFR", "EAS", "SAS"),
#'      fieldPopIn1KG="SuperPop",  fieldPopInfAnc="Pop", kList=seq_len(3),
#'      pcaList=seq_len(10), algorithm="exact", eigen.cnt=12L,
#'      missing.rate=0.02)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateComputePoolSyntheticAncestryGr <- function(gds, gdsSample, sampleRM,
        spRef, study.id.syn, np, listCatPop, fieldPopIn1KG,
        fieldPopInfAnc, kList, pcaList, algorithm, eigen.cnt, missing.rate) {

    ## The gds and gdsSample must be objects of class "gds.class"
    validateGDSClass(gds, "gds")
    validateGDSClass(gdsSample, "gdsSample")

    ## The parameter sampleRM must be a single positive integer
    if(!(is.character(sampleRM))) {
        stop("The \'sampleRM\' parameter must be a vector of character ",
                "strings.")
    }

    ## The parameter study.id.syn must be a character string
    if(!(is.character(study.id.syn))) {
        stop("The \'study.id.syn\' parameter must be a character string.")
    }

    ## The parameter np must be a single positive integer
    if(!(isSingleNumber(np) && (np > 0))) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## The parameter listCatPop must be a character string
    if(!(is.character(listCatPop))) {
        stop("The \'listCatPop\' parameter must be a vector of ",
                "character strings.")
    }

    ## The parameters must be vectors of positive integers
    validatePositiveIntegerVector(kList, "kList")
    validatePositiveIntegerVector(pcaList, "pcaList")

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    ## The parameter eigen.cnt must be a single integer
    if(!(isSingleNumber(eigen.cnt))) {
        stop("The \'eigen.cnt\' parameter must be a single integer.")
    }

    ## The parameter missing.rate must be a single numer [0,1] or NaN
    if(!(((isSingleNumber(missing.rate) && missing.rate >= 0.0 &&
                    missing.rate <= 1.0)) || is.nan(missing.rate)))  {
        stop("The \'missing.rate\' parameter must be a single positive ",
                "numeric between zero and one or NaN.")
    }

    return(0L)
}


#' @title Validate input parameters for estimateAllelicFraction() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{estimateAllelicFraction}} function.
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
#' the length of the chromosomes.
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "DNA" and "RNA". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"DNA"}.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage.
#'
#' @param minProb a single \code{numeric} between 0 and 1 representing TODO.
#'
#' @param eProb a single \code{numeric} between 0 and 1 representing the
#' probability of sequencing error.
#'
#' @param cutOffLOH a single \code{numeric} representing the cutoff, in log,
#' for the homozygote score to assign a region as LOH.
#'
#' @param cutOffHomoScore a single \code{numeric} representing the cutoff, in
#' log, that the SNVs in a block are called homozygote by error.
#'
#' @param wAR a single positive \code{integer} representing the size-1 of
#' the window used to compute an empty box.
#'
#' @param cutOffAR a single \code{numeric} representing the cutoff, in
#' log score, that the SNVs in a gene are allelic fraction different from 0.5.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file.
#'  **This parameter is RNA specific.**
#'
#' @param block.id a \code{character} string corresponding to the block
#' identifier in \code{gdsRefAnnot}. **This parameter is RNA specific.**
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required libraries
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## Get chromosome length information
#' chrInfo <- integer(25L)
#' for(i in seq_len(22L)){ chrInfo[i] <- length(Hsapiens[[paste0("chr", i)]])}
#' chrInfo[23] <- length(Hsapiens[["chrX"]])
#' chrInfo[24] <- length(Hsapiens[["chrY"]])
#' chrInfo[25] <- length(Hsapiens[["chrM"]])
#'
#' ## The validatiion should be successful
#' RAIDS:::validateEstimateAllelicFraction(gds=gds1KG, gdsSample=gdsSample,
#'     sampleCurrent="Sample01", study.id="Synthetic", chrInfo=chrInfo,
#'     studyType="DNA", minCov=10L, minProb=0.03, eProb=0.002, cutOffLOH=10,
#'     cutOffHomoScore=11, wAR=2, cutOffAR=10, gdsRefAnnot=gds1KG,
#'     block.id="1")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateEstimateAllelicFraction <- function(gds, gdsSample, sampleCurrent,
        study.id, chrInfo, studyType, minCov, minProb, eProb, cutOffLOH,
        cutOffHomoScore, wAR, cutOffAR, gdsRefAnnot, block.id) {

    ## The gds and gdsSample must be objects of class "gds.class"
    validateGDSClass(gds, "gds")
    validateGDSClass(gdsSample, "gdsSample")

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

    ## The minProb and eProb must be single positive numeric between 0 and 1
    validateSingleRatio(minProb, "minProb")
    validateSingleRatio(eProb, "eProb")

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

    return(0L)
}


#' @title Validate input parameters for createStudy2GDS1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{createStudy2GDS1KG}} function.
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
#'
#' @param PATHSAMPLEGDS a \code{character} string representing the path to
#' the directory where the GDS Sample files will be created.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' gds1KG <- file.path(data.dir, "gds1KG.gds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyInfo <- data.frame(study.id="Pancreatic.WES",
#'                 study.desc="Pancreatic study",
#'                 study.platform="WES",
#'                 stringsAsFactors=FALSE)
#'
#' ## PED Study
#' ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
#'             Case.ID=c("TCGA-H01", "TCGA-H02"),
#'             Sample.Type=c("DNA", "DNA"),
#'             Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))
#'
#' ## The validatiion should be successful
#' RAIDS:::validateCreateStudy2GDS1KG(pedStudy=ped, fileNameGDS=gds1KG,
#'             batch=1, studyDF=studyInfo,
#'             listSamples=c("Sample_01", "Sample_02"),
#'             PATHSAMPLEGDS=data.dir, verbose=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateCreateStudy2GDS1KG <- function(pedStudy, fileNameGDS, batch, studyDF,
                                        listSamples, PATHSAMPLEGDS, verbose) {

    ## The PED study must have the mandatory columns
    if (!(all(c("Name.ID", "Case.ID", "Sample.Type", "Diagnosis", "Source")
                        %in% colnames(pedStudy)))) {
        stop("The PED study data frame is incomplete. ",
                "One or more mandatory columns are missing.")
    }

    ## The fileNameGDS must be a character string and the file must exists
    if (!(is.character(fileNameGDS) && (file.exists(fileNameGDS)))) {
        stop("The \'fileNameGDS\' must be a character string representing ",
                "the GDS 1KG file. The file must exist.")
    }

    ## The batch must be a single numeric
    if(!(isSingleNumber(batch))) {
        stop("The \'batch\' must be a single integer.")
    }

    ## The listSamples must be a vector of character string
    if (!(is.character(listSamples) || is.null(listSamples))) {
        stop("The \'listSamples\' must be a vector ",
                "of character strings (1 entry or more) or NULL.")
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, "verbose")

    return(0L)
}


#' @title Validate input parameters for createStudy2GDS1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{createStudy2GDS1KG}} function.
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
#' threads.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations.
#'
#' @param fieldPopIn1KG a \code{character} string representing the name of TODO
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA.
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA}
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
#'
#' @return The function returns \code{0L} when successful.
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
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' listFiles <- file.path(data.dir,  "listSNPIndexes_Demo.rds")
#'
#' ## The validatiion should be successful
#' RAIDS:::validateComputeAncestryFromSyntheticFile(gds=gds1KG,
#'     gdsSample=gdsSample, listFiles=listFiles, sample.ana.id="sample01",
#'     spRef=NULL, study.id.syn="Synthetic", np=1L, listCatPop=c("AFR", "EUR"),
#'     fieldPopIn1KG="superpop", fieldPopInfAnc="Superpop", kList=c(2, 3, 4),
#'     pcaList=c(3, 4, 5), algorithm="exact", eigen.cnt=32L, missing.rate=0.2)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateComputeAncestryFromSyntheticFile <- function(gds, gdsSample,
                listFiles, sample.ana.id, spRef, study.id.syn, np, listCatPop,
                fieldPopIn1KG, fieldPopInfAnc, kList, pcaList,
                algorithm, eigen.cnt, missing.rate) {

    ## The gds and gdsSample must be objects of class "gds.class"
    validateGDSClass(gds, "gds")
    validateGDSClass(gdsSample, "gdsSample")

    ## The parameter np must be a single positive integer
    if(!(isSingleNumber(np) && (np > 0))) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## The parameters are character strings (vector of 1 entry)
    validateCharacterString(value=fieldPopIn1KG, "fieldPopIn1KG")
    validateCharacterString(value=fieldPopInfAnc, "fieldPopInfAnc")
    validateCharacterString(value=study.id.syn, "study.id.syn")

    ## The parameter listCatPop must be a vector of of character strings
    if(!(is.character(listCatPop))) {
        stop("The \'listCatPop\' parameter must be a vector of ",
                "character strings.")
    }

    ## The parameters must be vectors of positive integers
    validatePositiveIntegerVector(kList, "kList")
    validatePositiveIntegerVector(pcaList, "pcaList")

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    ## The missing.rate must be a positive numeric between zero and one or NaN
    if (!is.nan(missing.rate)) {
        if (!(isSingleNumber(missing.rate) && (missing.rate >= 0.0) &&
                    (missing.rate <= 1.0))) {
            stop("The \'missing.rate\' must be a single numeric positive ",
                        "value between 0 and 1 or NaN.")
        }
    }

    return(0L)
}


#' @title Validate input parameters for computePCARefSample() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{computePCARefSample}} function.
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
#' that will be used.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA.
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; if \code{NaN}, no missing
#' threshold.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## The validatiion should be successful
#' RAIDS:::validateComputePCARefSample(gdsSample=gdsSample, name.id="HCC01",
#'     study.id.ref="1KG", np=1L, algorithm="exact", eigen.cnt=32L,
#'     missing.rate=0.02)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateComputePCARefSample <- function(gdsSample, name.id, study.id.ref,
                                            np, algorithm,
                                            eigen.cnt, missing.rate) {

    ## The gdsSample must be object of class "gds.class"
    validateGDSClass(gdsSample, "gdsSample")

    ## Validate that name.id is a string
    if(!(is.character(name.id) && length(name.id) == 1)) {
        stop("The \'name.id\' parameter must be a single character string.")
    }

    ## Validate that study.id.ref is a string
    if(!(is.character(study.id.ref) && length(name.id) == 1)) {
        stop("The \'study.id.ref\' parameter must be a character string.")
    }

    ## The parameter np must be a single positive integer
    if(!(isSingleNumber(np) && (np > 0))) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    ## The parameter eigen.cnt must be a single integer
    if(!(isSingleNumber(eigen.cnt))) {
        stop("The \'eigen.cnt\' parameter must be a single integer.")
    }

    ## The missing.rate must be a positive numeric between zero and one or NaN
    if (!is.nan(missing.rate)) {
        if (!(isSingleNumber(missing.rate) && (missing.rate >= 0.0) &&
                (missing.rate <= 1.0))) {
            stop("The \'missing.rate\' must be a single numeric positive ",
                    "value between 0 and 1 or NaN.")
        }
    }

    return(0L)
}

#' @title Validate input parameters for appendStudy2GDS1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{appendStudy2GDS1KG}} function.
#'
#' @param PATHGENO a \code{character} string representing the path to the
#' directory containing the output of SNP-pileup, a VCF Sample file, for
#' each sample.
#'
#' @param fileNamePED a \code{character} string representing the path to the
#' RDS file that contains the information about the sample to analyse.
#'
#' @param fileNameGDS a \code{character} string representing the file name of
#' the 1KG GDS file. The file must exist.
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
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
#'
#' @param PATHSAMPLEGDS a \code{character} string representing the path to
#' the directory where the GDS Sample files will be created.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' gds1KG <- file.path(data.dir, "1KG_Demo.gds")
#' ped <- file.path(data.dir, "unrelatedPatientsInfo_Demo.rds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyInfo <- data.frame(study.id="Pancreatic.WES",
#'                 study.desc="Pancreatic study",
#'                 study.platform="WES",
#'                 stringsAsFactors=FALSE)
#'
#' ## The validatiion should be successful
#' RAIDS:::validateAppendStudy2GDS1KG(PATHGENO=data.dir,
#'     fileNamePED=ped, fileNameGDS=gds1KG,
#'     batch=1L, studyDF=studyInfo, listSamples=c("HC01", "HC02"),
#'     PATHSAMPLEGDS=data.dir, verbose=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateAppendStudy2GDS1KG <- function(PATHGENO, fileNamePED, fileNameGDS,
                                batch, studyDF, listSamples, PATHSAMPLEGDS,
                                verbose) {

    ## The PATHGENO must be a character string and the path must exists
    if (!(is.character(PATHGENO) && (dir.exists(PATHGENO)))) {
        stop("The \'PATHGENO\' must be a character string representing ",
                "a path. The path must exist.")
    }

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

    if (!(is.data.frame(studyDF) && all(c("study.id", "study.desc",
                        "study.platform") %in% colnames(studyDF)))) {
        stop("The \'studyDF\' must be a data.frame and contain those 3 ",
            "columns: \'study.id\', \'study.desc\' and \'study.platform\'.")
    }

    ## The listSamples must be a vector of character string
    if (!(is.character(listSamples) || is.null(listSamples))) {
        stop("The \'listSamples\' must be a vector ",
                "of character strings (1 entry or more) or NULL.")
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, name="verbose")

    return(0L)
}

#' @title Validate input parameters for add1KG2SampleGDS() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{add1KG2SampleGDS}} function.
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
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The validatiion should be successful
#' RAIDS:::validateAdd1KG2SampleGDS(gds=gds1KG,
#'     gdsSampleFile=file.path(data.dir, "GDS_Sample_with_study_demo.gds"),
#'     sampleCurrent="Sample01", study.id="Synthetic")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateAdd1KG2SampleGDS <- function(gds, gdsSampleFile, sampleCurrent,
                                       study.id) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gds, name="gds")

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

    return(0L)
}


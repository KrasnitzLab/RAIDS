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

    ## The parameter kList must be positive integer values
    if(!(is.numeric(kList) && is.vector(kList) && all(kList>0))) {
        stop("The \'kList\' parameter must be a vector of positive ",
             "integers.")
    }

    ## The parameter pcaList must be positive integer values
    if(!(is.numeric(pcaList) && is.vector(pcaList) && all(pcaList>0))) {
        stop("The \'pcaList\' parameter must be a vector of positive ",
             "integers.")
    }

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    ## The parameter eigen.cnt must be a single integer
    if(!(isSingleNumber(eigen.cnt))) {
        stop("The \'eigen.cnt\' parameter must be a single integer.")
    }

    ## The parameter missing.rate must be a single positive numeric between
    ## zero and one or NaN
    if(!(((isSingleNumber(missing.rate) && missing.rate >= 0.0 &&
           missing.rate <= 1.0)) || is.nan(missing.rate)))  {
        stop("The \'missing.rate\' parameter must be a single positive ",
             "numeric between zero and one or NaN.")
    }

    return(0L)
}



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
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'.")
    }

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
    if (!is.logical(keepGDSpruned)) {
        stop("The \'keepGDSpruned\' parameter must be a logical ",
                "(TRUE or FALSE).")
    }

    ## The parameter PATHSAMPLEGDS must be a character string representing an
    ## existing path
    if (!(is.character(PATHSAMPLEGDS) && dir.exists(PATHSAMPLEGDS))) {
        stop("The \'PATHSAMPLEGDS\' parameter must be a character string ",
                "representing an existing directory.")
    }

    ## The parameter keepFile must be a logical
    if (!is.logical(keepFile)) {
        stop("The \'keepFile\' parameter must be a logical (TRUE or FALSE).")
    }

    ## The parameter PATHPRUNED must be a character string representing an
    ## existing path
    if (!(is.character(PATHPRUNED) && dir.exists(PATHPRUNED))) {
        stop("The \'PATHPRUNED\' parameter must be a character string ",
                "representing an existing directory.")
    }

    return(0L)
}

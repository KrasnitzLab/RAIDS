#' @title Validate input parameters for getTableSNV() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{getTableSNV}} function.
#'
#' @param gds an object of class \code{\link[gdsfmt]{gds.class}} (a GDS file),
#' the opened 1KG GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened GDS Sample file.
#'
#' @param sampleCurrent a \code{character} string corresponding to
#' the sample identifier used in \code{\link{pruningSample}} function.
#'
#' @param study.id a \code{character} string corresponding to the study
#' identifier used in \code{\link{pruningSample}} function.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage.
#'
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#'
#' @param eProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error.
#'
#' @param verbose a \code{logicial} indicating if messages should be printed
#' when the function is running.
#'
#' @return The integer \code{0L} when the function is successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The GDS Sample file
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                          "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validateGetTableSNV(gds=gds1KG, gdsSample=gdsSample,
#'      sampleCurrent="A101TCGA", study.id="TCGA", minCov=10L,
#'      minProb=0.998, eProb=0.002, verbose=TRUE)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateGetTableSNV <- function(gds, gdsSample, sampleCurrent, study.id,
                                    minCov, minProb, eProb, verbose) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gds, name="gds")

    ## The gdsSample must be an object of class "gds.class"
    validateGDSClass(gds=gdsSample, name="gdsSample")

    ## The minCov must be a single positive number
    if (!(isSingleNumber(minCov) && (minCov >= 0)))  {
        stop("The \'minCov\' must be a single numeric positive value.")
    }

    ## The minProb must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(minProb) && (minProb >= 0.0) && (minProb <= 1.0)))  {
        stop("The \'minProb\' must be a single numeric positive ",
                "value between 0 and 1.")
    }

    ## The eProb must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(eProb) && (eProb >= 0.0) && (eProb <= 1.0)))  {
        stop("The \'eProb\' must be a single numeric positive ",
                "value between 0 and 1.")
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, "verbose")

    ## Successful
    return(0L)
}

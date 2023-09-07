


#' @title Validate that the input parameter is a GDS object
#'
#' @description This function validates that the input parameter inherits
#' the \link[gdsfmt]{gds.class} class.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1 KG GDS file.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(dataDir, "1KG_Demo.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validateGDSClass(gds=gds1KG, name="gds")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateGDSClass <- function(gds, name) {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'", name ,"\' must be an object of class \'gds.class\'.")
    }

    return(0L)
}


#' @title Validate that the input parameter is a logical
#'
#' @description This function validates that the input parameter is
#' a logical. If the parameter is not a logicial, the function generates
#' an error with a specific message.
#'
#' @param logical a \code{logical} that will be validated.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## The validation should be successful
#' RAIDS:::validateLogical(logical=TRUE, name="test")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateLogical <- function(logical, name) {

    ## The parameter must be a logical
    if (!is.logical(logical)) {
        stop("The \'", name, "\' parameter must be a logical (TRUE or FALSE).")
    }

    return(0L)
}


#' @title Validate that the input parameter is a single positive numeric
#' between zero and one (included)
#'
#' @description This function validates that the input parameter is
#' a single numeric between zero and one (included). If the parameter is not,
#' the function generates an error with a specific message.
#'
#' @param value a single positive \code{numeric} that will be validated.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## The validation should be successful
#' RAIDS:::validateSingleRatio(value=0.02, name="test")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateSingleRatio <- function(value, name) {

    ## The parameter must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(value) && (value >= 0.0) && (value <= 1.0))) {
        stop("The \'", name, "\' must be a single numeric positive ",
                "value between 0 and 1.")
    }

    return(0L)
}


#' @title Validate that the input parameter is a character sting
#'
#' @description This function validates that the input parameter is
#' a character string (vector of 1 entry). If the parameter is not a
#' character string, the function generates
#' an error with a specific message.
#'
#' @param value a \code{character} string that will be validated.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## The validation should be successful
#' RAIDS:::validateCharacterString(value="hi", name="test")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateCharacterString <- function(value, name) {

    ## The parameter value must be a character string
    if(!(is.character(value) && length(value) == 1)) {
        stop("The \'", name, "\' parameter must be a character string.")
    }

    return(0L)
}


#' @title Validate that the input parameter is a vector of positive
#' \code{numeric}
#'
#' @description This function validates that the input parameter is
#' a \code{vector} of positive \code{numeric} values (vector of 1 entry or
#' more). All values have to be positive (>\code{0}).
#' If the parameter is not respecting the validation, the function generates
#' an error with a specific message.
#'
#' @param value a \code{vector} of \code{numeric} that will be validated.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## The validation should be successful
#' RAIDS:::validatePositiveIntegerVector(value=c(1, 2 ,3, 5), name="test")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validatePositiveIntegerVector <- function(value, name) {

    ## The parameter value must be a character string
    if(!(is.numeric(value) && is.vector(value) && all(value > 0))) {
        stop("The \'", name, "\' parameter must be a vector of positive ",
                    "integers.")
    }

    return(0L)
}

#' @title Read a SNP-pileup file
#'
#' @description The function reads a generic SNP pileup file and
#' returns a data frame
#' containing the information about the read counts for the SNVs present in
#' the file.
#'
#' @param fileName a \code{character} string representing the name, including
#' the path, of a text file containing the SNV read counts as generated by
#' snp-pileup software. The text file must be comma separated.
#' The text file must contain those columns: Chromosome, Position, Ref, Alt,
#' File1R, File1A, File1E and File1D.
#'
#' @param offset a \code{integer} representing the offset to be added to the
#' position of the SNVs. The value of offset
#' is added to the position present in the file. Default: \code{0L}.
#'
#' @return the a \code{data.frame} containing at least:
#' \itemize{
#' \item{Chromosome}{ a \code{numeric} representing the name of
#' the chromosome}
#' \item{Position}{ a \code{numeric} representing the position on the
#' chromosome}
#' \item{Ref}{ a \code{character} string representing the reference nucleotide}
#' \item{Alt}{ a \code{character} string representing the alternative
#' nucleotide}
#' \item{File1R} { a \code{numeric} representing the count for
#' the reference nucleotide}
#' \item{File1A} { a \code{numeric} representing the count for the
#' alternative nucleotide}
#' \item{File1E} {a \code{numeric} representing the count for the
#' errors}
#' \item{File1D} {a \code{numeric} representing the count for the
#' deletions}
#' \item{count} { a \code{numeric} representing the total count}
#' }
#'
#' @examples
#'
#' ## Directory where demo SNP-pileup file
#' dataDir <- system.file("extdata/example/snpPileup", package="RAIDS")
#'
#' ## The SNP-pileup file
#' snpPileupFile <- file.path(dataDir, "ex1.txt.gz")
#'
#' info <- RAIDS:::readSNVPileupFile(fileName=snpPileupFile)
#' head(info)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
readSNVPileupFile <- function(fileName, offset = 0L) {

    matSample <- read.csv(fileName)

    # Check if the mandatory column are there
    if (!(all(c("Chromosome", "Position", "Ref", "Alt", "File1R",
                "File1A", "File1E", "File1D") %in% colnames(matSample)))) {
        stop("The SNP-pileup file must contain all those columns: ",
                "\'Chromosome\', \'Position\', \'Ref\', \'Alt\', \'File1R\', ",
                "\'File1A\', \'File1E\', \'File1D\'.")
    }

    matSample[, "Chromosome"] <- as.integer(gsub("chr", "",
                                            matSample[, "Chromosome"]))
    matSample[, "Position"] <- matSample[, "Position"] + offset
    matSample[, "count"] <- rowSums(matSample[, c("File1R", "File1A",
                                            "File1E", "File1D")])

    return(matSample)
}


#' @title Read a generic SNP pileup file
#'
#' @description The function reads a generic SNP pileup file and
#' returns a data frame
#' containing the information about the read counts for the SNVs present in
#' the file.
#'
#' @param fileName a \code{character} string representing the name, including
#' the path, of a text file containing the SNV read counts. The text file must
#' be comma separated. The text file must
#' contain those columns: Chromosome, Position, Ref, Alt, Count,
#' File1R and File1A.
#'
#' @param offset a \code{integer} representing the offset to be added to the
#' position of the SNVs. The value of offset
#' is added to the position present in the file. Default: \code{0L}.
#'
#' @return a \code{data.frame} containing at least:
#' \itemize{
#' \item{Chromosome}{ a \code{numeric} representing the name of
#' the chromosome}
#' \item{Position}{ a \code{numeric} representing the position on the
#' chromosome}
#' \item{Ref}{ a \code{character} string representing the reference nucleotide}
#' \item{Alt}{ a \code{character} string representing the alternative
#' nucleotide}
#' \item{File1R} { a \code{numeric} representing the count for
#' the reference nucleotide}
#' \item{File1A} { a \code{numeric} representing the count for the
#' alternative nucleotide}
#' \item{count} { a \code{numeric} representing the total count}
#' }
#'
#' @examples
#'
#'
#' ## Directory where demo SNP-pileup file
#' dataDir <- system.file("extdata/example/snpPileup", package="RAIDS")
#'
#' ## The SNP-pileup file
#' snpPileupFile <- file.path(dataDir, "ex1.generic.txt.gz")
#'
#' info <- RAIDS:::readSNVFileGeneric(fileName=snpPileupFile)
#' head(info)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
readSNVFileGeneric <- function(fileName, offset = 0L) {

    matSample <- read.csv(fileName)

    # Check if the mandatory column are there
    if (!(all(c("Chromosome", "Position", "Ref", "Alt", "File1R",
                    "File1A", "Count") %in% colnames(matSample)))) {
        stop("The generic SNP pileup file must contain all those columns: ",
                "\'Chromosome\', \'Position\', \'Ref\', \'Alt\', \'File1R\', ",
                "\'File1A\', \'Count\'.")
    }

    matSample[, "Chromosome"] <- as.integer(gsub("chr", "",
                                        matSample[, "Chromosome"]))
    matSample[, "Position"] <- matSample[, "Position"] + offset
    colnames(matSample)[colnames(matSample) == "Count"] <- "count"

    return(matSample)
}




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
#' gds1KG <- openfn.gds(file.path(dataDir, "gds1KG.gds"), readonly=TRUE)
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



#' @title TODO
#'
#' @description TODO
#'
#' @param snp.keep TODO
#'
#' @param PATHBLOCK TODO
#'
#' @param superPop TODO
#'
#' @param chr TODO
#'
#'
#' @return the a \code{array} with the sample from pedDF keept
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
processBlockChr <- function(snp.keep, PATHBLOCK, superPop, chr) {

    blockChr <- read.delim(file.path(PATHBLOCK,
                    paste0("block.sp.", superPop, ".f0.05.chr", chr,
                                                ".blocks.det")), sep="")

    z <- cbind(c(blockChr$BP1, snp.keep, blockChr$BP2+1),
                    c(seq_len(nrow(blockChr)),
                    rep(0, length(snp.keep)), -1*seq_len(nrow(blockChr))))

    z <- z[order(z[,1]),]
    block.snp <- cumsum(z[,2])[z[,2] == 0]

    curStart <- 0
    activeBlock <- 0
    blockState <- 0
    block.inter <- rep(0, length(which(block.snp == 0)))
    k <- 1
    for(i in seq_len(length(block.snp))){
        if(block.snp[i] == 0){
            if(activeBlock == 1){
                if(snp.keep[i] - curStart >= 10000) {
                    blockState <- blockState - 1

                    curStart <- snp.keep[i]
                }
            } else{
                blockState <- blockState - 1
                curStart <- snp.keep[i]
                curStart <- snp.keep[i]
                activeBlock <- 1
            }
            block.inter[k] <- blockState
            k <- k + 1
        }else{
            activeBlock <- 0
        }
    }
    block.snp[block.snp == 0] <- block.inter

    return(block.snp)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param fileName Output from snp-pileup
#' must csv with the columns:
#' Chromosome,Position,Ref,Alt,File1R,File1A,File1E,File1D
#'
#' @param offset TODO
#'
#' @return the a \code{data.frame} containing at least:
#' \itemize{
#' \item{Chromosome} {TODO}
#' \item{Position} {TODO}
#' \item{File1R} {TODO}
#' \item{File1A} {TODO}
#' \item{count} {TODO}
#' }
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
readSNVPileupFile <- function(fileName, offset = 0L) {

    matSample <- read.csv(fileName)


    matSample[, "Chromosome"] <- as.integer(gsub("chr", "",
                                            matSample[, "Chromosome"]))
    matSample[, "Position"] <- matSample[, "Position"] + offset
    matSample[, "count"] <- rowSums(matSample[, c("File1R", "File1A",
                                            "File1E", "File1D")])

    return(matSample)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param fileName File name with the path to a
#' csv with at least the columns:
#' Chromosome,Position,Ref,Alt,Count,File1R,File1A
#' where Count is the deep at the position,
#' FileR is the deep of the reference allele, and
#' File1A is the deep of the specific alternative allele
#'
#' @param offset TODO
#'
#' @return a \code{data.frame} containing at least:
#' \itemize{
#' \item{Chromosome} {TODO}
#' \item{Position} {TODO}
#' \item{Ref}
#' \item{Alt}
#' \item{File1R} {deep of the reference allele}
#' \item{File1A} {deep of the alternative allele}
#' \item{count} {Total deep at the position}
#' }
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
readSNVFileGeneric <- function(fileName, offset = 0L) {

    matSample <- read.csv(fileName)
    # Check if the mendatory column are there

    matSample[, "Chromosome"] <- as.integer(gsub("chr", "",
                                        matSample[, "Chromosome"]))
    matSample[, "Position"] <- matSample[, "Position"] + offset
    colnames(matSample)[colnames(matSample) == "Count"] <- "count"

    return(matSample)
}

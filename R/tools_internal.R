


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
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
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
    if(!(is.numeric(value) && all(value > 0))) {
        stop("The \'", name, "\' parameter must be a vector of positive ",
                    "integers.")
    }

    return(0L)
}





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

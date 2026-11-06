#' @title Initialization of the section related to the profile
#' information in the GDS file
#'
#' @description This function initializes the section related to the profile
#' information in the GDS file. The information is extracted from
#' the \code{data.frame} passed to the function. The nodes "sample.id" and
#' "sample.annot" are created in the GDS file.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param dfPedReference a \code{data.frame} containing the information
#' related to the
#' samples. It must have those columns: "sample.id", "Name.ID", "sex",
#' "pop.group", "superPop" and "batch". All columns, except "sex" and batch",
#' are \code{character} strings. The "batch" and "sex" columns are
#' \code{integer}. The unique identifier
#' of this \code{data.frame} is the "Name.ID" column. The row names of the
#' \code{data.frame} must correspond to the identifiers present in the
#' "Name.ID" column.
#'
#' @param listSamples a \code{vector} of \code{character} string representing
#' the identifiers of the selected profiles. If \code{NULL}, all profiles are
#' selected. Default: \code{NULL}.
#'
#' @return a \code{vector} of \code{character} string with the identifiers of
#' the profiles saved in the GDS file.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Temporary GDS file in current directory
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_10.gds")
#'
#' ## Create and open the GDS file
#' tmpGDS  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create "sample.annot" node (the node must be present)
#' pedInformation <- data.frame(sample.id=c("sample_01", "sample_02"),
#'         Name.ID=c("sample_01", "sample_02"),
#'         sex=c(1,1),  # 1:Male  2: Female
#'         pop.group=c("ACB", "ACB"),
#'         superPop=c("AFR", "AFR"),
#'         project = c("1000 Genomes Project", "1000 Genomes Project"),
#'         batch=c(1, 1),
#'         stringsAsFactors=FALSE)
#'
#' ## The row names must be the sample identifiers
#' rownames(pedInformation) <- pedInformation$Name.ID
#'
#' ## Add information about 2 samples to the GDS file
#' RAIDS:::generateGDSRefSample(gdsReference=tmpGDS,
#'         dfPedReference=pedInformation, listSamples=NULL)
#'
#' ## Read sample identifier list
#' read.gdsn(index.gdsn(node=tmpGDS, path="sample.id"))
#'
#' ## Read sample information from GDS file
#' read.gdsn(index.gdsn(node=tmpGDS, path="sample.annot"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=tmpGDS)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
generateGDSRefSampleV2 <- function(gdsReference, dfPedReference,
                                 listSamples=NULL) {

    if(!(is.null(listSamples))){
        dfPedReference <- dfPedReference[listSamples,]
    }

    add.gdsn(node=gdsReference, name="sample.id",
             val=dfPedReference[, "Name.ID"])

    ## Create a data.frame containing the information form the samples
    samp.annot <- data.frame(sex=dfPedReference[, "sex"],
                             pop.group=dfPedReference[, "pop.group"],
                             superPop=dfPedReference[, "superPop"],
                             project=dfPedReference[, "project"],
                             batch=dfPedReference[, "batch"],
                             superPop.Geo=dfPedReference[, "superPop.Geo"],
                             superPop.Alt=dfPedReference[, "superPop.Alt"],
                             superPop.NoAdm=dfPedReference[, "superPop.NoAdm"],
                             stringsAsFactors=FALSE)

    ## Add the data.frame to the GDS object
    add.gdsn(node=gdsReference, name="sample.annot", val=samp.annot)

    ## Return the vector of saved samples
    return(dfPedReference[, "Name.ID"])
}

#' @title Extract the specified column from the 1KG GDS 'sample.ref' node
#' for the reference profiles (real ancestry assignation)
#'
#' @description The function extract the specified column for the 'sample.ref'
#' node present in the Reference GDS file. The column must be present in the
#' \code{data.frame} saved in the 'sample.ref' node. Only the information for
#' the reference profiles is returned. The values
#' represent the known ancestry assignation.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{vector} of \code{character} strings representing the content
#' of the extracted column for the 1KG GDS 'sample.ref' node. The values
#' represent the known ancestry assignation. The profile
#' identifiers are used as names for the \code{vector}.
#'
#' @examples
#' 
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' 
#' ## Path for gds reference files
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG.gds")
#' 
#' pRAIDS <- RAIDSparam(fileReferenceGDS=fileGDS)
#' 
#' spRef <- getRefSuperPop2(pRAIDS)
#' 
#' head(spRef)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @export
getRefSuperPop2 <- function(pRAIDS) {
    ## , popName="superPop"
    ## The gdsReference must be an object of class "gds.class"
    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))
    
    popName <- fieldPopInRef(pRAIDS)
    ## The popName is a character string
    if (!is.character(popName)) {
        stop("The \'popName\' parameter must be a single character string.")
    }

    sample.ref <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    dataRef <- read.gdsn(index.gdsn(gdsReference,
                                    "sample.annot"))[which(sample.ref == TRUE),]

    if(! popName %in% colnames(dataRef)) {
        stop("The population ", popName, " is not supported ",
             "(not found in the reference GDS file).")
    }

    dataRef <- dataRef[, popName]
    names(dataRef) <- read.gdsn(index.gdsn(node=gdsReference,
                                    "sample.id"))[which(sample.ref == TRUE)]
    closefn.gds(gdsReference)
    return(dataRef)
}


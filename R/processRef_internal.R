#' @title Add the reference version in reference gds
#'
#' @description TODO
#'
#' @param nameVer the a \code{character} representing the name of the
#' reference version
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{integer} success is 0.
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#' ## See computeAdmixtureProportions in ancestryAnalysis.R
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom philentropy JSD
#' @encoding UTF-8
#' @keywords internal
referenceVersion <- function(nameVer, pRAIDS=pRAIDS) {

    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS), readonly=FALSE)


    nameAdmixtureDesc <- paste0("reference.version")

    if(! nameAdmixtureDesc %in% ls.gdsn(gdsReference)){
        # add

        add.gdsn(gdsReference, name=nameAdmixtureDesc , val = nameVer)

    }else{
        # append
        append.gdsn(index.gdsn(gdsReference, nameAdmixtureDesc),
                    nameVer, check=TRUE)

    }

    snpgdsClose(gdsReference)
    return(0L)
}

#' @title Add the information for admixture.desc and admixture.ag
#'
#' @description TODO
#'
#'
#' @param agGroups a \code{vector} of \code{character} representing the
#' parameters
#'
#' @param admixtureDesc the a \code{data.frame} containing at least:
#' \describe{
#' \item{id}{ a \code{character} representing the name of the admixture analysis}
#' \item{description}{ a \code{character} representing the description of
#' the admixture analysis}
#' \item{K}{ a \code{integer} integer representing the number of ancestral groups},
#' \item{anchorId}{ a \code{integer} integer representing the position in 
#' the second dimension of anchor.ref}
#' }
#'
#' @param pRAIDS a \code{parametersRAIDS} an object with all the RAIDS
#' parameters
#'
#' @return \code{integer} success is 0.
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#' ## See computeAdmixtureProportions in ancestryAnalysis.R
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom philentropy JSD
#' @encoding UTF-8
#' @keywords internal
admixtureMeta <- function(agGroups, admixtureDesc, pRAIDS=pRAIDS) {
    ## , popName="superPop"
    ## The gdsReference must be an object of class "gds.class"

    
    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS), readonly=FALSE)

    sampleRef <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[which(sampleRef==1)]

    # list of ancestry groups
    #
    # we need to add admixture description
    nameAdmixtureDesc <- paste0("admixture.desc")
    k <- length(agGroups)
    pos <- 1
    df <- data.frame(id.seq = as.integer(pos),
                     id = admixtureDesc$id[1],
                     description = admixtureDesc$description[1],
                     K = as.integer(k),
                     anchor.id = as.integer(admixtureDesc$anchorId),
                     stringsAsFactors = FALSE)
    if(! nameAdmixtureDesc %in% ls.gdsn(gdsReference)){
        # add

        add.gdsn(gdsReference, name=nameAdmixtureDesc , val = df,
                 compress="LZMA_RA")

    }else{
        # append

        pos <- max(read.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/id.seq")))) + 1
        df$id.seq <- as.integer(pos)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/id.seq")), compress="")
        append.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/id.seq")),
                    df$id.seq[1], check=TRUE)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/id.seq")), compress="LZMA_RA")
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/id")), compress="")
        append.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/id")),
                    df$id[1], check=TRUE)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/id")), compress="LZMA_RA")
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/description")), compress="")
        append.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/description")),
                    df$description[1], check=TRUE)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/description")), compress="LZMA_RA")
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/K")), compress="")
        append.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/K")),
                    df$K[1], check=TRUE)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/K")), compress="LZMA_RA")
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/anchor.id")), compress="")
        append.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/anchor.id")),
                    df$anchor.id[1], check=TRUE)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureDesc, "/anchor.id")), compress="LZMA_RA")

    }
    readmode.gdsn(index.gdsn(gdsReference, nameAdmixtureDesc))
    # we need pop which are represented by the column

    nameAdmixtureAG <- paste0("admixture.ag")
    dfAG <- data.frame(id.seq = as.integer(rep(pos, length(agGroups))),
                       column=seq_len(length(agGroups)),
                       ag = agGroups,
                       stringsAsFactors = FALSE)
    if(! nameAdmixtureAG %in% ls.gdsn(gdsReference)){
        # add
        add.gdsn(gdsReference, name=nameAdmixtureAG , val = dfAG,
                 compress="LZMA_RA")
    }else{
        # append
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/id.seq")), compress="")
        append.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/id.seq")),
                    dfAG$id.seq, check=TRUE)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/id.seq")), compress="LZMA_RA")
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/column")), compress="")
        append.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/column")),
                    dfAG$column, check=TRUE)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/column")), compress="LZMA_RA")
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/ag")), compress="")
        append.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/ag")),
                    dfAG$ag, check=TRUE)
        compression.gdsn(index.gdsn(gdsReference, paste0(nameAdmixtureAG, "/ag")), compress="LZMA_RA")

    }
    readmode.gdsn(index.gdsn(gdsReference, nameAdmixtureAG))

    snpgdsClose(gdsReference)
    return(0L)
}

#' @title Add the information for admixture.desc and admixture.ag
#'
#' @description TODO
#'
#'
#' @param matQ a \code{matrix} of \code{numeric} representing the
#' Q matrix from ADMIXTURE
#'
#' @param admixId the a \code{character} representing the name of the
#' admixture analysis
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{integer} success is 0.
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#' ## See computeAdmixtureProportions in ancestryAnalysis.R
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom philentropy JSD
#' @encoding UTF-8
#' @keywords internal
admixtureMatrix <- function(matQ, admixId, pRAIDS) {

    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS), readonly=FALSE)

    sampleRef <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[which(sampleRef==1)]
    
    if(length(sampleId) != nrow(matQ)){
        snpgdsClose(gdsReference)
        stop(paste0("Error sample.id not same size as matQ\n"))
    }
    admixPos <- -1
    if("admixture.desc" %in% ls.gdsn(gdsReference)){
        admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
        tmp <- which(admixDesc$id == admixId)
        if(length(tmp) == 1){
            admixPos <- admixDesc$id.seq[tmp]
        }else{
            snpgdsClose(gdsReference)
            stop(paste0("Error ", admixId," is not in admixture.desc\n"))
        }

    }else{
        snpgdsClose(gdsReference)
        stop(paste0("Error there is no field admixture.desc\n"))
    }

    # admixture matrix


    nameAdmixtureMatrix <- paste0("admixture.", ncol(matQ))
    if(! nameAdmixtureMatrix %in% ls.gdsn(gdsReference)){
        if(admixPos != 1){
            snpgdsClose(gdsReference)
            stop(paste0("Error the seq.id ", admixId, " doesn't match the posistion in the ", nameAdmixtureMatrix, " \n"))
        }
        mat3d <- array(c(unlist(matQ)), dim=c(nrow(matQ),ncol(matQ),1))
        add.gdsn(gdsReference, name=nameAdmixtureMatrix, val = mat3d, storage="float32",
                 compress="LZMA_RA")

    }else{
        admixMat <- index.gdsn(gdsReference, nameAdmixtureMatrix)
        # the number of 3th dimension of the matrix must be  admixPos-1
        # the new matrix will be append to the end at the positon admixPos
        if(objdesp.gdsn(admixMat)$dim[3] != (admixPos-1)){
            snpgdsClose(gdsReference)
            stop(paste0("Error the seq.id ", admixId, " doesn't match the posistion in the ", nameAdmixtureMatrix, " \n"))
        }
        mat3d <- array(c(unlist(matQ)), dim=c(nrow(matQ),ncol(matQ),1))
        compression.gdsn(admixMat, compress="")
        append.gdsn(admixMat,
                    mat3d, check=TRUE)
        compression.gdsn(admixMat, compress="LZMA_RA")

    }
    readmode.gdsn(index.gdsn(gdsReference, nameAdmixtureMatrix))

    snpgdsClose(gdsReference)
    return(0L)
}

#' @title Add the information for admixture.desc and admixture.ag
#'
#' @description TODO
#'
#'
#' @param anchor a \code{vector} of \code{character} representing the
#' sample of the reference with the names of the ancestral group when
#' the sample is considered not quasi admix.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{integer} success is 0.
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#' ## See computeAdmixtureProportions in ancestryAnalysis.R
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom philentropy JSD
#' @encoding UTF-8
#' @keywords internal
admixtureAnchor <- function(anchor, admixId, pRAIDS=pRAIDS) {
    ## , popName="superPop"
    ## The gdsReference must be an object of class "gds.class"

    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS), readonly=FALSE)

    sampleRef <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[which(sampleRef==1)]
    
    if(length(sampleId) != nrow(anchor)){
        snpgdsClose(gdsReference)
        stop(paste0("Error sample.id not same size as anchor\n"))
    }
    admixPos <- -1
    if("admixture.desc" %in% ls.gdsn(gdsReference)){
        admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
        tmp <- which(admixDesc$id == admixId)
        if(length(tmp) == 1){
            admixPos <- admixDesc$anchor.id[tmp]
        }else{
            snpgdsClose(gdsReference)
            stop(paste0("Error ", admixId," is not in admixture.desc\n"))
        }

    }else{
        snpgdsClose(gdsReference)
        stop(paste0("Error there is no field admixture.desc\n"))
    }


    # anchor matching the admixtureMatrix

    nameAdmixtureAnchor <- paste0("anchor.ref")
    if(! nameAdmixtureAnchor %in% ls.gdsn(gdsReference)){
        if(admixPos != 1){
            snpgdsClose(gdsReference)
            stop(paste0("Error the seq.id ", admixId, " doesn't match the posistion in the ", nameAdmixtureAnchor, " \n"))
        }
        matAnchor <- array(c(unlist(anchor)), dim=c(nrow(anchor),1))
        add.gdsn(gdsReference, name=nameAdmixtureAnchor, val = matAnchor, storage="string",
                 compress="LZMA_RA")

    }else{
        admixMat <- index.gdsn(gdsReference, nameAdmixtureAnchor)
        if(objdesp.gdsn(admixMat)$dim[2] != (admixPos-1)){
            snpgdsClose(gdsReference)
            stop(paste0("Error the seq.id ", admixId, " doesn't match the posistion in the ", nameAdmixtureAnchor, " \n"))
        }
        compression.gdsn(admixMat, compress="")
        matAnchor <- array(c(unlist(anchor)), dim=c(nrow(anchor),1))
        append.gdsn(admixMat, val=matAnchor)
        compression.gdsn(admixMat, compress="LZMA_RA")
    }
    readmode.gdsn(index.gdsn(gdsReference, nameAdmixtureAnchor))

    snpgdsClose(gdsReference)
    return(0L)
}

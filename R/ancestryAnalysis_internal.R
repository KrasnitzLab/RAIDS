#' @title Create the Profile GDS file(s) for one or multiple specific profiles
#' using the information from a RDS Sample description file and the 1KG
#' GDS file
#'
#' @description The function uses the information for the Reference GDS file
#' and the RDS Sample Description file to create the Profile GDS file. One
#' Profile GDS file is created per profile. One Profile GDS file will be
#' created for each entry present in the \code{listProfiles} parameter.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return The function returns a \code{RAIDSparam} an object with all the RAIDS
#' parameters updated.
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG.gds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyDF <- data.frame(study.id = "MYDATA",
#'                         study.desc = "Description",
#'                         study.platform = "PLATFORM",
#'                         stringsAsFactors = FALSE)
#'
#' ## The data.frame containing the information about the samples
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' samplePED <- data.frame(Name.ID=c("ex1"),
#'                     Case.ID=c("Patient_h11"),
#'                     Diagnosis=rep("Cancer"),
#'                     Sample.Type=c("Primary Tumor"),
#'                     Source=c("Databank B"), stringsAsFactors=FALSE,
#'                     drop=FALSE)
#' rownames(samplePED) <- samplePED$Name.ID
#'
#' pRAIDS <- RAIDS:::RAIDSparam(profileFile=file.path(dataDir, "ex1.txt.gz"),
#'     studyDF=studyDF,
#'     genoSource="snp-pileup",
#'     pedStudy=samplePED,
#'     pathProfileGDS=tempdir(),
#'     fileReferenceGDS=fileGDS
#'     )
#' ## Create the Profile GDS File for samples in 'listSamples' vector
#' ## (in this case, samples "ex1")
#' ## The Profile GDS file is created in the pathProfileGDS directory
#' pRAIDS <- RAIDS:::createProfile2(pRAIDS)
#'
#'
#' ## The Profile GDS file 'ex1.gds' has been created in the
#' ## specified directory
#' dir(pathProfileGDS(pRAIDS), "ex1.gds")
#'
#' ## Remove Profile GDS file (created for demo purpose)
#' unlink(file.path(pathProfileGDS(pRAIDS), "ex1.gds"), force=TRUE)
#' 
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @keywords internal

createProfile2 <- function(pRAIDS) {
    # profileFile, profileName,
    # filePedRDS=NULL, pedStudy=NULL, fileNameGDS,
    # batch=1, studyDF, listProfiles=NULL,
    # pathProfileGDS=NULL,
    # genoSource=c("snp-pileup", "generic", "VCF", "bam"),
    # paramProfile=list(ScanBamParam=NULL,
    #                   PileupParam=NULL,
    #                   yieldSize=10000000),
    # verbose=FALSE
    ## When filePedRDS is defined and pedStudy is null
    if (is.null(pedStudy(pRAIDS))) {
        stop("The parameter \'pedStudy\' must ",
             "be defined.")
    }

    ## Validate input parameters
    # validatecreateProfile( pedStudy=pedStudy(pRAIDS),
    #                        fileNameGDS=fileReferenceGDS(pRAIDS), batch=batch(pRAIDS), studyDF=studyDF(pRAIDS),
    #                        listProfiles=listProfiles(pRAIDS), pathProfileGDS=pathProfileGDS(pRAIDS),
    #                        genoSource=genoSource(pRAIDS), verbose=verbose(pRAIDS))

    

    ## Read the Reference GDS file
    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))

    ## Extract the chromosome and position information for all SNPs in 1KG GDS

    listPos <- NULL
    if( is.null(listPos(pRAIDS))){
        if(genoSource(pRAIDS) == "bam"){
            alDf <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))
            alDf <- matrix(unlist(strsplit(alDf,"\\/")),nrow=2)
            listPos(pRAIDS) <- data.frame(snp.chromosome = read.gdsn(index.gdsn(gdsReference, "snp.chromosome")),
                                  snp.position = read.gdsn(index.gdsn(gdsReference, "snp.position")),
                                  REF = alDf[1,],
                                  ALT = alDf[2,],
                                  stringsAsFactors = FALSE
            )
            # listChr <- unique(listPos$chr)
            # We can optimize
            # listPos <- lapply(listChr,
            #                 FUN=function(x, varDf){
            #                     return(varDf[which(varDf$chr == x),])
            #                 },
            #                 varDf=listPos)
            # names(listPos) <- paste0("chr", listChr)
            rm(alDf)
        } else{
            listPos(pRAIDS) <- data.frame(snp.chromosome=read.gdsn(index.gdsn(node=gdsReference, "snp.chromosome")),
                                  snp.position=read.gdsn(index.gdsn(node=gdsReference, "snp.position")))
        }
        
    }
    ## Create a data.frame containing the information

    # Need to reformat for bam in varDf

    if(verbose(pRAIDS)) {
        message("Start ", Sys.time())
        message("Sample info DONE ", Sys.time())
    }
    pRAIDS <- generateProfileRawGDS2(pRAIDS=pRAIDS)
    generateProfileGenoCall(pRAIDS=pRAIDS)
    # generateProfileGDS2(pRAIDS=pRAIDS)

    if(verbose(pRAIDS)) {
        message("Genotype DONE ", Sys.time())
    }

    ## Close 1KG GDS file
    closefn.gds(gdsReference)

    ## Return successful code
    return(pRAIDS)
}


#' @title Calculate Principal Component Analysis (PCA) on SNV genotype data set
#'
#' @description The functions calculates the principal component analysis (PCA)
#' for a list of pruned SNVs present in a Profile GDS file. The
#' \link[SNPRelate]{snpgdsPCA} function is used to do the calculation.
#'
#' @param listRM a \code{vector} of \code{character} strings containing the
#' identifiers for the reference samples that need to be removed for the
#' PCA analysis.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 2 entries:
#' \describe{
#' \item{pruned}{ a \code{vector} of SNV identifiers specifying selected SNVs
#' for the PCA analysis.}
#' \item{pca.unrel}{ a \code{snpgdsPCAClass} object containing the eigenvalues
#' as generated by \link[SNPRelate]{snpgdsPCA} function.}
#' }
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
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#' 
#' ## The Reference GDS file
#' path1KG <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#' 
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#' 
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#' 
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#' 
#' listPCARef <- RAIDS:::computePCARefRMMulti2(listRM = NULL, pRAIDS = pRAIDS)
#' 
#' names(listPCARef)
#' 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn openfn.gds
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASNPLoading
#' @encoding UTF-8
#' @keywords internal
computePCARefRMMulti2 <- function(listRM, pRAIDS) {
    # gdsProfile, refProfileIDs, listRM, np=1L,
    # algorithm="exact", eigenCount=32L,
    # missingRate=0.025, verbose
    # unrelatedSamples <- refProfileIDs[which(!(refProfileIDs %in% listRM))]

    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                                    paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(filename=fileGDSProfile)
    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))
    unrelatedSamples <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[
            read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
    
    if(! is.null(listRM)){
        unrelatedSamples <- unrelatedSamples[which(!(unrelatedSamples %in% 
                                                                    listRM))]
    }
    listPCA <- list()

    listPCA[["pruned"]] <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))

    ## Calculate Principal Component Analysis (PCA) on SNV genotype dataset
    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsobj=gdsReference,
                                        sample.id=unrelatedSamples,
                                        snp.id=listPCA[["pruned"]],
                                        num.thread=np(pRAIDS),
                                        missing.rate=PCAmissingRate(pRAIDS),
                                        algorithm=PCAalgorithm(pRAIDS),
                                        eigen.cnt=eigenCount(pRAIDS),
                                        verbose=verbose(pRAIDS))
    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                                    gdsobj=gdsReference,
                                                    num.thread=np(pRAIDS),
                                                    verbose=verbose(pRAIDS))
    closefn.gds(gdsReference)
    closefn.gds(gdsProfile)
    return(listPCA)
}

#' @title Calculate Principal Component Analysis (PCA) on SNV genotype data set
#'
#' @description The functions calculates the principal component analysis (PCA)
#' for a list of pruned SNVs present in a Profile GDS file. The
#' \link[SNPRelate]{snpgdsPCA} function is used to do the calculation.
#' 
#' @param listPCA a \code{list} containing the PCA \code{object} generated
#' with the 1KG reference profiles.
#'
#' @param profileId a \code{character} string representing the profile
#' identifier.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 2 entries:
#' \describe{
#' \item{pruned}{ a \code{vector} of SNV identifiers specifying selected SNVs
#' for the PCA analysis.}
#' \item{pca.unrel}{ a \code{snpgdsPCAClass} object containing the eigenvalues
#' as generated by \link[SNPRelate]{snpgdsPCA} function.}
#' }
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
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#' 
#' ## The Reference GDS file
#' path1KG <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#' 
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#' 
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#' 
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#' 
#' listPCARef <- RAIDS:::computePCARefRMMulti2(listRM = NULL, pRAIDS = pRAIDS)
#' listPCAProfile <- RAIDS:::computePCAProfile(listPCA=listPCARef, 
#'                       profileId=pedStudy(pRAIDS)$Name.ID[1], 
#'                       pRAIDS=pRAIDS)
#' 
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASampLoading snpgdsOpen snpgdsClose
#' @encoding UTF-8
#' @keywords internal
computePCAProfile <- function(listPCA, profileId, pRAIDS) {

    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                                    paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(filename=fileGDSProfile)
    listPCA[["samp.load"]] <- snpgdsPCASampLoading(listPCA[["snp.load"]],
            gdsobj=gdsProfile, sample.id=profileId,
            num.thread=np(pRAIDS), verbose=verbose(pRAIDS))

    rownames(listPCA[["pca.unrel"]]$eigenvect) <-
        listPCA[["pca.unrel"]]$sample.id
    rownames(listPCA[["samp.load"]]$eigenvect) <-
        listPCA[["samp.load"]]$sample.id

    listRes <- list(sample.id=listPCA[["samp.load"]]$sample.id,
                    eigenvector.ref=listPCA[["pca.unrel"]]$eigenvect,
                    eigenvector=listPCA[["samp.load"]]$eigenvect)

    snpgdsClose(gdsProfile)
    return(listRes)
}


#' @title Run a k-nearest neighbors analysis on a subset of the
#' synthetic dataset
#'
#' @description The function runs k-nearest neighbors analysis on a
#' subset of the synthetic data set. The function uses the 'knn' package.
#'
#' @param listEigenvector a \code{list} with 3 entries:
#' 'sample.id', 'eigenvector.ref' and 'eigenvector'. The \code{list} represents
#' the PCA done on the 1KG reference profiles and the synthetic profiles
#' projected onto it.
#'
#' @param K  a \code{integer} representing value tested for the  K parameter.
#' The K parameter represents the
#' number of neighbors used in the K-nearest neighbors analysis.
#'
#' @param D a \code{integer} representing  
#' value tested for the  D parameter. The D parameter represents the
#' number of dimensions used in the PCA analysis.
#' 
#' @param dfRef \code{data.frame} of \code{character}
#' a subset of \code{sample.annot} of the reference \code{GDS} 
#' the strings representing the known super population ancestry 
#' for a subset of the reference profiles. The 1KG profile  
#' identifiers are used as names for the \code{vector}. It must
#' be indexed by the sample.id.
#' 
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{vector} of \code{character} strings
#' representing the identifiers of the synthetic profiles analysed.}
#' \item{\code{sample1Kg}}{ a \code{vector} of \code{character} strings
#' representing the identifiers of the 1KG reference profiles used to
#' generate the synthetic profiles.}
#' \item{\code{sp}}{ a \code{vector} of \code{character} strings representing
#' the known super population ancestry of the 1KG reference profiles used
#' to generate the synthetic profiles.}
#' \item{\code{matKNN}}{ a \code{data.frame} containing the super population
#' inference for each synthetic profiles for different values of PCA
#' dimensions \code{D} and k-neighbors values \code{K}. The fourth column title
#' corresponds to the \code{fieldPopInfAnc} parameter.
#' The \code{data.frame} contains 4 columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing
#' the identifier of the synthetic profile analysed.}
#' \item{\code{D}}{ a \code{numeric} strings representing
#' the value of the PCA dimension used to infer the super population.}
#' \item{\code{K}}{ a \code{numeric} strings representing
#' the value of the k-neighbors used to infer the super population.}
#' \item{\code{fieldPopInfAnc} value}{ a \code{character} string representing
#' the inferred ancestry.}
#' }
#' }
#' }
#'
#' @examples
#' 
#' ## Required library
#' library(SNPRelate)
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' #' ## The Reference GDS file
#' path1KG <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#' gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))
#' dfRef <- read.gdsn(index.gdsn(gdsReference, "sample.annot"))
#' row.names(dfRef) <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
#' if(! is.null(sampleRef(pRAIDS))){
#'    dfRef <- dfRef[sampleRef(pRAIDS),]
#' }
#'
#' snpgdsClose(gdsReference)
#'
#' ## The PCA with 1 profile projected on the 1KG reference PCA
#' ## Only one profile is retained
#' pca <- demoPCASyntheticProfiles
#' pca$sample.id <- pca$sample.id[1]
#' pca$eigenvector <- pca$eigenvector[1, , drop=FALSE]
#'
#' ## Projects profile on 1KG PCA
#' results <- RAIDS:::computeKNNProfileSubSet(listEigenvector=pca,
#'     K=3,D=3, dfRef=dfRef, pRAIDS=pRAIDS)
#' head(results)
#' 
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom class knn
#' @encoding UTF-8
#' @keywords internal
computeKNNProfileSubSet <- function(listEigenvector, K, D, dfRef, pRAIDS) {

    
    # spRef <- read.gdsn(index.gdsn(gdsReference, "sample.annot/superPop"))
    # names(spRef) <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
    
    listCatPop <- unique(dfRef$superPop)
    # closefn.gds(gdsReference)
    listEigenvector$eigenvector.ref <- listEigenvector$eigenvector.ref[row.names(listEigenvector$eigenvector.ref) %in% row.names(dfRef),]

    resMat <- lapply(seq_len(length(listEigenvector$sample.id)),
                     FUN=function(i, listEigenvector, D, K, dfRef){
                         eigenvect <- rbind(listEigenvector$eigenvector.ref,
                                            listEigenvector$eigenvector[i,,drop=FALSE])

                         pcaND <- eigenvect[ ,seq_len(D)]
                         yPred <-
                             knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                                 test=pcaND[rownames(eigenvect)[nrow(eigenvect)],,
                                            drop=FALSE],
                                 cl=factor(dfRef[rownames(eigenvect)[-1*nrow(eigenvect)], "superPop"],
                                           levels=listCatPop, labels=listCatPop),
                                 k=K,
                                 prob=FALSE)


                         df <- data.frame(sample.id=listEigenvector$sample.id[i],
                                          D=D,
                                          K=K,
                                          stringsAsFactors=FALSE)
                         df[[fieldPopInfAnc(pRAIDS)]] <- listCatPop[as.integer(yPred)]
                         return(df)
                     },
                     listEigenvector=listEigenvector,
                     D=D,
                     K=K,
                     dfRef=dfRef
    )
    resMat <- do.call(rbind, resMat)
    return(resMat)
}


#' @title Create the P-matrices
#' 
#' @description TODO
#'
#' @param pathOut temporary path where the directory matP use to generate the P-matrix
#' 
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{vector} with the info of the index associated the synthetic
#' group with the matrix-P 
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' pathOut <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathOut)){
#'     dir.create(pathOut)
#' }
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#' matPIndex <- RAIDS:::computePmatrix(pathOut=pathOut, pRAIDS=pRAIDS)
#'
#' dir(file.path(pathOut,"matP"))
#'
#' o <- lapply(dir(file.path(pathOut,"matP")),
#'    FUN=function(x,pathC){unlink(file.path(pathC,x))}, 
#'    pathC = file.path(pathOut,"matP"))
#'
#' unlink(file.path(pathOut,"matP"), recursive=TRUE)
#' unlink(file.path(pathOut), recursive=TRUE)
#'
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @keywords internal
computePmatrix <- function(pathOut, pRAIDS) {
    
    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                            paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    
    
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    pruned <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    snpgdsClose(gdsProfile)
    
    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))
    
    admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
    pos <- admixDesc[admixDesc$id == anchorId(pRAIDS), "id.seq"]

    anchor <- read.gdsn(index.gdsn(gdsReference, "anchor.ref"),
                        start=c(1, pos), count=c(-1,1))
    snpgdsClose(gdsReference)
    
    # Bad order but ok
    matGr <- getMatrixPopSynthetic(pRAIDS)

    rowF <- prepPMatrix( anchor, snpId = pruned, matGr, pRAIDS)

    pathBaseAdmix <- file.path(pathOut,"matP")
    if(!dir.exists(file.path(pathBaseAdmix))){
        dir.create(file.path(pathBaseAdmix))
    }

    matPIndex <- writePMatrix(pathBaseAdmix, rowF, pRAIDS)

    return(matPIndex)
}


#' @title Generate the Plink files (BED, BIM, and FAM)
#' 
#' @description TODO
#'
#' @param pathOut a \code{character} representing temporary path 
#' where the directory profileAdmixture is create and where
#' the file (BED, BIM, and FAM) aree creted.
#' 
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{0L} when successful.
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' pathOut <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathOut)){
#'     dir.create(pathOut)
#' }
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#' RAIDS:::generateProfilePlinkFiles(pathOut, pRAIDS) 
#'
#' dir(file.path(pathOut,"profileAdmixture"))
#'
#' o <- lapply(dir(file.path(pathOut,"profileAdmixture")),
#'     FUN=function(x,pathC){unlink(file.path(pathC,x))}, 
#'     pathC = file.path(pathOut,"profileAdmixture"))
#'
#' unlink(file.path(pathOut, "profileAdmixture"), recursive = TRUE)
#' unlink(file.path(pathOut), recursive = TRUE)
#' 
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @keywords internal
generateProfilePlinkFiles <- function(pathOut, pRAIDS) {
    
    
    pathBaseMatP <- file.path(pathOut,"matP")
    pathOutP <- file.path(pathOut,  "profileAdmixture")
    if(!dir.exists(file.path(pathOutP))){
        dir.create(file.path(pathOutP))
    }

    
    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))
    
    admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
    pos <- admixDesc[admixDesc$id == anchorId(pRAIDS), "id.seq"]
    popSelected <- read.gdsn(index.gdsn(gdsReference, "admixture.ag"))
    popSelected <- popSelected[popSelected$id.seq == pos,]
    
    snpgdsClose(gdsReference)


    file.copy(file.path(pathBaseMatP, paste0( pedStudy(pRAIDS)$Name.ID[1], ".P.in")),
          file.path(pathOutP,
                    paste0(pedStudy(pRAIDS)$Name.ID[1], ".", nrow(popSelected), ".P.in")))

    writeBimPruned(fileOut=file.path(pathOutP, paste0(pedStudy(pRAIDS)$Name.ID[1], ".bim")),pRAIDS=pRAIDS)
    writeBedProfile(fileOut=file.path(pathOutP, paste0(pedStudy(pRAIDS)$Name.ID[1], ".bed")), listProfile=pedStudy(pRAIDS)$Name.ID[1], profileOnly = TRUE, pRAIDS=pRAIDS)
    writeFamProfile(fileOut=file.path(pathOutP, paste0(pedStudy(pRAIDS)$Name.ID[1], ".fam")), listProfile=pedStudy(pRAIDS)$Name.ID[1], profileOnly = TRUE, pRAIDS=pRAIDS)
    
    return(0L)
}

#' @title Compute JSD from the synthetic dataframe
#' 
#' @description TODO
#'
#' 
#' @param syntheticQ a \code{data.frame} an with the adixmture proportions 
#' results for the synthetic data.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#' 
#' @return \code{data.frame} with three more column for JSD.
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
syntheticJSD <- function(syntheticQ, pRAIDS) {
    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))
    
    admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
    pos <- admixDesc[admixDesc$id == anchorId(pRAIDS), "id.seq"]
    popSelected <- read.gdsn(index.gdsn(gdsReference, "admixture.ag"))
    popSelected <- popSelected[popSelected$id.seq == pos,]
    
    snpgdsClose(gdsReference)
    if(! nrow(syntheticQ) > 1){
        stop("Can't run syntheticJSD with only one synthetic data.")
    }
    matR <- rbind(as.matrix(syntheticQ[,popSelected$ag]),
                                            as.matrix(syntheticQ[,paste0(popSelected$ag, ".p")]))
    
    jsd.m <- JSD(as.matrix(matR))
    syntheticQ$JSDSynVsPos <- diag(jsd.m[1:nrow(syntheticQ), (nrow(syntheticQ) + 1):ncol(jsd.m)])
    matR <- rbind(as.matrix(syntheticQ[,paste0(popSelected$ag, ".p")]),
                                            as.matrix(syntheticQ[,paste0(popSelected$ag, ".gs")]))
    
    jsd.m <- JSD(as.matrix(matR))
    syntheticQ$JSDPosVsGs <- diag(jsd.m[1:nrow(syntheticQ), (nrow(syntheticQ) + 1):ncol(jsd.m)])
    matR <- rbind(as.matrix(syntheticQ[,popSelected$ag]),
                                            as.matrix(syntheticQ[,paste0(popSelected$ag, ".gs")]))
    
    jsd.m <- JSD(as.matrix(matR))
    syntheticQ$JSDSynVsGs <- diag(jsd.m[1:nrow(syntheticQ), (nrow(syntheticQ) + 1):ncol(jsd.m)])
    return(syntheticQ)
}

#' @title Add a column snvKeep to the \code{data.frame} keep 
#' in the specificSNV(pRAIDS)
#' 
#' @description TODO
#'
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters 
#' 
#' @return The function returns a \code{RAIDSparam} an object with all the RAIDS
#' parameters updated.
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#' fileGDS <- file.path(dataDir, "ex1_good_small_1KG.gds")
#' pRAIDS <- RAIDSparam( fileReferenceGDS=fileGDS)
#' 
#' gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))
#' snp.chromosome <- read.gdsn(node=index.gdsn(gdsReference, "snp.chromosome"))
#' snp.position <- read.gdsn(node=index.gdsn(gdsReference, "snp.position"))
#' snpgdsClose(gdsReference)
#' set.seed(4235)
#' tmp <- sample(seq_len(length(snp.position )),100)
#' 
#' df <- data.frame(snp.chromosome=c(snp.chromosome[tmp], rep(1,10)),
#'                  snp.position=c(snp.position[tmp], 1987486 +10+seq_len(10))  )
#' 
#' specificSNV(pRAIDS) <- df
#' pRAIDS <- RAIDS:::specificSNVColKeep(pRAIDS)
#' head(specificSNV(pRAIDS))
#' tail(specificSNV(pRAIDS))
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @encoding UTF-8
#' @keywords internal
specificSNVColKeep <- function(pRAIDS){
    if(!is.null(specificSNV(pRAIDS))){
        if(! "snvKeep" %in% colnames(specificSNV(pRAIDS))){
            gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))
            snp.chromosome <- read.gdsn(node=index.gdsn(gdsReference, "snp.chromosome"))
            keepPos <- seq_len(length(snp.chromosome))
            if("snp.KeepDefault" %in% ls.gdsn(gdsReference) ){
                keepPos <- read.gdsn(index.gdsn(gdsReference, "snp.KeepDefault"))
                snp.chromosome <- snp.chromosome[keepPos]
            }
            snp.position <- read.gdsn(node=index.gdsn(gdsReference, "snp.position"))[keepPos]
            
            snpgdsClose(gdsReference)
            df <- specificSNV(pRAIDS)
            df$snp.chromosome <- as.integer(gsub("chr", "", df$snp.chromosome))

            df <- df[!duplicated(df[,c("snp.chromosome", "snp.position")]),]
            z <- cbind(c(df$snp.chromosome,
                         snp.chromosome,
                         df$snp.chromosome),
                       c(df$snp.position,
                         snp.position,
                         df$snp.position),
                       c(-1 * seq_len(nrow(df)),
                         rep(0, length(snp.position)),
                         seq_len(nrow(df))),
                       c(rep(0, nrow(df)),
                         seq_len(length(snp.position)),
                         rep(0, nrow(df))))
            z <- z[order(z[,1], z[,2], z[,3]),]
            
            df$snvKeep <- as.integer(rep(-1, nrow(df)))
            df$snvKeep[-1*cumsum(z[,3])[z[,3] == 0]] <- as.integer(z[cumsum(z[,3]) < 0 & z[,3] == 0, 4])
            specificSNV(pRAIDS) <- df
        }
    }

    return(pRAIDS)
}

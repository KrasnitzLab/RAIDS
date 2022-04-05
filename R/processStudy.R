#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened
#'
#' @param method the parameter method in SNPRelate::snpgdsLDpruning
#'
#' @param sampleCurrent A \code{vector} of \code{string} corresponding to
#' the sample.ids
#' use in LDpruning
#'
#' @param listSNP the list of snp.id keep
#'
#' @param slide.max.bp.v TODO
#'
#' @param ld.threshold.v TODO
#'
#' @param np TODO
#'
#' @param verbose.v TODO
#'
#' @param chr TODO
#'
#' @param minAF.SuperPop TODO
#'
#' @param keepGDSpruned TODO
#'
#' @param PATHSAMPLEGDS TODO
#'
#' @param keepFile TODO

#' @param PATHPRUNED TODO
#'
#' @param outPref TODO
#'
#'
#' @return TODO a \code{vector} of \code{string}
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
pruningSample <- function(gds,
                          method="corr",
                          sampleCurrent,
                          listSNP = NULL,
                          slide.max.bp.v = 5e5,
                          ld.threshold.v=sqrt(0.1),
                          np = 1,
                          verbose.v=FALSE,
                          chr = NULL,
                          minAF.SuperPop = NULL,
                          keepGDSpruned = TRUE,
                          PATHSAMPLEGDS = NULL,
                          keepFile = FALSE,
                          PATHPRUNED = ".",
                          outPref = "pruned"
                          ){

    filePruned <- file.path(PATHPRUNED, paste0(outPref, ".rds"))
    fileObj <- file.path(PATHPRUNED, paste0(outPref, ".Obj.rds"))

    snp.id <- read.gdsn(index.gdsn(gds, "snp.id"))

    sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))

    posSample <- which(sample.id == sampleCurrent)
    if(length(posSample) != 1){
        stop(paste0("In pruningSample the sample ", sampleCurrent, " doesn't exists\n"))
    }
    # Get the genotype for sampleCurrent
    g <- read.gdsn(index.gdsn(gds, "genotype"), start = c(1, posSample), count = c(-1,1))

    listGeno <- which(g != 3)
    rm(g)

    listKeepPos <- listGeno

    if(!is.null(chr)){
        snpCHR <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
        listKeepPos <- intersect(which(snpCHR == chr), listKeepPos)
    }

    if(!is.null(minAF.SuperPop)){
        snpAF <- read.gdsn(index.gdsn(gds, "snp.AF"))
        listKeepPos <- intersect(which(snpCHR == chr), listKeepPos)
        snpAF <- read.gdsn(index.gdsn(gds, "snp.EAS_AF"))
        listTMP <- NULL
        for(sp in c("EAS", "EUR", "AFR", "AMR", "SAS")){
            snpAF <- read.gdsn(index.gdsn(gds, paste0("snp.", sp, "_AF") ))
            listTMP <- union(listTMP, which(snpAF >= minAF.SuperPop & snpAF <= 1 - minAF.SuperPop))
        }
        listKeepPos <- intersect(listTMP, listKeepPos)
    }

    if(length(listKeepPos) == 0){
        stop(paste0("In pruningSample the sample ", sampleCurrent, " doesn't snp after filters\n"))
    }
    listKeep <- snp.id[listKeepPos]

    sample.ref <- read.gdsn(index.gdsn(gds, "sample.ref"))
    listSamples <- sample.id[which(sample.ref == 1)]

    snpset <- runLDPruning(gds,
                           method,
                           listSamples=listSamples,
                           listKeep=listKeep,
                           slide.max.bp.v = slide.max.bp.v,
                           ld.threshold.v=ld.threshold.v)

    pruned <- unlist(snpset, use.names=FALSE)
    if(keepFile){
        saveRDS(pruned, filePruned)
        saveRDS(snpset, fileObj)
    }
    if(keepGDSpruned){
        fileGDSSample <- file.path(PATHSAMPLEGDS, paste0(sampleCurrent, ".gds"))
        gdsSample <- openfn.gds(fileGDSSample, readonly = FALSE)
        addGDSStudyPruning(gdsSample, pruned, sampleCurrent)
        closefn.gds(gdsSample)
    }

    return(0L)
}



#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened
#'
#' @param gdsSampleFile the path of an object of class \code{gds} related to
#' the sample
#'
#'
#' @return TODO a \code{vector} of \code{string}
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn objdesp.gdsn
#' @encoding UTF-8
#' @export

add1KG2SampleGDS <- function(gds,
                             gdsSampleFile){

    gdsSample <- openfn.gds(gdsSampleFile, readonly = FALSE)

    snp.id <- read.gdsn(index.gdsn(gds,"snp.id"))
    pruned <- read.gdsn(index.gdsn(gdsSample, "pruned.study"))
    listSNP <- which(snp.id %in% pruned)
    listRef <- which(read.gdsn(index.gdsn(gds, "sample.ref"))==1)
    sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))
    sampleCur <- read.gdsn(index.gdsn(gdsSample, "sampleStudy"))
    posCur <- which(sample.id == sampleCur)


    snp.chromosome <- read.gdsn(index.gdsn(gds,"snp.chromosome"))[listSNP]
    snp.position <-  read.gdsn(index.gdsn(gds,"snp.position"))[listSNP]

    add.gdsn(gdsSample, "sample.id", c(sample.id[listRef], sampleCur) )

    add.gdsn(gdsSample, "snp.id", snp.id[listSNP])
    add.gdsn(gdsSample, "snp.chromosome", snp.chromosome)
    add.gdsn(gdsSample, "snp.position", snp.position)
    add.gdsn(gdsSample, "snp.index", listSNP)


    var.geno <- NULL

    j <- 1
    for(i in listRef){
        g <- read.gdsn(index.gdsn(gds, "genotype"), start=c(1,i), count = c(-1,1))[listSNP]

        if(! ("genotype" %in% ls.gdsn(gdsSample))){
            var.geno <- add.gdsn(gdsSample, "genotype",
                                 valdim=c(length(listSNP),
                                          1),
                                 g,
                                 storage="bit2")

        }else{
            if(is.null(var.geno)){
                var.geno <- index.gdsn(gdsSample, "genotype")
            }
            append.gdsn(var.geno, g)
        }
        if(j %% 5 == 0){
            sync.gds(gdsSample)
        }
        j <- j + 1
    }



    add.gdsn(gdsSample, "SamplePos", objdesp.gdsn(index.gdsn(gdsSample, "genotype"))$dim[2] + 1)

    g <- read.gdsn(index.gdsn(gds, "genotype"), start=c(1, posCur), count = c(-1,1))[listSNP]
    append.gdsn(var.geno, g)

    add.gdsn(gdsSample, "lap",
             rep(0.5, objdesp.gdsn(index.gdsn(gdsSample, "genotype"))$dim[1]),
             storage="packedreal8")


    closefn.gds(gdsSample)


}

#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened
#'
#' @param PATHSAMPLEGDS the path of an object of class \code{gds} related to
#' the sample
#'
#' @param PATHGENO TODO
#'
#' @param fileLSNP TODO
#'
#' @return TODO a \code{vector} of \code{string}
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
addPhase1KG2SampleGDSFromFile <- function(gds, PATHSAMPLEGDS,
                                            PATHGENO, fileLSNP) {

    listSel <- readRDS(fileLSNP)

    listGDSSample <- dir(PATHSAMPLEGDS, pattern = ".+.gds")


    indexAll <- NULL
    for(gdsSampleFile in listGDSSample){
        gdsSample <- openfn.gds(file.path(PATHSAMPLEGDS, gdsSampleFile))

        snp.index <- read.gdsn(index.gdsn(gdsSample,"snp.index"))

        indexAll <- union(indexAll, snp.index)
        closefn.gds(gdsSample)
    }

    gdsSample <- createfn.gds(file.path(PATHSAMPLEGDS, "phase1KG.gds"))
    indexAll <- indexAll[order(indexAll)]
    snp.id <- read.gdsn(index.gdsn(gds,"snp.id"))[indexAll]
    add.gdsn(gdsSample, "snp.id", snp.id)
    add.gdsn(gdsSample, "snp.index", indexAll)
    listRef <- which(read.gdsn(index.gdsn(gds, "sample.ref"))==1)
    listSample <- read.gdsn(index.gdsn(gds, "sample.id"))[listRef]
    listSNP <- readRDS(fileLSNP)
    i<-1
    for(sample1KG in listSample){
        print(paste0("P ", i, " ", Sys.time()))
        i<-i+1
        file1KG <- file.path(PATHGENO, paste0(sample1KG,".csv.bz2"))
        matSample <- read.csv2( file1KG,
                                row.names = NULL)
        matSample <- matSample[listSNP[indexAll],, drop=FALSE]
        matSample <- matrix(as.numeric(unlist(strsplit( matSample[,1], "\\|"))),nrow=2)[1,]
        var.phase <- NULL
        if(! ("phase" %in% ls.gdsn(gdsSample))){
            var.phase <- add.gdsn(gdsSample, "phase",
                                 valdim=c(length(indexAll),
                                          1),
                                 matSample,
                                 storage="bit2")

        }else{
            if(is.null(var.phase)){
                var.phase <- index.gdsn(gdsSample, "phase")
            }
            append.gdsn(var.phase, matSample)
        }
    }

    closefn.gds(gdsSample)
    return(0L)

}

#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param listPCA  a \code{list}  with with two objects
#' pca.unrel -> \code{snpgdsPCAClass}
#' and a snp.load -> \code{snpgdsPCASNPLoading}
#'
#' @param sample.current the sample.id to project in the PCA
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1}.
#'
#' @return a \code{snpgdsPCAClass} object, a \code{list} that contains:
#' \itemize{
#'    \item{sample.id} {the sample ids used in the analysis}
#'    \item{snp.id} {the SNP ids used in the analysis}
#'    \item{eigenvalues} {eigenvalues}
#'    \item{eigenvect} {eigenvactors, “# of samples” x “eigen.cnt”}
#'    \item{TraceXTX} {the trace of the genetic covariance matrix}
#'    \item{Bayesian} {whether use bayerisan normalization}
#'}
#'
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom SNPRelate snpgdsPCASampLoading
#' @encoding UTF-8
#' @export
projectSample2PCA <- function(gds,
                              listPCA,
                              sample.current,
                              np = 1) {

    samplePCA <- snpgdsPCASampLoading(listPCA[["snp.load"]],
                                gdsobj=gds, sample.id=sample.current,
                                num.thread=1, verbose=TRUE)

    return(samplePCA)
}

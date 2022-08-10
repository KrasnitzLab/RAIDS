#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \link[gdsfmt]{gdsn.class} (a GDS node), or
#' \link[gdsfmt]{gds.class} (a GDS file) containing the information about
#' 1000 Genome (1KG GDS file).
#'
#' @param nbSamples a single positive \code{integer} representing the number
#' of samples that will be selected for each super-population present in the
#' 1KG GDS file. If the number of samples in a specific super-population is
#' smaller than the \code{nbSamples}, the number of samples in the
#' super-population will correspond to the size of the super-population.
#'
#' @return a \code{data.frame} containing those columns:
#' \itemize{
#' \item{sample.id} { a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group} { a \code{character} string representing the
#' super-population assigned to the sample. }
#' \item{superPop} { a \code{character} string representing the
#' super-population assigned to the sample.
#' }
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
select1KGPop <- function(gds, nbSamples) {

    ## Validate that nbSamples parameter is a single positive numeric
    if(! (isSingleNumber(nbSamples) && nbSamples > 0)) {
        stop("The \'nbSamples\' parameter must be a single positive integer.")
    }

    ## Select reference samples
    listRef <- read.gdsn(index.gdsn(gds, "sample.ref"))
    listKeep <- which(listRef == 1)
    rm(listRef)

    # Extract information about the selected reference samples
    # Including all the super-population classes represented
    sample.annot <- read.gdsn(index.gdsn(gds, "sample.annot"))[listKeep,]
    sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))[listKeep]
    listPop <- unique(sample.annot$pop.group)
    listSel <- list()

    ## For each super-population class, randomly select a fixed number of
    ## samples
    for(i in seq_len(length(listPop))) {
        listGroup <- which(sample.annot$pop.group == listPop[i])
        tmp <- sample(listGroup, min(nbSamples, length(listGroup)))
        listSel[[i]] <- data.frame(sample.id=sample.id[tmp],
                                    pop.group=sample.annot$pop.group[tmp],
                                    superPop=sample.annot$superPop[tmp],
                                    stringsAsFactors=FALSE)
    }

    df <- do.call(rbind, listSel)

    return(df)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param dataRef a \code{data.frame} containing those columns:
#' \itemize{
#' \item{sample.id} { TODO }
#' \item{pop.group} { TODO }
#' \item{superPop} { TODO }
#' }
#'
#' @return a \code{data.frame} containing TODO
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @encoding UTF-8
#' @export
splitSelectByPop <- function(dataRef){

    tmp <- table(dataRef$pop.group)
    if(length(which(tmp != tmp[1])) != 0){
        stop("splitSelectByPop with dataRef with different number of pop\n")
    }

    listPOP <- unique(dataRef$pop.group)
    sampleRM <- vapply(listPOP, function(x, dataRef){
        return(dataRef[which(dataRef$pop.group == x), "sample.id"])
    }, FUN.VALUE = character(tmp[1]), dataRef = dataRef)

    return(sampleRM)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param gdsSampleFile a \code{character} string representing the file name
#' of the GDS file containing the information about the samples to be opened.
#'
#' @param listSampleRef a \code{character} string TODO
#'
#' @param data.id.profile a \code{character} string TODO
#'
#' @param studyDF TODO
#'
#' @param nbSim a single positive \code{integer} representing the number of
#' simulations per combination of sample and profile. Default: \code{1}.
#'
#' @param prefId a \code{character} string TODO
#'
#' @param pRecomb a single \code{numeric} between 0 and 1 TODO.
#' Default: \code{0.01}.
#'
#' @param minProb a single \code{numeric} TODO. Default: \code{0.999}.
#'
#' @param seqError a single positive \code{numeric} between 0 and 1
#' representing the sequencing error rate. Default: \code{0.001}.
#'
#' @return \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
prepSynthetic <- function(gdsSampleFile,
                            listSampleRef,
                            data.id.profile,
                            studyDF,
                            nbSim = 1,
                            prefId = "",
                            pRecomb = 0.01,
                            minProb=0.999,
                            seqError=0.001) {

    gdsSample <- openfn.gds(gdsSampleFile, readonly = FALSE) #

    study.SRC <- read.gdsn(index.gdsn(gdsSample, "study.annot"))
    posStudy <- which(study.SRC$data.id == data.id.profile)
    if(length(posStudy) != 1){
        closefn.gds(gdsSample)
        stop("Error with the data.id of the profile for synthetic data ",
                data.id.profile, "\n")
    }

    sampleSim <- paste(paste0(prefId, ".",data.id.profile),
                        paste(rep(listSampleRef,each=nbSim),
                                seq_len(nbSim), sep="."), sep = ".")

    if(length(which(sampleSim %in% study.SRC$data.id)) > 0) {
        closefn.gds(gdsSample)
        stop("Error data.id of the simulation exists change prefId\n")
    }

    study.list <- data.frame(study.id=studyDF$study.id,
                     study.desc=studyDF$study.desc,
                     study.platform="Synthetic",
                     stringsAsFactors=FALSE)

    pedSim <- data.frame(Name.ID=sampleSim,
                Case.ID=rep(listSampleRef, each=nbSim),
                Sample.Type=rep("Synthetic", length(listSampleRef) * nbSim),
                Diagnosis=rep(study.SRC$diagnosis[posStudy],
                                length(listSampleRef) * nbSim),
                Source=rep("Synthetic", length(listSampleRef) * nbSim),
                stringsAsFactors=FALSE)

    addStudyGDSSample(gdsSample, pedSim, batch=1, listSamples=NULL, study.list)

    closefn.gds(gdsSample)
    return(0L)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened
#'
#' @param gdsRefAnnot an object of class \code{gds} opened
#'
#' @param gdsSampleFile a \code{character} string representing the file name of
#' the GDS file containing the information about the samples.
#'
#' @param listSampleRef a \code{character} string TODO
#'
#' @param data.id.profile a \code{character} string TODO
#'
#' @param nbSim a single positive \code{integer} TODO . Default: \code{1}.
#'
#' @param prefId a \code{string} TODO . Default: \code{""}.
#'
#' @param pRecomb a \code{numeric} between 0 and 1 TODO. Default: \code{0.01}.
#'
#' @param minProb a single \code{numeric}  TODO. Default: \code{0.999}.
#'
#' @param seqError a single positive \code{numeric} between 0 and 1
#' representing the sequencing error rate. Default: \code{0.001}.
#'
#' @return \code{OL} when the function is successful.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats rmultinom
#' @encoding UTF-8
#' @export
syntheticGeno <- function(gds,
                            gdsRefAnnot,
                            gdsSampleFile,
                            data.id.profile,
                            listSampleRef,
                            nbSim=1,
                            prefId="",
                            pRecomb=0.01,
                            minProb=0.999,
                            seqError=0.001) {


    gdsSample <- openfn.gds(filename=gdsSampleFile, readonly=FALSE)


    sampleSim <- paste(paste0(prefId, ".", data.id.profile),
                        paste(rep(listSampleRef,each=nbSim),
                                seq_len(nbSim), sep="."), sep = ".")

    sample.id <- read.gdsn(index.gdsn(gdsSample, "sample.id"))

    if(length(which(sampleSim %in% sample.id)) > 0){
        closefn.gds(gdsSample)
        stop("Error data.id of the simulation exists change prefId\n")
    }
    sample.1kg <- read.gdsn(index.gdsn(gds, "sample.id"))
    listPosRef <- which(sample.id %in% listSampleRef)
    listPosRef.1kg <- which(sample.1kg %in% listSampleRef)


    superPop <- read.gdsn(index.gdsn(gds, "sample.annot/superPop"))[listPosRef.1kg]
    if(all.equal(sample.id[listPosRef], sample.1kg[listPosRef.1kg]) != TRUE) {
        stop("Order between 1kg and the sample are not the same\n")
    }

    # Get the index of the snp.id from gdsSample in the gds1KG
    list1KG <- read.gdsn(index.gdsn(gdsSample, "snp.index"))



    infoSNV <- data.frame(count.tot=read.gdsn(index.gdsn(gdsSample,
                                            "Total.count"))[list1KG],
                          lap=read.gdsn(index.gdsn(gdsSample, "lap")))

    nbSNV <- nrow(infoSNV)

    # Define a table for each "count.tot","lap" and, Freq (number of occurence)
    # to reduce the numbe of sampling call later
    df <- as.data.frame(table(infoSNV[,c("count.tot","lap")]))
    df <- df[df$Freq > 0, ]
    df <- df[order(df$count.tot, df$lap), ]
    # order of SNV relatively to df
    listOrderSNP <- order(infoSNV$count.tot, infoSNV$lap)

    # pos in listOrderSNP of each entries of df
    # df[i, ] where i < nrow(df) link to
    # infoSNV[listOrderSNP[(posDF[i]+1):(posDF[i+1])],]
    posDF <- c(0,cumsum(df$Freq))

    block.Annot <- read.gdsn(index.gdsn(gdsRefAnnot, "block.annot"))

    posSP <- data.frame(EAS=which(block.Annot$block.id == "EAS.0.05.500k"),
                        EUR=which(block.Annot$block.id == "EUR.0.05.500k"),
                        AFR=which(block.Annot$block.id == "AFR.0.05.500k"),
                        AMR=which(block.Annot$block.id == "AMR.0.05.500k"),
                        SAS=which(block.Annot$block.id == "SAS.0.05.500k"))

    #g <- read.gdsn(index.gdsn(gds, "genotype"), start=c(1,i), count = c(-1,1))[listSNP]
    blockDF <- data.frame(EAS=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                            start=c(1,posSP$EAS), count = c(-1,1))[list1KG],
                          EUR=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                            start=c(1,posSP$EUR), count = c(-1,1))[list1KG],
                          AFR=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                            start=c(1,posSP$AFR), count = c(-1,1))[list1KG],
                          AMR=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                            start=c(1,posSP$AMR), count = c(-1,1))[list1KG],
                          SAS=read.gdsn(index.gdsn(gdsRefAnnot, "block"),
                            start=c(1,posSP$SAS), count = c(-1,1))[list1KG])
    # prec <- -1
    # minCur <- 3

    # Better to have a field with the seg in the GDS
    #infoSNV$zone <- read.gdsn(index.gdsn(gdsSample, "segment"))


    # For each reference simulate
    for(r in seq_len(length(listPosRef))) {

        curSynt <- listPosRef[r]
        r.1kg <- which(sample.id[listPosRef[r]] == sample.1kg)
        # get the genotype of the sample r
        g <- read.gdsn(index.gdsn(gdsSample, "genotype"),
                        start=c(1, curSynt),
                        count=c(-1, 1))

        # Order the SNV by count.tot and, lap
        gOrder <- g[listOrderSNP]


        matSim1 <- matrix(nrow=sum(df$Freq), ncol=nbSim)
        matSim2 <- matrix(nrow=sum(df$Freq), ncol=nbSim)


        # Loop on the read.count and lap
        # Faster to group the read.count and lap
        # to run rmultinom
        for(i in seq_len(nrow(df))){

            lap <- as.numeric(as.character(df$lap[i]))

            # number of SNV heterozygote corresponding to
            # df$count.tot[i] and df$lap[i]
            hetero <- which(gOrder[(posDF[i]+1):(posDF[i+1])] == 1)
            nbHetero <- length(hetero)
            # Define the tree prob for the muultinomial
            p1 <- lap * (1- 3 *seqError) +
                (1 - lap) * seqError

            p2 <- (1 - lap) * (1- 3 *seqError) +
                lap * seqError

            p3 <- 2 * seqError

            tmp <- rmultinom(nbHetero * nbSim,
                             as.numeric(as.character(df$count.tot[i])),
                             c(p1, p2, p3))
            # depht of allele 1
            matSim1[listOrderSNP[hetero + posDF[i]],] <- matrix(tmp[1,],
                                                                    ncol=nbSim)
            # depht of allele 2
            matSim2[listOrderSNP[hetero + posDF[i]],] <- matrix(tmp[2,],
                                                                    ncol=nbSim)

            # number of SNV homozygote corresponding to
            # df$count.tot[i] and df$lap[i]
            nbHomo <- df$Freq[i] - nbHetero
            homo <- which(gOrder[(posDF[i]+1):(posDF[i+1])] != 1)



            tmpHomo <- rmultinom(nbHomo * nbSim,
                                 as.numeric(as.character(df$count.tot[i])),
                                 c(1- 3 * seqError,
                                   seqError,
                                   2*seqError))

            # The order between between ref and alt is done with the phase
            # later
            matSim1[listOrderSNP[homo + posDF[i]],] <- matrix(tmpHomo[1,],
                                                              ncol=nbSim)
            matSim2[listOrderSNP[homo + posDF[i]],] <- matrix(tmpHomo[2,],
                                                              ncol=nbSim)
        }

        # superPop of the 1kg sample r is the same
        # for 1kg in list listPosRef.1kg and listPosRef
        curSP <- superPop[r]
        # define a negative block for SNV not in block
        #blockDF[,curSP][which(blockDF[,curSP] == 0)] <- -1*seq_len(length(which(blockDF[,curSP] == 0)))
        if(length(which(blockDF[,curSP] == 0)) > 0){
            stop("There is block set to 0\n")
        }

        listB <- unique(blockDF[,curSP])

        # block where the phase switch
        recombSwitch <- matrix(sample(x=c(0, 1), size=nbSim *(length(listB)),
                                replace=TRUE,
                                prob=c(1-pRecomb, pRecomb)), ncol=nbSim)

        #rownames(recombSwitch) <- listB

        # indice for each zone with the same phase
        blockZone <- apply(recombSwitch, 2, cumsum)

        rownames(blockZone) <- listB


        # We have to manage multipple simulation which mean
        # different number of zone for the different simulation
        LAPparent <- matrix(nrow = nbSNV, ncol = nbSim)
        for(i in seq_len(nbSim)){
            # list of zone with the same phase relatively to 1KG
            listZone <- unique(blockZone[,i])

            # matrix if the lap is the first entry in the phase or the second for each zone
            lapPos <- matrix(sample(x=c(0,1), size=1 *(length(listZone)),
                                        replace=TRUE), ncol=1)

            rownames(lapPos) <- listZone
            #LAPparent <- matrix(nr=nbSNV, nc=nbSim)
            LAPparent[, i] <- lapPos[as.character(blockZone[as.character(blockDF[,curSP]),i]),]
        }

        phaseVal <- read.gdsn(index.gdsn(gdsRefAnnot, "phase"),
                                start = c(1,listPosRef.1kg[r]),
                                count = c(-1,1))[list1KG]


        # mat1 is lap mat2 is 1-lap
        # LAPparent if 0 lap left and 1 lap is right

        # Ok note phaseVal must be the first allele
        tmp <- phaseVal + g * LAPparent
        refC <- matSim1 * ((tmp+1) %% 2) + matSim2 * ((tmp) %% 2)
        altC <- matSim1 * ((tmp) %% 2) + matSim2 * ((tmp+1) %% 2)
        rm(phaseVal, tmp)

        # infoSNV$count.tot
        listCount <- table(infoSNV$count.tot)
        cutOffA <- data.frame(count = unlist(vapply(as.integer(names(listCount)),
                                                    FUN=function(x, minProb, eProb){return(max(2,qbinom(minProb, x,eProb)))},
                                                    FUN.VALUE = numeric(1), minProb=minProb, eProb= 2 * seqError )),
                              allele = unlist(vapply(as.integer(names(listCount)),
                                                     FUN=function(x, minProb, eProb){return(max(2,qbinom(minProb, x,eProb)))},
                                                     FUN.VALUE = numeric(1), minProb=minProb, eProb=seqError)))
        row.names(cutOffA) <- names(listCount)

        gSyn <- matrix(rep(-1, nbSim * nrow(infoSNV)), nrow = nrow(infoSNV))

        # g <- -1 if infoSNV$count.tot - (refC + altC) >= cutOffA[as.character(infoSNV$count.tot), "count"]
        # g <- 0 if altC == 0 & infoSNV$count.tot - (refC + altC) < cutOffA[as.character(infoSNV$count.tot), "count"]
        # g <- 1 if altR > 0 & infoSNV$count.tot - (refC + altC) < cutOffA[as.character(infoSNV$count.tot), "count"]
        # g <- 2 if altR == 0 & infoSNV$count.tot - (refC + altC) < cutOffA[as.character(infoSNV$count.tot), "count"]

        gSyn <- gSyn +
            (infoSNV$count.tot - (refC + altC) < cutOffA[as.character(infoSNV$count.tot), "count"]) * # multiply by 0 if too much error
            ( (refC == 0 |  altC == 0) + # 1 if homozygot
                  (refC >= cutOffA[as.character(infoSNV$count.tot), "allele"]) *
                  (altC >= cutOffA[as.character(infoSNV$count.tot), "allele"]) # 1 if both allele are higher than cutoff hetero
            ) * # 1 if homozygote or hetero and 0 if both > 0 both can't decide if error or hetero
            (1 + (altC > 0) * (1 + (refC == 0) ) ) # if altC == 0 than 1, altC > 0 and refC == 0 than 3

        appendGDSSampleOnly(gdsSample, paste(paste0(prefId, ".",
                                                    data.id.profile),
                                             paste(rep(sample.id[curSynt],
                                                       each=nbSim),
                                                   seq_len(nbSim), sep="."),
                                             sep = "."))
        appendGDSgenotypeMat(gdsSample, gSyn)
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
#' @param gdsSample TODO
#'
#' @param study.id TODO
#'
#' @param popName TODO
#'
#' @return \code{data.frame} TODO study.annot with study.annot == study.id and
#' with the column popName.
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats rmultinom
#' @encoding UTF-8
#' @export
#'
prepPedSynthetic1KG <- function(gds, gdsSample, study.id, popName){


    study.annot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    studyCur <- study.annot[which(study.annot$study.id == study.id),]
    rm(study.annot)
    dataRef <- read.gdsn(index.gdsn(node=gds, "sample.annot"))
    if(! popName %in% colnames(dataRef)) {
        stop("The population ", popName, " is not supported")
    }
    row.names(dataRef) <- read.gdsn(index.gdsn(node=gds, "sample.id"))

    studyCur[[popName]] <- dataRef[studyCur$case.id, popName]
    rownames(studyCur) <- studyCur$data.id
    return(studyCur)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param matKNN TODO
#'
#' @param pedCall TODO see return of prepPedSynthetic1KG
#'
#' @param refCall TODO column name in pedCall with the call
#'
#' @param predCall TODO column name in matKNN with the call
#'
#' @param listCall TODO array of the possible call
#'
#' @return \code{list} TODO
#' with the column popName.
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @encoding UTF-8
#' @export
#'
computeSyntheticConfMat <- function(matKNN, pedCall, refCall, predCall, listCall){

    matAccuracy <- data.frame(pcaD = matKNN$D[1],
                              K = matKNN$K[1],
                              Accu.CM = numeric(1),
                              CM.CI =  numeric(1),
                              N = nrow(matKNN),
                              NBNA = length(which(is.na(matKNN[[predCall]])) ))
    i <- 1
    if(length(unique(matKNN$D)) != 1 | length(unique(matKNN$K)) != 1){
        stop("Compute synthetic accuracy with different pca dimension or K\n")
    }

    #matCur <- matKNN[which(matKNN$D == pcaD & matKNN$K == k),]
    listKeep <- which(!(is.na(matKNN[[predCall]])) )

    fCall <- factor(pedCall[matKNN$sample.id[listKeep], refCall],
                    levels = listCall,
                    labels = listCall)

    fP <- factor(matKNN[[predCall]][listKeep],
                 levels = listCall,
                 labels = listCall)

    cm <- table(fCall,
                fP)


    matAccuracy[i, 3] <- sum(diag(cm[rownames(cm) %in% listCall,
                                     colnames(cm) %in% listCall])) /
        nrow(pedCall[matKNN$sample.id, ][listKeep,])

    matAccuracy[i, 4] <- 1.96 * (matAccuracy[i, 3] * (1 - matAccuracy[i, 3]) /
                                     nrow(pedCall[matKNN$sample.id, ][listKeep,]))^0.5

    res <- list(confMat = cm,
                matAccuracy = matAccuracy)
    return(res)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param matKNN TODO
#'
#' @param pedCall TODO see return of prepPedSynthetic1KG
#'
#' @param refCall TODO column name in pedCall with the call
#'
#' @param predCall TODO column name in matKNN with the call
#'
#' @param listCall TODO array of the possible call
#'
#' @return \code{list} TODO
#' with the column popName.
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom pROC multiclass.roc roc
#' @encoding UTF-8
#' @export
#'
computeSyntheticROC <- function(matKNN, pedCall, refCall, predCall, listCall){

    matAccuracy <- data.frame(pcaD = matKNN$D[1],
                              K = matKNN$K[1],
                              ROC.AUC = numeric(1),
                              ROC.CI = numeric(1),
                              N = nrow(matKNN),
                              NBNA = length(which(is.na(matKNN[[predCall]])) ) )
    i <- 1
    if(length(unique(matKNN$D)) != 1 | length(unique(matKNN$K)) != 1){
        stop("Compute synthetic accuracy with different pca dimension or K\n")
    }

    #matCur <- matKNN[which(matKNN$D == pcaD & matKNN$K == k),]
    listKeep <- which(!(is.na(matKNN[[predCall]])) )
    #listKeep <- which(!(is.na(pedCall[matKNN$sample.id, refCall])) & fCall %in% listCall)

    fCall <- factor(pedCall[matKNN$sample.id[listKeep], refCall],
                    levels = listCall,
                    labels = listCall)



    predMat <- t(vapply(matKNN[[predCall]][listKeep], FUN=function(x, listCall){
        p <- numeric(length(listCall))
        p[which(listCall == x)] <- 1
        return(p)
    },
    FUN.VALUE=numeric(length(listCall)),
    listCall = listCall))
    colnames(predMat) <- listCall

    listAccuPop <- list()


    df <- data.frame(pcaD = matKNN$D[1],
                     K = matKNN$K[1],
                     Call = listCall,
                     L = NA,
                     AUC = NA,
                     H = NA,
                     stringsAsFactors = FALSE)



    resROC <- multiclass.roc(fCall[listKeep], predMat)
    matAccuracy[i, 3] <- as.numeric(resROC$auc)
    matAccuracy[i, 4] <- 0

    # matAccuracy[i, 6] <- ciBS(fCall[listKeep], predMat, 1,100)
    listROC <- list()
    for(j in seq_len(length(listCall))){
        fCur <- rep(0, length(listKeep))
        fCur[fCall[listKeep] == listCall[j]] <- 1

        if(length(which(fCur == 1))>0){
            listROC[[listCall[j]]] <- roc(fCur ~ predMat[,j], ci=TRUE)
            pos <- which(df$Call == listCall[j])
            for(r in seq_len(3)){
                df[pos, r + 3] <- as.numeric(listROC[[j]]$ci[r])
            }
        }else{
            listROC[[listCall[j]]] <- NA
        }
    }
    res <- list(matAUROC.All = matAccuracy,
                matAUROC.Call = df,
                listROC.Call = listROC)
    return(res)
}


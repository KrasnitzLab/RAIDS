#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \link[gdsfmt]{gdsn.class} (a GDS node), or
#' \link[gdsfmt]{gds.class} (a GDS file) containing the information about
#' 1000 Genome (1kG).
#'
#' @param nbSamples a single positive \code{integer}.
#'
#' @return a \code{data.frame) containing those columns:
#' \itemize{
#' \item{sample.id} { TODO }
#' \item{pop.group} { TODO }
#' \item{superPop} { TODO }
#' }
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleInteger
#' @encoding UTF-8
#' @export
select1KGPop <- function(gds, nbSamples) {

    ## Validate that nbSamples is a single positive integer
    if(! (isSingleInteger(nbSamples) && nbSamples > 0)) {
        stop("The \'nbSamples\' parameter must be a single positive integer.")
    }

    listRef <- read.gdsn(index.gdsn(gds, "sample.ref"))
    listKeep <- which(listRef == 1)
    rm(listRef)

    sample.annot <- read.gdsn(index.gdsn(gds, "sample.annot"))[listKeep,]
    sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))[listKeep]
    listPop <- unique(sample.annot$pop.group)
    listSel <- list()

    for(i in seq_len(length(listPop))) {
        listGroup <- which(sample.annot$pop.group == listPop[i])
        tmp <- sample(listGroup, min(nbSamples, length(listGroup)) )
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
#' @param seqError a single \code{numeric} TODO. Default: \code{0.001}.
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
#' @param nbSim a single positive \code{integer} TODO
#'
#' @param prefId a \code{string} TODO
#'
#' @param pRecomb a \code{numeric} between 0 and 1 TODO. Default: \code{0.01}.
#'
#' @param minProb a single \code{numeric}  TODO. Default: \code{0.999}.
#'
#' @param seqError a single \code{numeric}  TODO. Default: \code{0.001}.
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
                            nbSim = 1,
                            prefId = "",
                            pRecomb = 0.01,
                            minProb=0.999,
                            seqError=0.001) {

    # if(nbSim != 1){
    #     stop("Just 1 simulation is manage yet")
    # }

    gdsSample <- openfn.gds(filename=gdsSampleFile, readonly=FALSE) #


    sampleSim <- paste(paste0(prefId, ".", data.id.profile),
                        paste(rep(listSampleRef,each=nbSim),
                                seq_len(nbSim), sep="."), sep = ".")

    sample.id <- read.gdsn(index.gdsn(gdsSample, "sample.id"))

    if(length(which(sampleSim %in% sample.id)) > 0){
        closefn.gds(gdsSample)
        stop("Error data.id of the simulation exists change prefId\n")
    }

    listPosRef <- which(sample.id %in% listSampleRef)


    superPop <- read.gdsn(index.gdsn(gds, "sample.annot/superPop"))[listPosRef]

    # Get the index of the snp.id from gdsSample in the gds1KG
    list1KG <- read.gdsn(index.gdsn(gdsSample, "snp.index"))



    infoSNV <- data.frame(count.tot = read.gdsn(index.gdsn(gdsSample,
                                            "Total.count"))[list1KG],
                          lap = read.gdsn(index.gdsn(gdsSample, "lap")))

    nbSNV <- nrow(infoSNV)

    # Define a table for each "count.tot","lap" and, Freq (number of occurence)
    df <- as.data.frame(table(infoSNV[,c("count.tot","lap")]))
    df <- df[df$Freq > 0,]
    df <- df[order(df$count.tot, df$lap),]
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
        # get the genotype of the sample r
        g <- read.gdsn(index.gdsn(gdsSample, "genotype"),
                       start = c(1,curSynt),
                       count = c(-1,1))

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
            # depht of allele 1 this is the allele homozygote
            matSim1[listOrderSNP[homo + posDF[i]],] <- matrix(tmpHomo[1,],
                                                                    ncol=nbSim)
            # depht of allele 2 (the depth by error of the other allele )
            matSim2[listOrderSNP[homo + posDF[i]],] <- matrix(tmpHomo[2,],
                                                                    ncol=nbSim)
        }

        # superPop of the 1kg sample
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


        # TODO we have to manage multipple simulation which mean
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
                                start = c(1,curSynt),
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




#' @title Random selection of a specific number of reference samples in each
#' subcontinental population present in the 1KG GDS file
#'
#' @description The function randomly selects a fixed number of reference
#' for each subcontinental population present in the 1KG GDS file. When a
#' subcontinental population has less samples than the fixed number, all
#' samples from the subcontinental population are selected.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param nbSamples a single positive \code{integer} representing the number
#' of samples that will be selected for each subcontinental population present
#' in the 1KG GDS file. If the number of samples in a specific subcontinental
#' population is smaller than the \code{nbSamples}, the number of samples
#' selected in this
#' subcontinental population will correspond to the size of this population.
#'
#' @return a \code{data.frame} containing those columns:
#' \itemize{
#' \item{sample.id} { a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group} { a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop} { a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @examples
#'
#' ## The number of samples needed by subcontinental population
#' ## The number is small for demonstration purpose
#' nbSamples <- 5L
#'
#' ## Open 1KG GDS Demo file
#' ## This file only one superpopulation (for demonstration purpose)
#' data.dir <- system.file("extdata", package="RAIDS")
#' gdsFile <- file.path(data.dir, "gds1KG.gds")
#' gdsFileOpen <- openfn.gds(gdsFile)
#'
#' ## Extract a selected number of random samples
#' ## for each subcontinental population
#' ## In the 1KG GDS Demo file, there is one subcontinental population
#' dataR <- select1KGPop(gds=gdsFileOpen, nbSamples=nbSamples)
#'
#' ## Important to close the 1KG GDS Demo file
#' closefn.gds(gdsFileOpen)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
select1KGPop <- function(gds, nbSamples) {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'")
    }

    ## Validate that nbSamples parameter is a single positive numeric
    if(! (isSingleNumber(nbSamples) && nbSamples > 0)) {
        stop("The \'nbSamples\' parameter must be a single positive integer.")
    }

    ## Select reference samples
    listRef <- read.gdsn(index.gdsn(gds, "sample.ref"))
    listKeep <- which(listRef == 1)
    rm(listRef)

    # Extract information about the selected reference samples
    # Including all the subcontinental population classes represented
    sample.annot <- read.gdsn(index.gdsn(gds, "sample.annot"))[listKeep,]
    sample.id <- read.gdsn(index.gdsn(gds, "sample.id"))[listKeep]
    listPop <- unique(sample.annot$pop.group)
    listSel <- list()

    ## For each subcontinental population, randomly select a fixed number of
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


#' @title Group samples per subcontinental population
#'
#' @description The function groups the samples per subcontinental
#' population and generates a matrix containing the sample identifiers and
#' where each column is a subcontinental
#' population.
#'
#' @param dataRef a \code{data.frame} containing those columns:
#' \itemize{
#' \item{sample.id} { a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group} { a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop} { a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @return a \code{matrix} containing the sample identifiers and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @examples
#'
#' ## A data.frame containing samples from 2 subcontinental populations
#' demo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC", "SampleD"),
#'     pop.group=c("TSI", "TSI", "YRI", "YRI"),
#'     superPop=c("EUR", "EUR", "AFR", "AFR"))
#'
#' ## Generate a matrix populated with the sample identifiers and where
#' ## each row is a subcontinental population
#' splitSelectByPop(dataRef=demo)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
splitSelectByPop <- function(dataRef) {

    ## The dataRef must be an data.frame object
    if (!is.data.frame(dataRef)) {
        stop("The \'dataRef\' must be a data.frame object.")
    }

    ## The dataRef must have a pop.group column
    if (!("pop.group" %in% colnames(dataRef))) {
        stop("The \'dataRef\' must have a column named \'pop.group\'.")
    }

    ## The dataRef must have a sample.id column
    if (!("sample.id" %in% colnames(dataRef))) {
        stop("The \'dataRef\' must have a column named \'sample.id\'.")
    }

    tmp <- table(dataRef$pop.group)

    if(length(which(tmp != tmp[1])) != 0) {
        stop("The number of samples in each subcontinental population ",
                "has to be the same.\n")
    }

    listPOP <- unique(dataRef$pop.group)

    ## Generate a matrix where each column is a subcontinental population
    sampleRM <- vapply(listPOP, function(x, dataRef){
        return(dataRef[which(dataRef$pop.group == x), "sample.id"])
    }, FUN.VALUE = character(tmp[1]), dataRef = dataRef)

    return(sampleRM)
}


#' @title Add information related to the synthetic samples
#' (study and sample information) into a GDS Sample file
#'
#' @description This function add entries related to synthetic samples into
#' a GDS Sample file. The entries are 1) a specific synthetic study and
#' 2) specific synthetic sample information into a GDS Sample file.
#' The study information is appended to the GDS Sample file "study.list" node.
#' In this case, the "study.platform" entry, in the study information node,
#' is always "Synthetic".
#' The sample information for all selected samples is appended to the GDS
#' Sample file "study.annot" node. In this case, the "Source" and the
#' "Sample.Type" entries are always "Synthetic".
#' The synthetic samples are assigned unique names by combining :
#' "prefId", ".", "data.id.profile", ".", "listSampleRef", ".",
#' simulation number (between 1 and nbSim)
#'
#' @param gdsSampleFile a \code{character} string representing the file name
#' of the GDS Sample file containing the information about the sample
#' used to generate the synthetic profiles.
#'
#' @param listSampleRef a \code{vector} of \code{character} string
#' representing the
#' identifiers of the selected 1KG samples that will be used as reference to
#' generate the synthetic profiles.
#'
#' @param data.id.profile a \code{character} string TODO
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 2 columns: "study.id" and "study.desc". Those 2 columns
#' must be in \code{character} strings (no factor). Other columns can be
#' present, such as "study.platform", but won't be used.
#'
#' @param nbSim a single positive \code{integer} representing the number of
#' simulations per combination of sample and 1KG reference. Default: \code{1L}.
#'
#' @param prefId a single \code{character} string representing the prefix that
#' is going to be added to the name of the synthetic profile. The prefix
#' enables the creation of multiple synthetic profile using the same
#' combination of sample and 1KG reference. Default: \code{""}.
#'
#' @return \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
prepSynthetic <- function(gdsSampleFile,
                            listSampleRef,
                            data.id.profile,
                            studyDF,
                            nbSim=1L,
                            prefId="") {

    ## The gdsSampleFile must be a character string and the file must exists
    if (!(is.character(gdsSampleFile) && (file.exists(gdsSampleFile)))) {
        stop("The \'gdsSampleFile\' must be a character string representing ",
                "the GDS Sample information file. The file must exist.")
    }

    ## The listSampleRef must be character string
    if (!is.character(listSampleRef)) {
        stop("The \'listSampleRef\' must be a vector of character strings.")
    }

    ## The study.id must have the 2 mandatory columns
    if(sum(c("study.id", "study.desc") %in% colnames(studyDF)) != 2 ) {
        stop("The \'studyDF\' data frame is incomplete. ",
             "One or more mandatory column is missing.\n")
    }

    ## The nbSim must be a single positive numeric
    if (!(isSingleNumber(nbSim) && nbSim > 0)) {
        stop("The \'nbSim\' must be a single positive integer.")
    }

    ## The prefId must be a single character String
    if (!(is.character(prefId) && length(prefId) == 1)) {
        stop("The \'prefId\' must be a single character string.")
    }

    ## Open the GDS Sample file
    gdsSample <- openfn.gds(gdsSampleFile, readonly=FALSE)

    ## Extract information about the samples listed in the GDS Samples
    study.SRC <- read.gdsn(index.gdsn(gdsSample, "study.annot"))
    posStudy <- which(study.SRC$data.id == data.id.profile)
    if(length(posStudy) != 1) {
        closefn.gds(gdsSample)
        stop("Error with the data.id of the profile for synthetic data ",
                data.id.profile, "\n")
    }

    ## Assign unique names to synthetic samples using
    ## the same of the sample, the name of the 1KG reference sample,
    ## the number of simulations and the prefix
    sampleSim <- paste(paste0(prefId, ".", data.id.profile),
                        paste(rep(listSampleRef,each=nbSim),
                                seq_len(nbSim), sep="."), sep = ".")

    if(length(which(sampleSim %in% study.SRC$data.id)) > 0) {
        closefn.gds(gdsSample)
        stop("Error data.id of the simulation exists change prefId\n")
    }

    ## Create a study information data frame using the information passed
    ## as parameters
    ## The synthetic samples will be associated to this information
    ## The study platform is always set to "Synthetic"
    study.list <- data.frame(study.id=studyDF$study.id,
                     study.desc=studyDF$study.desc,
                     study.platform="Synthetic",
                     stringsAsFactors=FALSE)

    ## Create Pedigree information data frame for the synthetic samples
    ## The sample type is always set to "Synthetic" (idem for the source)
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

#' @title Generate synthetic profiles for each cancer sample and 1KG
#' reference sample combination and add them to the GDS Sample file
#'
#' @description The functions uses one cancer sample in combination with one
#' 1KG reference sample to generate an synthetic profile that is saved in the
#' GDS Sample file. When more than one 1KG reference samples are specified,
#' the function recursively generates synthetic profiles for
#' each cancer sample + 1KG reference sample combination. The number of
#' synthetic profiles generated by combination if specified by the number of
#' simulation requested.
#'
#' @param gds an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file.
#'
#' @param gdsSampleFile a \code{character} string representing the file name of
#' the GDS Sample file containing the information about the sample.
#' The file must exist.
#'
#' @param listSampleRef a \code{vector} of \code{character} strings
#' representing the sample identifiers of the 1KG selected reference samples.
#'
#' @param data.id.profile a \code{character} string representing the unique
#' identifier of the cancer sample.
#'
#' @param nbSim a single positive \code{integer} representing the number of
#' simulations that will be generated per sample + 1KG reference combination.
#' Default: \code{1L}.
#'
#' @param prefId a \code{character} string that represent the prefix that will
#' be added to the name of the synthetic profiles generated by the function.
#' Default: \code{""}.
#'
#' @param pRecomb a single positive \code{numeric} between 0 and 1 that
#' represents the frequency of phase switching in the synthetic profiles,
#' Default: \code{0.01}.
#'
#' @param minProb a single positive \code{numeric} between 0 and 1 that
#' represents the probability that the genotype is correct. TODO.
#' Default: \code{0.999}.
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
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats rmultinom
#' @encoding UTF-8
#' @export
syntheticGeno <- function(gds, gdsRefAnnot,
                            gdsSampleFile,
                            data.id.profile,
                            listSampleRef,
                            nbSim=1L,
                            prefId="",
                            pRecomb=0.01,
                            minProb=0.999,
                            seqError=0.001) {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'.")
    }

    ## The gdsRefAnnot must be an object of class "gds.class"
    if (!inherits(gdsRefAnnot, "gds.class")) {
        stop("The \'gdsRefAnnot\' must be an object of class \'gds.class\'.")
    }

    ## The gdsSampleFile must be an character string and the file must exist
    if (!(is.character(gdsSampleFile) && file.exists(gdsSampleFile))) {
        stop("The \'gdsSampleFile\' must be a character string and the file ",
                "must exist.")
    }

    ## The data.id.profile must be a character string
    if (!(is.character(data.id.profile) && length(data.id.profile) == 1)) {
        stop("The \'data.id.profile\' must be a character string.")
    }

    ## The listSampleRef must be a character string
    if(!is.character(listSampleRef)) {
        stop("The \'listSampleRef\' must be a vector of character strings.")
    }

    ## The parameter nbSim must be a single positive integer
    if(!(isSingleNumber(nbSim) && (nbSim >= 0))) {
        stop("The \'nbSim\' parameter must be a single positive ",
             "numeric value.")
    }

    ## The parameter prefId must be a single character string
    if(!(is.character(prefId) && (length(prefId) == 1))) {
        stop("The \'prefId\' parameter must be a single character ",
             "string.")
    }

    ## The parameter pRecomb must be a single positive integer
    if(!(isSingleNumber(pRecomb) && (pRecomb >= 0.0) && (pRecomb <= 1.0))) {
        stop("The \'pRecomb\' parameter must be a single positive ",
             "numeric value between 0 and 1.")
    }

    ## The parameter minProb must be a single positive integer
    if(!(isSingleNumber(minProb) && (minProb >= 0.0) && (minProb <= 1.0))) {
        stop("The \'minProb\' parameter must be a single positive ",
             "numeric value between 0 and 1.")
    }

    ## The parameter seqError must be a single positive integer
    if(!(isSingleNumber(seqError) && (seqError >= 0.0) && (seqError <= 1.0))) {
        stop("The \'seqError\' parameter must be a single positive ",
             "numeric value between 0 and 1.")
    }

    ## Open the GDS Sample file
    gdsSample <- openfn.gds(filename=gdsSampleFile, readonly=FALSE)


    ## The name of the simulated profiles
    sampleSim <- paste(paste0(prefId, ".", data.id.profile),
                        paste(rep(listSampleRef,each=nbSim),
                                seq_len(nbSim), sep="."), sep = ".")

    sample.id <- read.gdsn(index.gdsn(gdsSample, "sample.id"))

    if(length(which(sampleSim %in% sample.id)) > 0) {
        closefn.gds(gdsSample)
        stop("Error data.id of the simulation exists change prefId\n")
    }

    ## Find information about the 1KG reference samples used to generate
    ## the simulated profiles
    sample.1kg <- read.gdsn(index.gdsn(gds, "sample.id"))
    listPosRef <- which(sample.id %in% listSampleRef)
    listPosRef.1kg <- which(sample.1kg %in% listSampleRef)


    superPop <- read.gdsn(index.gdsn(gds,
                                "sample.annot/superPop"))[listPosRef.1kg]

    if(! all.equal(sample.id[listPosRef], sample.1kg[listPosRef.1kg])) {
        stop("The order between 1kg and the list of samples is not ",
                "the same.\n")
    }

    ## Get indexes of the SNV associated to the sample from the GDS Sample file
    list1KG <- read.gdsn(index.gdsn(gdsSample, "snp.index"))

    ## Create a table with the coverage and low allelic fraction information
    infoSNV <- data.frame(count.tot=read.gdsn(index.gdsn(gdsSample,
                                            "Total.count"))[list1KG],
                            lap=read.gdsn(index.gdsn(gdsSample, "lap")))

    nbSNV <- nrow(infoSNV)

    ## Define a table for each "count.tot","lap" and Freq (number of occurrence)
    ## to reduce the number of sampling calls later
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

    #g <- read.gdsn(index.gdsn(gds, "genotype"), start=c(1,i),
    #          count = c(-1,1))[listSNP]

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

        # Order the SNV by count.tot and, lap (low allelic proportion)
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
        # blockDF[,curSP][which(blockDF[,curSP] == 0)] <-
        #     -1*seq_len(length(which(blockDF[,curSP] == 0)))
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

            ## matrix if the lap is the first entry in the phase or
            ## the second for each zone
            lapPos <- matrix(sample(x=c(0,1), size=1 *(length(listZone)),
                                        replace=TRUE), ncol=1)

            rownames(lapPos) <- listZone
            #LAPparent <- matrix(nr=nbSNV, nc=nbSim)
            LAPparent[, i] <-
                        lapPos[as.character(blockZone[as.character(blockDF[,
                                                                curSP]),i]),]
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
        cutOffA <- data.frame(count=unlist(vapply(as.integer(names(listCount)),
                FUN=function(x, minProb, eProb){
                        return(max(2,qbinom(minProb, x,eProb))) },
                FUN.VALUE=numeric(1), minProb=minProb,
                                eProb= 2 * seqError )),
                allele=unlist(vapply(as.integer(names(listCount)),
                FUN=function(x, minProb, eProb){
                        return(max(2,qbinom(minProb, x,eProb))) },
                FUN.VALUE=numeric(1), minProb=minProb, eProb=seqError)))
        row.names(cutOffA) <- names(listCount)

        gSyn <- matrix(rep(-1, nbSim * nrow(infoSNV)), nrow=nrow(infoSNV))

        # g <- -1 if infoSNV$count.tot - (refC + altC) >=
        #                   cutOffA[as.character(infoSNV$count.tot), "count"]
        # g <- 0 if altC == 0 & infoSNV$count.tot - (refC + altC) <
        #                   cutOffA[as.character(infoSNV$count.tot), "count"]
        # g <- 1 if altR > 0 & infoSNV$count.tot - (refC + altC) <
        #                   cutOffA[as.character(infoSNV$count.tot), "count"]
        # g <- 2 if altR == 0 & infoSNV$count.tot - (refC + altC) <
        #                   cutOffA[as.character(infoSNV$count.tot), "count"]

        gSyn <- gSyn +
            (infoSNV$count.tot - (refC + altC) <
                cutOffA[as.character(infoSNV$count.tot), "count"]) * # multiply by 0 if too much error
            ((refC == 0 |  altC == 0) + # 1 if homozygot
                  (refC >= cutOffA[as.character(infoSNV$count.tot), "allele"]) *
                  (altC >= cutOffA[as.character(infoSNV$count.tot), "allele"]) # 1 if both allele are higher than cutoff hetero
            ) * # 1 if homozygote or hetero and 0 if both > 0 both can't decide if error or hetero
            (1 + (altC > 0) * (1 + (refC == 0) ) ) # if altC == 0 than 1, altC > 0 and refC == 0 than 3

        ## Append the name of the samples to the GDS Sample file "sample.id"
        ## node
        appendGDSSampleOnly(gds=gdsSample,
            listSamples=paste(paste0(prefId, ".", data.id.profile),
                paste(rep(sample.id[curSynt], each=nbSim),
                seq_len(nbSim), sep="."), sep = "."))

        ## Append the genotype matrix to the GDS Sample file "genotype" node
        appendGDSgenotypeMat(gdsSample, gSyn)
    }

    ## Close GDS Sample file
    closefn.gds(gdsSample)

    return(0L)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class
#' \code{\link[gdsfmt:gds.class]{gdsfmt::gds.class}}, the opened 1 KG GDS file.
#'
#' @param gdsSample an object of class
#' \code{\link[gdsfmt:gds.class]{gdsfmt::gds.class}}, the opened GDS Sample
#' file.
#'
#' @param study.id a \code{character} string
#'
#' @param popName a \code{character} string representing the name of the
#' column from the \code{data.frame} stored in the "sample.annot" node of the
#' 1KG GDS file. The column must be present in the \code{data.frame}.
#'
#'
#' @return \code{data.frame} TODO study.annot with study.annot == study.id and
#' with the column popName.
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats rmultinom
#' @encoding UTF-8
#' @export
prepPedSynthetic1KG <- function(gds, gdsSample, study.id, popName) {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds,  "gds.class")) {
        stop("The \'gds\' must be an object of class \'gds.class\'.")
    }

    ## The gdsSample must be an object of class "gds.class"
    if (!inherits(gdsSample, "gds.class")) {
        stop("The \'gdsSample\' must be an object of class ",
                "\'gds.class\'.")
    }

    ## The study.id must be a character string
    if (!is.character(study.id)) {
        stop("The \'study.id\' parameter must be a character string.")
    }

    ## Extract study information from the GDS Sample file
    study.annot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    ## Retain the information associated to the current study
    studyCur <- study.annot[which(study.annot$study.id == study.id),]
    rm(study.annot)

    ## Get the information from 1KG GDS file
    dataRef <- read.gdsn(index.gdsn(node=gds, "sample.annot"))

    if(! popName %in% colnames(dataRef)) {
        stop("The population ", popName, " is not supported.")
    }

    ## Assign sample names to the information
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
#' @return \code{list} containing 2 entries:
#' \itemize{
#' \item{confMat} { TODO }
#' \item{matAccuracy} { a \code{data.frame} TODO}
#' }
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @encoding UTF-8
#' @export
computeSyntheticConfMat <- function(matKNN, pedCall, refCall,
                                        predCall, listCall) {

    matAccuracy <- data.frame(pcaD=matKNN$D[1],
                        K=matKNN$K[1],
                        Accu.CM=numeric(1),
                        CM.CI=numeric(1),
                        N=nrow(matKNN),
                        NBNA=length(which(is.na(matKNN[[predCall]]))))
    i <- 1
    if(length(unique(matKNN$D)) != 1 | length(unique(matKNN$K)) != 1){
        stop("Compute synthetic accuracy with different pca dimension or K\n")
    }

    #matCur <- matKNN[which(matKNN$D == pcaD & matKNN$K == k),]
    listKeep <- which(!(is.na(matKNN[[predCall]])) )

    fCall <- factor(pedCall[matKNN$sample.id[listKeep], refCall],
                        levels=listCall, labels=listCall)

    fP <- factor(matKNN[[predCall]][listKeep],
                        levels = listCall, labels = listCall)

    cm <- table(fCall, fP)


    matAccuracy[i, 3] <- sum(diag(cm[rownames(cm) %in% listCall,
                                        colnames(cm) %in% listCall])) /
        nrow(pedCall[matKNN$sample.id, ][listKeep,])

    matAccuracy[i, 4] <- 1.96 * (matAccuracy[i, 3] * (1 - matAccuracy[i, 3]) /
                            nrow(pedCall[matKNN$sample.id, ][listKeep,]))^0.5

    ## Generate list that will be returned
    ##
    res <- list(confMat=cm, matAccuracy=matAccuracy)

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
#' @return \code{list} containing 3 entries:
#' \itemize{
#' \item{matAUROC.All}{ TODO }
#' \item{matAUROC.Call}{ TODO }
#' \item{listROC.Call}{a \code{list} TODO}
#' }
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom pROC multiclass.roc roc
#' @encoding UTF-8
#' @export
computeSyntheticROC <- function(matKNN, pedCall, refCall, predCall, listCall) {

    matAccuracy <- data.frame(pcaD=matKNN$D[1],
                             K=matKNN$K[1],
                             ROC.AUC=numeric(1),
                             ROC.CI=numeric(1),
                             N=nrow(matKNN),
                             NBNA=length(which(is.na(matKNN[[predCall]]))))

    i <- 1

    if(length(unique(matKNN$D)) != 1 | length(unique(matKNN$K)) != 1) {
        stop("Compute synthetic accuracy with different pca dimension or K\n")
    }

    #matCur <- matKNN[which(matKNN$D == pcaD & matKNN$K == k),]
    listKeep <- which(!(is.na(matKNN[[predCall]])) )
    #listKeep <- which(!(is.na(pedCall[matKNN$sample.id, refCall])) &
    #               fCall %in% listCall)

    fCall <- factor(pedCall[matKNN$sample.id[listKeep], refCall],
                        levels=listCall, labels=listCall)

    predMat <- t(vapply(matKNN[[predCall]][listKeep], FUN=function(x, listCall){
        p <- numeric(length(listCall))
        p[which(listCall == x)] <- 1
        return(p)
    },
    FUN.VALUE=numeric(length(listCall)),
    listCall = listCall))
    colnames(predMat) <- listCall

    listAccuPop <- list()


    df <- data.frame(pcaD=matKNN$D[1], K=matKNN$K[1], Call=listCall,
                        L=NA, AUC=NA, H=NA, stringsAsFactors=FALSE)

    resROC <- suppressWarnings(multiclass.roc(fCall[listKeep], predMat))
    matAccuracy[i, 3] <- as.numeric(resROC$auc)
    matAccuracy[i, 4] <- 0

    # matAccuracy[i, 6] <- ciBS(fCall[listKeep], predMat, 1,100)
    listROC <- list()
    for(j in seq_len(length(listCall))) {
        fCur <- rep(0, length(listKeep))
        fCur[fCall[listKeep] == listCall[j]] <- 1

        if(length(which(fCur == 1))>0) {
            listROC[[listCall[j]]] <- suppressWarnings(roc(fCur ~ predMat[,j],
                                                            ci=TRUE))
            pos <- which(df$Call == listCall[j])
            for(r in seq_len(3)){
                df[pos, r + 3] <- as.numeric(listROC[[j]]$ci[r])
            }
        }else{
            listROC[[listCall[j]]] <- NA
        }
    }

    res <- list(matAUROC.All=matAccuracy,
                    matAUROC.Call=df,
                    listROC.Call=listROC)

    return(res)
}


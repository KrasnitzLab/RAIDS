#' @title Random selection of a specific number of reference profiles in each
#' subcontinental population present in the 1KG GDS file
#'
#' @description The function randomly selects a fixed number of reference
#' for each subcontinental population present in the 1KG GDS file. When a
#' subcontinental population has less samples than the fixed number, all
#' samples from the subcontinental population are selected.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param nbProfiles a single positive \code{integer} representing the number
#' of samples that will be selected for each subcontinental population present
#' in the 1KG GDS file. If the number of samples in a specific subcontinental
#' population is smaller than the \code{nbProfiles}, the number of samples
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
#' ## Required library
#' library(gdsfmt)
#'
#' ## The number of samples needed by subcontinental population
#' ## The number is small for demonstration purpose
#' nbProfiles <- 5L
#'
#' ## Open 1KG GDS Demo file
#' ## This file only one superpopulation (for demonstration purpose)
#' dataDir <- system.file("extdata", package="RAIDS")
#' fileGDS <- file.path(dataDir, "PopulationReferenceDemo.gds")
#' gdsFileOpen <- openfn.gds(fileGDS, readonly=TRUE)
#'
#' ## Extract a selected number of random samples
#' ## for each subcontinental population
#' ## In the 1KG GDS Demo file, there is one subcontinental population
#' dataR <- select1KGPop(gdsReference=gdsFileOpen, nbProfiles=nbProfiles)
#'
#' ## Close the 1KG GDS Demo file (important)
#' closefn.gds(gdsFileOpen)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @export
select1KGPop <- function(gdsReference, nbProfiles) {

    ## The gdsReference must be an object of class "gds.class"
    if (!inherits(gdsReference, "gds.class")) {
        stop("The \'gdsReference\' must be an object of class \'gds.class\'")
    }

    ## Validate that nbProfiles parameter is a single positive numeric
    if(! (isSingleNumber(nbProfiles) && nbProfiles > 0)) {
        stop("The \'nbProfiles\' parameter must be a single positive integer.")
    }

    ## Select reference samples
    listRef <- read.gdsn(index.gdsn(gdsReference, "sample.ref"))
    listKeep <- which(listRef == 1)
    rm(listRef)

    # Extract information about the selected reference samples
    # Including all the subcontinental population classes represented
    sample.annot <- read.gdsn(index.gdsn(gdsReference,
                                            "sample.annot"))[listKeep,]
    sample.id <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[listKeep]
    listPop <- unique(sample.annot$pop.group)


    ## For each subcontinental population, randomly select a fixed number of
    ## samples
    dfAll <- lapply(seq_len(length(listPop)), function(i) {
            listGroup <- which(sample.annot$pop.group == listPop[i])
            tmp <- sample(listGroup, min(nbProfiles, length(listGroup)))
            return(data.frame(sample.id=sample.id[tmp],
                                pop.group=sample.annot$pop.group[tmp],
                                superPop=sample.annot$superPop[tmp],
                                stringsAsFactors=FALSE)) })

    df <- do.call(rbind, dfAll)

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


#' @title Add information related to the synthetic profiles
#' (study and synthetic reference profiles information) into a Profile GDS file
#'
#' @description This function add entries related to synthetic profiles into
#' a Profile GDS file. The entries are related to two types of information:
#' the synthetic study and the synthetic profiles.
#'
#' The study information is appended to the Profile GDS file "study.list" node.
#' The "study.platform" entry is always set to 'Synthetic'.
#'
#' The profile information, for all selected synthetic profiles, is appended to
#' the Profile GDS file "study.annot" node. Both the "Source" and the
#' "Sample.Type" entries are always set to 'Synthetic'.
#'
#' The synthetic profiles are assigned unique names by combining:
#' \code{prefix}.\code{data.id.profile}.\code{listSampleRef}.\code{simulation number(1 to nbSim)}
#'
#' @param fileProfileGDS a \code{character} string representing the file name
#' of the Profile GDS file containing the information about the reference
#' profiles used to generate the synthetic profiles.
#'
#' @param listSampleRef a \code{vector} of \code{character} string
#' representing the
#' identifiers of the selected 1KG profiles that will be used as reference to
#' generate the synthetic profiles.
#'
#' @param profileID a \code{character} string representing the profile
#' identifier present in the \code{fileProfileGDS} that will be used to
#' generate synthetic profiles.
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
#' @param prefix a single \code{character} string representing the prefix that
#' is going to be added to the name of the synthetic profile. The prefix
#' enables the creation of multiple synthetic profile using the same
#' combination of sample and 1KG reference. Default: \code{""}.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{FALSE}.
#'
#' @return \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Temporary Profile GDS file
#' fileNameGDS <- file.path(tempdir(), "ex1.gds")
#'
#' ## Copy the Profile GDS file demo that has been pruned and annotated
#' file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                  fileNameGDS)
#'
#' ## Information about the synthetic data set
#' syntheticStudyDF <- data.frame(study.id="MYDATA.Synthetic",
#'         study.desc="MYDATA synthetic data", study.platform="PLATFORM",
#'         stringsAsFactors=FALSE)
#'
#' ## Add information related to the synthetic profiles into the Profile GDS
#' prepSynthetic(fileProfileGDS=fileNameGDS,
#'         listSampleRef=c("HG00243", "HG00150"), profileID="ex1",
#'         studyDF=syntheticStudyDF, nbSim=1L, prefix="synthetic",
#'         verbose=FALSE)
#'
#' ## Open Profile GDS file
#' profileGDS <- openfn.gds(fileNameGDS)
#'
#' ## The synthetic profiles should be added in the 'study.annot' entry
#' tail(read.gdsn(index.gdsn(profileGDS, "study.annot")))
#'
#' ## The synthetic study information should be added to
#' ## the 'study.list' entry
#' tail(read.gdsn(index.gdsn(profileGDS, "study.list")))
#'
#' ## Close GDS file (important)
#' closefn.gds(profileGDS)
#'
#' ## Remove Profile GDS file (created for demo purpose)
#' unlink(fileNameGDS, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @export
prepSynthetic <- function(fileProfileGDS, listSampleRef,
                            profileID, studyDF, nbSim=1L,
                            prefix="", verbose=FALSE) {

    validatePepSynthetic(fileProfileGDS=fileProfileGDS,
        listSampleRef=listSampleRef, profileID=profileID,
        studyDF=studyDF, nbSim=nbSim, prefix=prefix, verbose=verbose)

    ## Open the GDS Sample file
    gdsSample <- openfn.gds(fileProfileGDS, readonly=FALSE)

    ## Extract information about the samples listed in the GDS Samples
    studySRC <- read.gdsn(index.gdsn(gdsSample, "study.annot"))
    posStudy <- which(studySRC$data.id == profileID)
    if(length(posStudy) != 1) {
        closefn.gds(gdsSample)
        stop("The profile \'", profileID, "\' is not found in the annoted ",
                    "study present in the Profile GDS file.\n")
    }

    ## Assign unique names to synthetic profiles using
    ## the same name of the profile, the name of the 1KG reference profile,
    ## the number of simulations and the prefix
    sampleSim <- paste(paste0(prefix, ".", profileID),
                        paste(rep(listSampleRef,each=nbSim),
                                seq_len(nbSim), sep="."), sep = ".")

    if(length(which(sampleSim %in% studySRC$data.id)) > 0) {
        closefn.gds(gdsSample)
        stop("Error data.id of the simulation exists change prefix\n")
    }

    ## Create a study information data frame using the information passed
    ## as parameters
    ## The synthetic samples will be associated to this information
    ## The study platform is always set to "Synthetic"
    study <- data.frame(study.id=studyDF$study.id,
                        study.desc=studyDF$study.desc,
                        study.platform="Synthetic", stringsAsFactors=FALSE)

    ## Create Pedigree information data frame for the synthetic samples
    ## The sample type is always set to "Synthetic" (idem for the source)
    pedSim <- data.frame(Name.ID=sampleSim,
                Case.ID=rep(listSampleRef, each=nbSim),
                Sample.Type=rep("Synthetic", length(listSampleRef) * nbSim),
                Diagnosis=rep(studySRC$diagnosis[posStudy],
                                length(listSampleRef) * nbSim),
                Source=rep("Synthetic", length(listSampleRef) * nbSim),
                stringsAsFactors=FALSE)

    ## Add information to the Profile GDS file
    addStudyGDSSample(gdsProfile=gdsSample, pedProfile=pedSim, batch=1,
        listSamples=NULL, studyDF=study, verbose=verbose)

    closefn.gds(gdsSample)

    return(0L)
}

#' @title Generate synthetic profiles for each cancer profile and 1KG
#' reference profile combination and add them to the Profile GDS file
#'
#' @description The functions uses one cancer profile in combination with one
#' 1KG reference profile to generate an synthetic profile that is saved in the
#' Profile GDS file.
#'
#' When more than one 1KG reference profiles are specified,
#' the function recursively generates synthetic profiles for
#' each cancer profile + 1KG reference profile combination.
#'
#' The number of
#' synthetic profiles generated by combination is specified by the number of
#' simulation requested.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened 1KG GDS file.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened 1KG SNV Annotation GDS file.
#'
#' @param fileProfileGDS a \code{character} string representing the file name
#' of Profile GDS file containing the information about the sample.
#' The file must exist.
#'
#' @param listSampleRef a \code{vector} of \code{character} strings
#' representing the sample identifiers of the 1KG selected reference samples.
#'
#' @param profileID a \code{character} string representing the unique
#' identifier of the cancer profile.
#'
#' @param nbSim a single positive \code{integer} representing the number of
#' simulations that will be generated per sample + 1KG reference combination.
#' Default: \code{1L}.
#'
#' @param prefix a \code{character} string that represent the prefix that will
#' be added to the name of the synthetic profiles generated by the function.
#' Default: \code{""}.
#'
#' @param pRecomb a single positive \code{numeric} between 0 and 1 that
#' represents the frequency of phase switching in the synthetic profiles,
#' Default: \code{0.01}.
#'
#' @param minProb a single positive \code{numeric} between 0 and 1 that
#' represents the probability that the genotype is correct.
#' Default: \code{0.999}.
#'
#' @param seqError a single positive \code{numeric} between 0 and 1
#' representing the sequencing error rate. Default: \code{0.001}.
#'
#' @return The integer \code{OL} when the function is successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Profile GDS file (temporary)
#' fileNameGDS <- file.path(tempdir(), "ex1.gds")
#'
#' ## Copy the Profile GDS file demo that has been pruned and annotated
#' file.copy(file.path(dataDir, "ex1_demo_with_pruning_and_1KG_annot.gds"),
#'                  fileNameGDS)
#'
#' ## Information about the synthetic data set
#' syntheticStudyDF <- data.frame(study.id="MYDATA.Synthetic",
#'         study.desc="MYDATA synthetic data", study.platform="PLATFORM",
#'         stringsAsFactors=FALSE)
#'
#' ## Add information related to the synthetic profiles into the Profile GDS
#' prepSynthetic(fileProfileGDS=fileNameGDS,
#'         listSampleRef=c("HG00243", "HG00150"), profileID="ex1",
#'         studyDF=syntheticStudyDF, nbSim=1L, prefix="synthTest",
#'         verbose=FALSE)
#'
#' ## The 1KG files
#' gds1KG <- snpgdsOpen(file.path(dataDir,
#'                             "ex1_good_small_1KG_GDS.gds"))
#' gds1KGAnnot <- openfn.gds(file.path(dataDir,
#'                             "ex1_good_small_1KG_Annot_GDS.gds"))
#'
#' ## Generate the synthetic profiles and add them into the Profile GDS
#' syntheticGeno(gdsReference=gds1KG, gdsRefAnnot=gds1KGAnnot,
#'         fileProfileGDS=fileNameGDS, profileID="ex1",
#'         listSampleRef=c("HG00243", "HG00150"), nbSim=1,
#'         prefix="synthTest",
#'         pRecomb=0.01, minProb=0.999, seqError=0.001)
#'
#' ## Open Profile GDS file
#' profileGDS <- openfn.gds(fileNameGDS)
#'
#' tail(read.gdsn(index.gdsn(profileGDS, "sample.id")))
#'
#' ## Close GDS files (important)
#' closefn.gds(profileGDS)
#' closefn.gds(gds1KG)
#' closefn.gds(gds1KGAnnot)
#'
#' ## Remove Profile GDS file (created for demo purpose)
#' unlink(fileNameGDS, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom stats rmultinom
#' @encoding UTF-8
#' @export
syntheticGeno <- function(gdsReference, gdsRefAnnot, fileProfileGDS, profileID,
    listSampleRef, nbSim=1L, prefix="", pRecomb=0.01, minProb=0.999,
    seqError=0.001) {

    ## Validate the input parameters
    validateSyntheticGeno(gdsReference=gdsReference, gdsRefAnnot=gdsRefAnnot,
                fileProfileGDS=fileProfileGDS, profileID=profileID,
                listSampleRef=listSampleRef, nbSim=nbSim, prefix=prefix,
                pRecomb=pRecomb, minProb=minProb, seqError=seqError)

    ## Open the GDS Sample file
    gdsSample <- openfn.gds(filename=fileProfileGDS, readonly=FALSE)

    ## The name of the simulated profiles
    sampleSim <- paste(paste0(prefix, ".", profileID),
                        paste(rep(listSampleRef,each=nbSim),
                                seq_len(nbSim), sep="."), sep = ".")

    sample.id <- read.gdsn(index.gdsn(gdsSample, "sample.id"))

    if (length(which(sampleSim %in% sample.id)) > 0) {
        closefn.gds(gdsSample)
        stop("Error data.id of the simulation exists change prefix\n")
    }

    ## Find information about the 1KG reference samples used to generate
    ## the simulated profiles
    sample.1kg <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
    listPosRef <- which(sample.id %in% listSampleRef)
    listPosRef.1kg <- which(sample.1kg %in% listSampleRef)

    superPop <- read.gdsn(index.gdsn(gdsReference,
                                "sample.annot/superPop"))[listPosRef.1kg]

    if (! all.equal(sample.id[listPosRef], sample.1kg[listPosRef.1kg])) {
        stop("The order between 1kg and the list of samples is not the same.\n")
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

    # For each reference simulate
    for(r in seq_len(length(listPosRef))) {

        curSynt <- listPosRef[r]
        r.1kg <- which(sample.id[listPosRef[r]] == sample.1kg)
        # get the genotype of the sample r
        g <- read.gdsn(index.gdsn(gdsSample, "genotype"),
                        start=c(1, curSynt), count=c(-1, 1))

        # Order the SNV by count.tot and, lap (low allelic proportion)
        gOrder <- g[listOrderSNP]

        matSim1 <- matrix(nrow=sum(df$Freq), ncol=nbSim)
        matSim2 <- matrix(nrow=sum(df$Freq), ncol=nbSim)

        # Loop on the read.count and lap
        # Faster to group the read.count and lap to run rmultinom
        for(i in seq_len(nrow(df))){

            lap <- as.numeric(as.character(df$lap[i]))

            ## Number of SNV heterozygote corresponding to
            ## df$count.tot[i] and df$lap[i]
            hetero <- which(gOrder[(posDF[i]+1):(posDF[i+1])] == 1)
            nbHetero <- length(hetero)
            # Define the tree prob for the muultinomial
            p1 <- lap * (1- 3 *seqError) + (1 - lap) * seqError

            p2 <- (1 - lap) * (1- 3 *seqError) + lap * seqError

            p3 <- 2 * seqError

            tmp <- rmultinom(nbHetero * nbSim,
                    as.numeric(as.character(df$count.tot[i])),  c(p1, p2, p3))
            # depth of allele 1
            matSim1[listOrderSNP[hetero + posDF[i]],] <- matrix(tmp[1,],
                                                                    ncol=nbSim)
            # depth of allele 2
            matSim2[listOrderSNP[hetero + posDF[i]],] <- matrix(tmp[2,],
                                                                    ncol=nbSim)

            # number of SNV homozygote corresponding to
            # df$count.tot[i] and df$lap[i]
            nbHomo <- df$Freq[i] - nbHetero
            homo <- which(gOrder[(posDF[i]+1):(posDF[i+1])] != 1)

            tmpHomo <- rmultinom(nbHomo * nbSim,
                                    as.numeric(as.character(df$count.tot[i])),
                                    c(1- 3 * seqError, seqError, 2*seqError))

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

        # indice for each zone with the same phase
        blockZone <- apply(recombSwitch, 2, cumsum)

        rownames(blockZone) <- listB


        ## FOR_LOOP modification to be validated by Pascal
        ## Remove commented code and this text after validation

        # We have to manage multiple simulation which means
        # different number of zone for the different simulations
        LAPparent <- matrix(nrow = nbSNV, ncol = nbSim)
        # for(i in seq_len(nbSim)){
        #     # list of zone with the same phase relatively to 1KG
        #     listZone <- unique(blockZone[,i])
        #
        #     ## matrix if the lap is the first entry in the phase or
        #     ## the second for each zone
        #     lapPos <- matrix(sample(x=c(0,1), size=1 *(length(listZone)),
        #                                 replace=TRUE), ncol=1)
        #
        #     rownames(lapPos) <- listZone
        #
        #     LAPparent[, i] <-
        #                 lapPos[as.character(blockZone[as.character(blockDF[,
        #                                                         curSP]),i]),]
        # }


        # We have to manage multiple simulations which means
        # different number of zones for the different simulations
        lapValues <- vapply(seq_len(nbSim), function(i) {
            # list of zone with the same phase relatively to 1KG
            listZone <- unique(blockZone[,i])

            ## matrix if the lap is the first entry in the phase or
            ## the second for each zone
            lapPos <- matrix(sample(x=c(0,1), size=1 *(length(listZone)),
                                    replace=TRUE), ncol=1)

            rownames(lapPos) <- listZone

            return(lapPos[as.character(blockZone[as.character(blockDF[,
                                            curSP]),i]),])
        }, double(nbSNV))
        LAPparent[, seq_len(nbSim)] <- lapValues

        phaseVal <- read.gdsn(index.gdsn(gdsRefAnnot, "phase"),
            start=c(1,listPosRef.1kg[r]), count=c(-1,1))[list1KG]

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

        # count total multiply by 0 if too much error
        gSyn <- gSyn +
            (infoSNV$count.tot - (refC + altC) <
                cutOffA[as.character(infoSNV$count.tot), "count"]) *
            ((refC == 0 |  altC == 0) + # 1 if homozygot
                (refC >= cutOffA[as.character(infoSNV$count.tot), "allele"]) *
                (altC >= cutOffA[as.character(infoSNV$count.tot), "allele"])
                # 1 if both allele are higher than cutoff hetero
            ) * # 1 if homozygote or hetero and 0 if both > 0 both can't
                # decide if error or hetero
            (1 + (altC > 0) * (1 + (refC == 0)))
            # if altC == 0 than 1, altC > 0 and refC == 0 than 3

        ## Append the profile names of to the Profile GDS file "sample.id" node
        appendGDSSampleOnly(gds=gdsSample,
            listSamples=paste(paste0(prefix, ".", profileID),
                paste(rep(sample.id[curSynt], each=nbSim),
                seq_len(nbSim), sep="."), sep = "."))

        ## Append the genotype matrix to the GDS Sample file "genotype" node
        appendGDSgenotypeMat(gdsSample, gSyn)
    }

    ## Close GDS Sample file
    closefn.gds(gdsSample)

    return(0L)
}


#' @title Calculate the AUROC of the inferences for specific
#' values of D and K using the inferred ancestry results from the synthetic
#' profiles.
#'
#' @description The function calculates the AUROC of the inferences for specific
#' values of D and K using the inferred ancestry results from the synthetic
#' profiles. The calculations are done on each super-population separately as
#' well as on all the results together.
#'
#' @param matKNN a \code{data.frame} containing the inferred ancestry results
#' for fixed values of _D_ and _K_. On of the column names of the
#' \code{data.frame} must correspond to the \code{matKNNAncestryColumn}
#' argument.
#'
#' @param matKNNAncestryColumn  a \code{character} string
#' representing the
#' name of the column that contains the inferred ancestry for the specified
#' synthetic profiles. The column must be present in the \code{matKNN}
#' argument.
#'
#' @param pedCall a \code{data.frame} containing the information about
#' the super-population information from the 1KG GDS file
#' for profiles used to generate the synthetic profiles. The \code{data.frame}
#' must contained a column named as the \code{pedCallAncestryColumn} argument.
#' The row names must correspond to the sample identifiers (mandatory).
#'
#' @param pedCallAncestryColumn a \code{character} string representing the
#' name of the column that contains the known ancestry for the reference
#' profiles in the Reference GDS file. The column must be present in
#' the \code{pedCall} argument.
#'
#' @param listCall a \code{vector} of \code{character} strings representing
#' the list of all possible ancestry assignations.
#' Default: \code{c("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @return \code{list} containing 3 entries:
#' \itemize{
#' \item{\code{matAUROC.All}}{ a \code{data.frame} containing the AUROC for all
#' the ancestry results.  }
#' \item{\code{matAUROC.Call}}{ a \code{data.frame} containing the AUROC
#' information for each super-population.}
#' \item{\code{listROC.Call}}{ a \code{list} containing the output from the
#' \code{roc} function for each super-population.}
#' }
#'
#' @examples
#'
#' ## Loading demo dataset containing pedigree information for synthetic
#' ## profiles and known ancestry of the profiles used to generate the
#' ## synthetic profiles
#' data(pedSynthetic)
#'
#' dataDirRes <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' ## The inferred ancestry results for the synthetic data using
#' ## values of D=6 and K=5
#' matKNN <- readRDS(file.path(dataDirRes, "matKNN.RDS"))
#' matKNN <- matKNN[matKNN$K == 6 & matKNN$D == 5, ]
#'
#' ## Compile statistics from the
#' ## synthetic profiles for fixed values of D and K
#' results <- RAIDS:::computeSyntheticROC(matKNN=matKNN,
#'     matKNNAncestryColumn="SuperPop",
#'     pedCall=pedSynthetic, pedCallAncestryColumn="superPop",
#'     listCall=c("EAS", "EUR", "AFR", "AMR", "SAS"))
#'
#' results$matAUROC.All
#' results$matAUROC.Call
#' results$listROC.Call
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom pROC multiclass.roc roc
#' @encoding UTF-8
#' @export
computeSyntheticROC <- function(matKNN, matKNNAncestryColumn, pedCall,
                        pedCallAncestryColumn,
                        listCall=c("EAS", "EUR", "AFR", "AMR", "SAS")) {

    validateComputeSyntheticRoc(matKNN=matKNN,
        matKNNAncestryColumn=matKNNAncestryColumn, pedCall=pedCall,
        pedCallAncestryColumn=pedCallAncestryColumn, listCall=listCall)

    matAccuracy <- data.frame(pcaD=matKNN$D[1], K=matKNN$K[1],
            ROC.AUC=numeric(1), ROC.CI=numeric(1), N=nrow(matKNN),
            NBNA=length(which(is.na(matKNN[[matKNNAncestryColumn]]))))
    i <- 1

    listKeep <- which(!(is.na(matKNN[[matKNNAncestryColumn]])))

    fCall <- factor(pedCall[matKNN$sample.id[listKeep], pedCallAncestryColumn],
                        levels=listCall, labels=listCall)

    predMat <- t(vapply(matKNN[[matKNNAncestryColumn]][listKeep],
        FUN=function(x, listCall){
            p <- numeric(length(listCall))
            p[which(listCall == x)] <- 1
            return(p)
        }, FUN.VALUE=numeric(length(listCall)), listCall = listCall))
    colnames(predMat) <- listCall

    listAccuPop <- list()

    df <- data.frame(pcaD=matKNN$D[1], K=matKNN$K[1], Call=listCall,
                        L=NA, AUC=NA, H=NA, stringsAsFactors=FALSE)

    resROC <- suppressWarnings(multiclass.roc(fCall[listKeep], predMat))
    matAccuracy[i, 3] <- as.numeric(resROC$auc)
    matAccuracy[i, 4] <- 0

    listROC <- list()
    for (j in seq_len(length(listCall))) {
        fCur <- rep(0, length(listKeep))
        fCur[fCall[listKeep] == listCall[j]] <- 1

        if (length(which(fCur == 1))>0) {
            listROC[[listCall[j]]] <- invisible(roc(fCur ~ predMat[,j],
                                                        ci=TRUE, quiet=TRUE))
            pos <- which(df$Call == listCall[j])
            for (r in seq_len(3)) {
                df[pos, r + 3] <- as.numeric(listROC[[j]]$ci[r])
            }
        } else {
            listROC[[listCall[j]]] <- NA
        }
    }

    res <- list(matAUROC.All=matAccuracy,  matAUROC.Call=df,
                        listROC.Call=listROC)
    return(res)
}


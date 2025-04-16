#' @title Validate input parameters for pruningSample() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{pruningSample}} function.
#'
#' @param gdsReference an object of class \link[gdsfmt]{gds.class}
#' (a GDS file), the Population Reference GDS file.
#'
#' @param method a \code{character} string that represents the method that will
#' be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function.
#'
#' @param currentProfile  a \code{character} string
#' corresponding to the profile identifier used in LD pruning done by the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. A Profile GDS file
#' corresponding to the profile identifier must exist and be located in the
#' \code{pathProfileGDS} directory.
#'
#' @param studyID a \code{character} string corresponding to the study
#' identifier used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#' The study identifier must be present in the Profile GDS file.
#'
#' @param listSNP a \code{vector} of SNVs identifiers specifying selected to
#' be passed the the pruning function;
#' if \code{NULL}, all SNVs are used in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#'
#' @param slideWindowMaxBP a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done in the \code{\link[SNPRelate]{snpgdsLDpruning}}
#' function.
#'
#' @param thresholdLD a single \code{numeric} value that represents the LD
#' threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used.
#'
#' @param verbose a \code{logicial} indicating if information is shown
#' during the process in the \code{\link[SNPRelate]{snpgdsLDpruning}}
#' function.
#'
#' @param chr a \code{character} string representing the chromosome where the
#' selected SNVs should belong. Only one chromosome can be handled. If
#' \code{NULL}, the chromosome is not used as a filtering criterion.
#'
#' @param superPopMinAF a single positive \code{numeric} representing the
#' minimum allelic frequency used to select the SNVs. If \code{NULL}, the
#' allelic frequency is not used as a filtering criterion.
#'
#' @param keepPrunedGDS a \code{logicial} indicating if the information about
#' the pruned SNVs should be added to the GDS Sample file.
#'
#' @param pathProfileGDS a \code{character} string representing the directory
#' where the Profile GDS files will be created. The directory must exist.
#'
#' @param keepFile a \code{logical} indicating if RDS files containing the
#' information about the pruned SNVs must be created.
#'
#' @param pathPrunedGDS a \code{character} string representing an existing
#' directory. The directory must exist.
#'
#' @param outPrefix a \code{character} string that represents the prefix of the
#' RDS files that will be generated. The RDS files are only generated when
#' the parameter \code{keepFile}=\code{TRUE}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(dataDir,
#'             "PopulationReferenceDemo.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validatePruningSample(gdsReference=gds1KG, method="corr",
#'      currentProfile="TGCA_01", studyID="TCGA",
#'      listSNP=c("sr10103", "sr10202"), slideWindowMaxBP=1000L,
#'      thresholdLD=0.008, np=1L, verbose=TRUE, chr=1,
#'      superPopMinAF=0.002, keepPrunedGDS=TRUE, pathProfileGDS=dataDir,
#'      keepFile=FALSE, pathPrunedGDS=dataDir, outPrefix="test")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validatePruningSample <- function(gdsReference, method, currentProfile, studyID,
                        listSNP, slideWindowMaxBP, thresholdLD, np, verbose,
                        chr, superPopMinAF, keepPrunedGDS, pathProfileGDS,
                        keepFile, pathPrunedGDS, outPrefix) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gdsReference, name="gdsReference")

    ## The parameter currentProfile must be a character string
    if (!(is.character(currentProfile))) {
        stop("The \'currentProfile\' parameter must be a character string.")
    }

    ## The parameter method must be a character string
    if(!(is.character(method))) {
        stop("The \'method\' parameter must be a character string.")
    }

    ## The parameter thresholdLD must be a single positive integer
    if (!(isSingleNumber(thresholdLD) && (thresholdLD >= 0.0))) {
        stop("The \'thresholdLD\' parameter must be a single positive ",
                "numeric value.")
    }

    ## The parameter slideWindowMaxBP must be a single positive integer
    if (!(isSingleNumber(slideWindowMaxBP) && (slideWindowMaxBP >= 0.0))) {
        stop("The \'slideWindowMaxBP\' parameter must be a single positive ",
                "numeric value.")
    }

    ## The parameter np must be a single positive integer
    if (!(isSingleNumber(np) && (np >= 0.0))) {
        stop("The \'np\' parameter must be a single positive numeric value.")
    }

    ## The parameter keepPrunedGDS must be a logical
    validateLogical(logical=keepPrunedGDS, "keepPrunedGDS")

    ## The parameter verbose must be a logical
    validateLogical(logical=verbose, "verbose")

    ## The parameter pathProfileGDS must be a character string for existing path
    if (!(is.character(pathProfileGDS) && dir.exists(pathProfileGDS))) {
        stop("The \'pathProfileGDS\' parameter must be a character string ",
                "representing an existing directory.")
    }

    ## The parameter keepFile must be a logical
    validateLogical(logical=keepFile, "keepFile")

    ## The parameter pathPrunedGDS must be a character string for existing path
    if (!(is.character(pathPrunedGDS) && dir.exists(pathPrunedGDS))) {
        stop("The \'pathPrunedGDS\' parameter must be a character string ",
                "representing an existing directory.")
    }

    return(0L)
}


#' @title Validate input parameters for computePoolSyntheticAncestryGr()
#' function
#'
#' @description This function validates the input parameters for the
#' \code{\link{computePoolSyntheticAncestryGr}} function.
#'
#' @param gdsProfile an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, the
#' opened Profile GDS file.
#'
#' @param sampleRM a \code{vector} of \code{character} strings representing
#' the identifiers of the population reference samples that should not
#' be used to create the reference PCA.
#'
#' @param spRef a \code{vector} of \code{character} strings representing the
#' known super population ancestry for the population reference profiles.
#' The population reference profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier.
#' The study identifier must be present in the GDS Sample file.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations.
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016).
#'
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned.
#'
#' @param missingRate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA} only
#' with "<= missingRate" only; if \code{NaN}, no missing threshold.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
#'
#' @return The function returns \code{0L} when successful.
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
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The Profile GDS (opened)
#' gdsSample <- openfn.gds(file.path(dataDir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## The known super population ancestry for the population reference profiles
#' spRef <- c("EUR", "SAS", "EAS", "EUR", "AFR")
#' names(spRef) <- c("HG00100", "HG00101", "HG00102", "HG00103", "HG00104")
#'
#' ## The validation should be successful
#' RAIDS:::validateComputePoolSyntheticAncestryGr(gdsProfile=gdsSample,
#'      sampleRM="TGCA_01", spRef=spRef,
#'      studyIDSyn="TCGA", np=1L, listCatPop=c("AFR", "EAS", "SAS"),
#'      fieldPopInfAnc="Pop", kList=seq_len(3),
#'      pcaList=seq_len(10), algorithm="exact", eigenCount=12L,
#'      missingRate=0.02, verbose=FALSE)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateComputePoolSyntheticAncestryGr <- function(gdsProfile, sampleRM,
        spRef, studyIDSyn, np, listCatPop,
        fieldPopInfAnc, kList, pcaList, algorithm, eigenCount, missingRate,
        verbose) {

    ## The gdsProfile must be objects of class "gds.class"
    validateGDSClass(gdsProfile, "gdsProfile")

    ## The parameter sampleRM must be character strings
    if(!(is.character(sampleRM))) {
        stop("The \'sampleRM\' parameter must be a vector of character ",
                "strings.")
    }

    ## the parameter spRef must be a vector of character strings with names
    if(!(is.character(spRef) && !is.null(names(spRef)))) {
        stop("The \'spRef\' parameter must be a vector of character ",
                "strings with profile identifiers as names.")
    }

    ## The parameter studyIDSyn must be a character string
    validateCharacterString(studyIDSyn, "studyIDSyn")

    ## The parameter fieldPopInfAnc must be a character string
    validateCharacterString(fieldPopInfAnc, "fieldPopInfAnc")

    ## The parameter np must be a single positive integer
    if(!(isSingleNumber(np) && (np > 0))) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## The parameter listCatPop must be a character string
    if(!(is.character(listCatPop))) {
        stop("The \'listCatPop\' parameter must be a vector of ",
                "character strings.")
    }

    ## The parameters must be vectors of positive integers
    validatePositiveIntegerVector(kList, "kList")
    validatePositiveIntegerVector(pcaList, "pcaList")

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    if(!(isSingleNumber(eigenCount))) {
        stop("The \'eigenCount\' parameter must be a single integer.")
    }

    ## The parameter missingRate must be a single number [0,1] or NaN
    if(!(((isSingleNumber(missingRate) && missingRate >= 0.0 &&
                    missingRate <= 1.0)) || is.nan(missingRate)))  {
        stop("The \'missingRate\' parameter must be a single positive ",
                "numeric between zero and one or NaN.")
    }

    validateLogical(logical=verbose, name="verbose")
    return(0L)
}


#' @title Validate input parameters for estimateAllelicFraction() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{estimateAllelicFraction}} function.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the Population Reference GDS file.
#'
#' @param gdsProfile an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the Profile GDS file.
#'
#' @param currentProfile a \code{character} string corresponding to
#' the sample identifier as used in \code{\link{pruningSample}} function.
#'
#' @param studyID a \code{character} string corresponding to the name of
#' the study as
#' used in \code{\link{pruningSample}} function.
#'
#' @param chrInfo a \code{vector} of \code{integer} values representing
#' the length of the chromosomes.
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "DNA" and "RNA". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"DNA"}.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage.
#'
#' @param minProb  a single positive \code{numeric} between 0 and 1 that
#' represents the probability that the genotype is correct.
#'
#' @param eProb a single \code{numeric} between 0 and 1 representing the
#' probability of sequencing error.
#'
#' @param cutOffLOH a single \code{numeric} representing the cutoff, in log,
#' for the homozygote score to assign a region as LOH.
#'
#' @param cutOffHomoScore a single \code{numeric} representing the cutoff, in
#' log, that the SNVs in a block are called homozygote by error.
#'
#' @param wAR a single positive \code{integer} representing the size-1 of
#' the window used to compute an empty box.
#'
#' @param cutOffAR a single \code{numeric} representing the cutoff, in
#' log score, that the SNVs in a gene are allelic fraction different from 0.5.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 Population Reference SNV Annotation GDS file.
#' This parameter is RNA specific.
#'
#' @param blockID a \code{character} string corresponding to the block
#' identifier in \code{gdsRefAnnot}. This parameter is RNA specific.
#'
#' @param verbose a \code{logicial} indicating if the function should print
#' message when running.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG Population Reference GDS Demo file (opened)
#' gds1KG <- openfn.gds(file.path(dataDir,
#'                 "PopulationReferenceDemo.gds"), readonly=TRUE)
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(dataDir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#'
#'     ## The validation should be successful
#'     RAIDS:::validateEstimateAllelicFraction(gdsReference=gds1KG,
#'         gdsProfile=gdsSample,
#'         currentProfile="Sample01", studyID="Synthetic", chrInfo=chrInfo,
#'         studyType="DNA", minCov=10L, minProb=0.03, eProb=0.002, cutOffLOH=10,
#'         cutOffHomoScore=11, wAR=2, cutOffAR=10, gdsRefAnnot=gds1KG,
#'         blockID="1", verbose=FALSE)
#'
#'     ## All GDS file must be closed
#'     closefn.gds(gdsfile=gds1KG)
#'     closefn.gds(gdsfile=gdsSample)
#'
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateEstimateAllelicFraction <- function(gdsReference, gdsProfile,
        currentProfile,
        studyID, chrInfo, studyType, minCov, minProb, eProb, cutOffLOH,
        cutOffHomoScore, wAR, cutOffAR, gdsRefAnnot, blockID, verbose) {

    ## The gds and gdsSample must be objects of class "gds.class"
    validateGDSClass(gdsReference, "gdsReference")
    validateGDSClass(gdsProfile, "gdsProfile")

    ## The currentProfile must be a character string
    if (!is.character(currentProfile)) {
        stop("The \'currentProfile\' must be a character string.")
    }

    ## The studyID must be a character string
    if (!is.character(studyID)) {
        stop("The \'studyID\' must be a character string.")
    }

    ## The studyType must be a character string
    if (!is.character(studyType)) {
        stop("The \'studyType\' must be a character string.")
    }

    ## The minCov parameter must be a single positive integer
    if (!(isSingleNumber(minCov) && (minCov >= 0.0))) {
        stop("The \'minCov\' must be a single numeric positive value")
    }

    ## The minProb and eProb must be single positive numeric between 0 and 1
    validateSingleRatio(minProb, "minProb")
    validateSingleRatio(eProb, "eProb")

    ## The wAR parameter must be a single positive numeric superior to 1
    if (!(isSingleNumber(wAR) && (wAR >= 1))) {
        stop("The \'wAR\' parameter must be a single numeric positive value.")
    }

    ## The cutOffLOH parameter must be a single numeric
    if (!(isSingleNumber(cutOffLOH))) {
        stop("The \'cutOffLOH\' parameter must be a single numeric value.")
    }

    ## The cutOffHomoScore parameter must be a single numeric
    if (!(isSingleNumber(cutOffHomoScore))) {
        stop("The \'cutOffHomoScore\' parameter must be a single ",
                "numeric value.")
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, name="verbose")

    return(0L)
}


#' @title Validate input parameters for createStudy2GDS1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{createStudy2GDS1KG}} function.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the VCF output of SNP-pileup for each sample. The
#' SNP-pileup files must be compressed (gz files) and have the name identifiers
#' of the samples. A sample with "Name.ID" identifier would have an
#' associated SNP-pileup file called "Name.ID.txt.gz". The directory must
#' exist.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined.
#'
#' @param fileNameGDS a \code{character} string representing the file name of
#' the Population Reference GDS file. The file must exist.
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param listProfiles a \code{vector} of \code{character} string corresponding
#' to the profile identifiers that will have a GDS Sample file created. The
#' profile identifiers must be present in the "Name.ID" column of the RDS file
#' passed to the \code{filePedRDS} parameter.
#' If \code{NULL}, all profiles in the \code{filePedRDS} are selected.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the Profile GDS files will be created.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo 1KG Population Reference GDS file
#' gds1KG <- file.path(dataDir, "PopulationReferenceDemo.gds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyInfo <- data.frame(study.id="Pancreatic.WES",
#'                 study.desc="Pancreatic study",
#'                 study.platform="WES",
#'                 stringsAsFactors=FALSE)
#'
#' ## PED Study
#' ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
#'             Case.ID=c("TCGA-H01", "TCGA-H02"),
#'             Sample.Type=c("DNA", "DNA"),
#'             Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))
#'
#' ## The validation should be successful
#' RAIDS:::validateCreateStudy2GDS1KG(pathGeno=dataDir, pedStudy=ped,
#'             fileNameGDS=gds1KG, batch=1, studyDF=studyInfo,
#'             listProfiles=c("Sample_01", "Sample_02"),
#'             pathProfileGDS=dataDir,
#'             genoSource="snp-pileup", verbose=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateCreateStudy2GDS1KG <- function(pathGeno, pedStudy, fileNameGDS, batch,
                            studyDF, listProfiles, pathProfileGDS,
                            genoSource, verbose) {

    ## The pathGeno must be a existing directory
    if (!dir.exists(pathGeno)) {
        stop("The \'pathGeno\' must be an existing directory.")
    }

    ## The PED study must have the mandatory columns
    validatePEDStudyParameter(pedStudy=pedStudy)

    ## The fileNameGDS must be a character string and the file must exists
    if (!(is.character(fileNameGDS) && (file.exists(fileNameGDS)))) {
        stop("The \'fileNameGDS\' must be a character string representing ",
                "the Population Reference GDS file. The file must exist.")
    }

    ## The batch must be a single numeric
    if(!(isSingleNumber(batch))) {
        stop("The \'batch\' must be a single integer.")
    }

    ## The Study DF must have the mandatory columns
    validateStudyDataFrameParameter(studyDF=studyDF)

    ## The listProfiles must be a vector of character string
    if (!(is.character(listProfiles) || is.null(listProfiles))) { #
        stop("The \'listProfiles\' must be a vector ",
                "of character strings (1 entry or more) or NULL.")
    }


    ## The pathProfileGDS must be a character string
    if (!is.character(pathProfileGDS)) {
        stop("The \'pathProfileGDS\' must be a character string representing",
                " the path where the Profile GDS files will be generated.")
    }
    if(is.character(listProfiles)){
        for(profileCur in listProfiles){
            if(file.exists(file.path(pathProfileGDS, paste0(profileCur, ".gds")))){
                stop(paste0("The gds file for ", profileCur, " already exist."))
            }
        }
    }
    ## The genoSource must be a character string
    if(!(is.character(genoSource))) {
        stop("The \'genoSource\' parameter must be a character string.")
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, "verbose")

    return(0L)
}

#' @title Validate input parameters for createProfile() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{createStudy2GDS1KG}} function.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined.
#'
#' @param fileNameGDS a \code{character} string representing the file name of
#' the Population Reference GDS file. The file must exist.
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param listProfiles a \code{vector} of \code{character} string corresponding
#' to the profile identifiers that will have a GDS Sample file created. The
#' profile identifiers must be present in the "Name.ID" column of the RDS file
#' passed to the \code{filePedRDS} parameter.
#' If \code{NULL}, all profiles in the \code{filePedRDS} are selected.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the Profile GDS files will be created.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo 1KG Population Reference GDS file
#' gds1KG <- file.path(dataDir, "PopulationReferenceDemo.gds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyInfo <- data.frame(study.id="Pancreatic.WES",
#'                 study.desc="Pancreatic study",
#'                 study.platform="WES",
#'                 stringsAsFactors=FALSE)
#'
#' ## PED Study
#' ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
#'             Case.ID=c("TCGA-H01", "TCGA-H02"),
#'             Sample.Type=c("DNA", "DNA"),
#'             Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))
#'
#' ## The validation should be successful
#' RAIDS:::validateCreateStudy2GDS1KG(pathGeno=dataDir, pedStudy=ped,
#'             fileNameGDS=gds1KG, batch=1, studyDF=studyInfo,
#'             listProfiles=c("Sample_01", "Sample_02"),
#'             pathProfileGDS=dataDir,
#'             genoSource="snp-pileup", verbose=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validatecreateProfile <- function(pedStudy,
                                  fileNameGDS, batch,
                                  studyDF, listProfiles, pathProfileGDS,
                                  genoSource, verbose) {
    # profileFile, profileName,
    # filePedRDS=NULL, pedStudy=NULL, fileNameGDS,
    # batch=1, studyDF, listProfiles=NULL,
    # pathProfileGDS=NULL,
    # genoSource=c("snp-pileup", "generic", "VCF", "bam"),
    # paramProfile=list(ScanBamParam=NULL,
    #                   PileupParam=NULL,
    #                   yieldSize=10000000),
    # verbose=FALSE
    ## The pathGeno must be a existing directory


    ## The PED study must have the mandatory columns
    validatePEDStudyParameter(pedStudy=pedStudy)

    ## The fileNameGDS must be a character string and the file must exists
    if (!(is.character(fileNameGDS) && (file.exists(fileNameGDS)))) {
        stop("The \'fileNameGDS\' must be a character string representing ",
             "the Population Reference GDS file. The file must exist.")
    }

    ## The batch must be a single numeric
    if(!(isSingleNumber(batch))) {
        stop("The \'batch\' must be a single integer.")
    }

    ## The Study DF must have the mandatory columns
    validateStudyDataFrameParameter(studyDF=studyDF)

    ## The listProfiles must be a vector of character string
    if (!(is.character(listProfiles) || is.null(listProfiles))) { #
        stop("The \'listProfiles\' must be a vector ",
             "of character strings (1 entry or more) or NULL.")
    }


    ## The pathProfileGDS must be a character string
    if (!is.character(pathProfileGDS)) {
        stop("The \'pathProfileGDS\' must be a character string representing",
             " the path where the Profile GDS files will be generated.")
    }

    ## The genoSource must be a character string
    if(!(is.character(genoSource))) {
        stop("The \'genoSource\' parameter must be a character string.")
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, "verbose")

    return(0L)
}


#' @title Validate input parameters for computeAncestryFromSyntheticFile()
#' function
#'
#' @description This function validates the input parameters for the
#' \code{\link{computeAncestryFromSyntheticFile}} function.
#'
#' @param gdsReference an object of class \link[gdsfmt]{gds.class} (a GDS
#' file), the opened Population Reference GDS file.
#'
#' @param gdsProfile an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Profile GDS file.
#'
#' @param listFiles a \code{vector} of \code{character} strings representing
#' the name of files that contain the results of ancestry inference done on
#' the synthetic profiles for multiple values of _D_ and _K_. The files must
#' exist.
#'
#' @param currentProfile a \code{character} string representing the profile
#' identifier of the current profile on which ancestry will be inferred.
#'
#' @param spRef a \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The Population
#' Reference profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier. The study identifier must be present in the GDS Sample file.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations.
#'
#' @param fieldPopIn1KG a \code{character} string representing the name of the
#' column that contains the known ancestry for the reference profiles in
#' the Reference GDS file.
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA.
#'
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned.
#'
#' @param missingRate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA}
#' with "<= missingRate" only; if \code{NaN}, no missing threshold.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function.
#'
#' @return The function returns \code{0L} when successful.
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
#' ## Required library
#' library(gdsfmt)
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG Population Reference GDS demo file (opened)
#' gds1KG <- openfn.gds(file.path(dataDir,
#'             "PopulationReferenceDemo.gds"), readonly=TRUE)
#'
#' ## The Profile GDS (opened)
#' gdsSample <- openfn.gds(file.path(dataDir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' listFiles <- file.path(dataDir,  "listSNPIndexes_Demo.rds")
#'
#' ## The validation should be successful
#' RAIDS:::validateComputeAncestryFromSyntheticFile(gdsReference=gds1KG,
#'     gdsProfile=gdsSample, listFiles=listFiles, currentProfile="sample01",
#'     spRef=c("EUR", "AFR"), studyIDSyn="Synthetic", np=1L,
#'     listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
#'     fieldPopIn1KG="superpop", fieldPopInfAnc="Superpop", kList=c(2, 3, 4),
#'     pcaList=c(3, 4, 5), algorithm="exact", eigenCount=32L, missingRate=0.2,
#'     verbose=FALSE)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber isSingleString
#' @encoding UTF-8
#' @keywords internal
validateComputeAncestryFromSyntheticFile <- function(gdsReference, gdsProfile,
                listFiles, currentProfile, spRef, studyIDSyn, np, listCatPop,
                fieldPopIn1KG, fieldPopInfAnc, kList, pcaList,
                algorithm, eigenCount, missingRate, verbose) {

    ## The gdsReference and gdsProfile must be objects of class "gds.class"
    validateGDSClass(gdsReference, "gdsReference")
    validateGDSClass(gdsProfile, "gdsProfile")

    if(!isSingleString(currentProfile)) {
        stop("The \'currentProfile\' parameter must be a character string.")
    }

    if(!(is.character(spRef))) {
        stop("The \'spRef\' parameter must be a vector of character strings.")
    }

    ## The parameter np must be a single positive integer
    if(!(isSingleNumber(np) && (np > 0))) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## The parameters are character strings (vector of 1 entry)
    validateCharacterString(value=fieldPopIn1KG, "fieldPopIn1KG")
    validateCharacterString(value=fieldPopInfAnc, "fieldPopInfAnc")
    validateCharacterString(value=studyIDSyn, "studyIDSyn")

    ## The parameter listCatPop must be a vector of of character strings
    if(!(is.character(listCatPop))) {
        stop("The \'listCatPop\' parameter must be a vector of ",
                "character strings.")
    }

    ## The parameters must be vectors of positive integers
    validatePositiveIntegerVector(kList, "kList")
    validatePositiveIntegerVector(pcaList, "pcaList")

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    if (!(isSingleNumber(eigenCount))) {
        stop("The \'eigenCount\' parameter must be a single integer.")
    }

    ## The missingRate must be a positive numeric between zero and one or NaN
    if (!is.nan(missingRate)) {
        if (!(isSingleNumber(missingRate) && (missingRate >= 0.0) &&
                    (missingRate <= 1.0))) {
            stop("The \'missingRate\' must be a single numeric positive ",
                        "value between 0 and 1 or NaN.")
        }
    }

    validateLogical(verbose, "verbose")

    return(0L)
}


#' @title Validate input parameters for computePCARefSample() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{computePCARefSample}} function.
#'
#' @param gdsProfile an object of class \link[gdsfmt]{gds.class},
#' an opened Profile GDS file.
#'
#' @param currentProfile a single \code{character} string representing the
#' profile identifier.
#'
#' @param studyIDRef a single \code{character} string representing the
#' study identifier.
#'
#' @param np a single positive \code{integer} representing the number of CPU
#' that will be used.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA.
#'
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#'
#' @param missingRate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; if \code{NaN}, no missing
#' threshold.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(dataDir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validateComputePCARefSample(gdsProfile=gdsSample,
#'     currentProfile="HCC01", studyIDRef="1KG", np=1L, algorithm="exact",
#'     eigenCount=32L, missingRate=0.02, verbose=FALSE)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateComputePCARefSample <- function(gdsProfile, currentProfile, studyIDRef,
                                            np, algorithm, eigenCount,
                                            missingRate, verbose) {

    ## The gdsProfile must be object of class "gds.class"
    validateGDSClass(gdsProfile, "gdsProfile")

    ## Validate that currentProfile is a string
    validateCharacterString(currentProfile, "currentProfile")

    ## Validate that studyIDRef is a string
    validateCharacterString(studyIDRef, "studyIDRef")

    ## The parameter np must be a single positive integer
    if(!(isSingleNumber(np) && (np > 0))) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    ## The parameter eigenCount must be a single integer
    if(!(isSingleNumber(eigenCount))) {
        stop("The \'eigenCount\' parameter must be a single integer.")
    }

    ## The missing.rate must be a positive numeric between zero and one or NaN
    if (!is.nan(missingRate)) {
        if (!(isSingleNumber(missingRate) && (missingRate >= 0.0) &&
                (missingRate <= 1.0))) {
            stop("The \'missingRate\' must be a single numeric positive ",
                    "value between 0 and 1 or NaN.")
        }
    }

    validateLogical(verbose, "verbose")

    return(0L)
}


#' @title Validate input parameters for add1KG2SampleGDS() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{add1KG2SampleGDS}} function.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Population Reference
#' GDS file.
#'
#' @param gdsProfileFile a \code{character} string representing the path and
#' file name of the Profile GDS file. The Profile GDS file must exist.
#'
#' @param currentProfile a \code{character} string corresponding to the profile
#' identifier associated to the current list of pruned SNVs.
#'
#' @param studyID a \code{character} string corresponding to the study
#' identifier associated to the current list of pruned SNVs.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG Population Reference GDS demo file (opened)
#' gds1KG <- openfn.gds(file.path(dataDir,
#'                 "PopulationReferenceDemo.gds"), readonly=TRUE)
#'
#' ## The validatiion should be successful
#' RAIDS:::validateAdd1KG2SampleGDS(gdsReference=gds1KG,
#'     gdsProfileFile=file.path(dataDir, "GDS_Sample_with_study_demo.gds"),
#'     currentProfile="Sample01", studyID="Synthetic")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateAdd1KG2SampleGDS <- function(gdsReference, gdsProfileFile,
                                        currentProfile, studyID) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gdsReference, name="gdsReference")

    ## The gdsProfileFile must be a character string and the file must exists
    if(!(is.character(gdsProfileFile) && (file.exists(gdsProfileFile)))) {
        stop("The \'fileProfileGDS\' must be a character string representing ",
                "the GDS Sample file. The file must exist.")
    }

    ## The currentProfile must be a character string
    if(!(is.character(currentProfile))) {
        stop("The \'currentProfile\' must be a character string.")
    }

    ## The studyID must be a character string
    if(!(is.character(studyID))) {
        stop("The \'studyID\' must be a character string.")
    }

    return(0L)
}




#' @title Validate the parameters of the runExomeAncestry() function
#'
#' @description The function validates the input parameters for the
#' \code{\link{runExomeAncestry}} function. When a parameter is not as
#' expected, an error message is generated.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Profile files will be created.
#' Default: \code{NULL}.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the VCF output of SNP-pileup for each sample. The
#' SNP-pileup files must be compressed (gz files) and have the name identifiers
#' of the samples. A sample with "Name.ID" identifier would have an
#' associated SNP-pileup file called "Name.ID.txt.gz". The directory must
#' exist.
#'
#' @param pathOut a \code{character} string representing the path to
#' the directory where the output files are created.
#'
#' @param fileReferenceGDS  a \code{character} string representing the file
#' name of the Population Reference GDS file. The file must exist.
#'
#' @param fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Population Reference GDS annotation file.
#' The file must exist.
#'
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section.
#'
#' @param syntheticRefDF a \code{data.frame} containing those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{FALSE}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Path where the output file will be generated
#' pathOut <- file.path(dataDir, "example", "res.out")
#'
#' ## Study data frame
#' study <- data.frame(study.id = "MYDATA",
#'                       study.desc = "Description",
#'                       study.platform = "PLATFORM",
#'                       stringsAsFactors = FALSE)
#'
#' ## Population Reference GDS demo file
#' gdsRef <- file.path(dataDir, "PopulationReferenceDemo.gds")
#'
#' gdsAnnotRef <- file.path(dataDir, "PopulationReferenceSNVAnnotationDemo.gds")
#'
#' ## Pedigree Study data frame
#' ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
#'             Case.ID=c("TCGA-H01", "TCGA-H02"),
#'             Sample.Type=c("DNA", "DNA"),
#'             Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#'
#'     ## Profiles used for synthetic data set
#'     syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
#'         "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
#'         superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)
#'
#'     ## Returns OL when all parameters are valid
#'     RAIDS:::validateRunExomeOrRNAAncestry(pedStudy=ped, studyDF=study,
#'         pathProfileGDS=dataDir, pathGeno=dataDir, pathOut=pathOut,
#'         fileReferenceGDS=gdsRef, fileReferenceAnnotGDS=gdsAnnotRef,
#'         chrInfo=chrInfo, syntheticRefDF=syntheticRefDF,
#'         genoSource="snp-pileup", verbose=FALSE)
#'
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateRunExomeOrRNAAncestry <- function(pedStudy, studyDF, pathProfileGDS,
    pathGeno, pathOut, fileReferenceGDS, fileReferenceAnnotGDS,
    chrInfo, syntheticRefDF, genoSource, verbose) {

    ## The PED study must have the mandatory columns
    validatePEDStudyParameter(pedStudy=pedStudy)

    ## The study data frame must have the mandatory columns
    validateStudyDataFrameParameter(studyDF=studyDF)

    ## The pathGeno must be a existing directory
    if (!dir.exists(pathGeno)) {
        stop("The \'pathGeno\' must be an existing directory.")
    }

    ## The pathOut must be a character string
    if (!is.character(pathOut)) {
        stop("The \'pathOut\' must be a character string representing",
                " the path where the output files will be generated.")
    }

    ## The fileReferenceGDS must be a character string and the file must exists
    if (!(is.character(fileReferenceGDS) && (file.exists(fileReferenceGDS)))) {
        stop("The \'fileReferenceGDS\' must be a character string ",
                "representing the Reference GDS file. The file must exist.")
    }

    ## The fileReferenceAnnotGDS must be a character string and
    ## the file must exist
    if (!(is.character(fileReferenceAnnotGDS) &&
                (file.exists(fileReferenceAnnotGDS)))) {
        stop("The \'fileReferenceAnnotGDS\' must be a character string ",
                "representing the Reference Annotation GDS file. ",
                "The file must exist.")
    }

    ## The chrInfo must be a vector of integer
    validatePositiveIntegerVector(chrInfo, "chrInfo")

    ## The syntheticRefDF must have the mandatory columns
    validateDataRefSynParameter(syntheticRefDF=syntheticRefDF)

    ## The genoSource must be a character string
    if(!(is.character(genoSource))) {
        stop("The \'genoSource\' parameter must be a character string.")
    }

    validateLogical(verbose, "verbose")

    return(0L)
}


#' @title Validate that the PED study has the mandatory columns
#'
#' @description The function validates the input PED study. The PED study
#' must be a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor).
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor).
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Pedigree Study data frame
#' ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
#'             Case.ID=c("TCGA-H01", "TCGA-H02"),
#'             Sample.Type=c("DNA", "DNA"),
#'             Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))
#'
#' ## Return 0L when PED is valid
#' RAIDS:::validatePEDStudyParameter(pedStudy=ped)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validatePEDStudyParameter <- function(pedStudy) {

    ## The PED study must have the mandatory columns
    if (!(is.data.frame(pedStudy) &&
            all(c("Name.ID", "Case.ID", "Sample.Type", "Diagnosis", "Source")
                            %in% colnames(pedStudy)))) {
        stop("The PED study data frame is incomplete. ",
                        "One or more mandatory columns are missing.")
    }

    return(0L)
}


#' @title Validate that the reference profile data set has
#' the mandatory columns
#'
#' @description The function validates the input reference profile data set.
#' The reference profile data set
#' must be a \code{data.frame} with those mandatory columns:
#' "sample.id", "pop.group", "superPop". All columns must be in
#' \code{character} strings (no factor).
#'
#' @param syntheticRefDF a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The mandatory columns are:
#' "sample.id", "pop.group", "superPop". All columns must be in
#' \code{character} strings (no factor).
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Profiles used for synthetic data set
#' syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
#'     "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
#'     superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)
#'
#' ## Return 0L when the reference profile data set is valid
#' RAIDS:::validateDataRefSynParameter(syntheticRefDF=syntheticRefDF)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateDataRefSynParameter <- function(syntheticRefDF) {

    ## The reference profile data.frame must have the mandatory columns
    if (!(is.data.frame(syntheticRefDF) &&
                all(c("sample.id", "pop.group", "superPop")
                    %in% colnames(syntheticRefDF)))) {
        stop("The reference profile data frame \'syntheticRefDF\' is ",
                "incomplete. One or more mandatory columns are missing. The ",
                "mandatory columns are: \'sample.id\', \'pop.group\', ",
                "\'superPop\'.")
    }

    return(0L)
}


#' @title Validate that the study data set has
#' the mandatory columns
#'
#' @description The function validates the input dtudy data set.
#' The study data set
#' must be a \code{data.frame} with those mandatory columns:
#' "studyID", "study.desc", "study.platform". All columns must be in
#' \code{character} strings (no factor).
#'
#' @param studyDF a \code{data.frame} containing the study information.
#' The mandatory columns are:
#' "study.id", "study.desc", "study.platform". All columns must be in
#' \code{character} strings (no factor).
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Study data frame
#' study <- data.frame(study.id = "MYDATA",
#'                       study.desc = "Description",
#'                       study.platform = "PLATFORM",
#'                       stringsAsFactors = FALSE)
#'
#' ## Return 0L when the study data set is valid
#' RAIDS:::validateStudyDataFrameParameter(studyDF=study)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateStudyDataFrameParameter <- function(studyDF) {

    ## The study data.frame must have the mandatory columns
    if (!(is.data.frame(studyDF) &&
            all(c("study.id", "study.desc", "study.platform")
                                    %in% colnames(studyDF)))) {
        stop("The study data frame \'studyDF\' is incomplete. ",
            "One or more mandatory columns are missing. The mandatory ",
            "columns are: \'study.id\', \'study.desc\', \'study.platform\'.")
    }

    return(0L)
}


#' @title Validate the parameters for the addStudy1Kg() function
#'
#' @description The function validates the input parameters for the
#' \code{\link{addStudy1Kg}} function. When a parameter is not as
#' expected, an error message is generated.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param fileProfileGDS a \code{character} string representing the path and
#' file name of the GDS Sample file. The GDS Sample file must exist.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#' fileReferenceGDS <- file.path(dataDir, "PopulationReferenceDemo.gds")
#' gds1KG <- snpgdsOpen(fileReferenceGDS)
#'
#' ## Path to demo Profile GDS file
#' fileProfileGDS  <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")
#'
#' ## Returns OL when all parameters are valid
#' RAIDS:::validateAddStudy1Kg(gdsReference=gds1KG,
#'     fileProfileGDS=fileProfileGDS, verbose=TRUE)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateAddStudy1Kg <- function(gdsReference, fileProfileGDS, verbose) {

    ## The gdsReference must be an object of class "gds.class"
    validateGDSClass(gds=gdsReference, "gdsReference")

    ## The fileProfileGDS must be a character string and the file must exists
    if(!(is.character(fileProfileGDS) && (file.exists(fileProfileGDS)))) {
        stop("The \'fileProfileGDS\' must be a character string representing ",
                "the GDS Sample file. The file must exist.")
    }

    ## Validate that verbose is a logical
    validateLogical(logical=verbose, name="verbose")

    return(0L)
}



#' @title Validate that the Profile GDS file exists for the specified profile
#'
#' @description The function validates that the Profile GDS file associated
#' to a profile identifier exists in the specified directory.
#'
#' @param pathProfile a \code{character} string representing the directory
#' where the Profile GDS files will be created. The directory must exist.
#'
#' @param profile  a \code{character} string
#' corresponding to the profile identifier. A Profile GDS file
#' corresponding to the profile identifier must exist and be located in the
#' \code{pathProfile} directory.
#'
#' @return a \code{character} string representing the path to the existing
#' Profile GDS file.
#'
#' @examples
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata/tests", package="RAIDS")
#'
#' ## The function returns the path to the existing Profile GDS file
#' RAIDS:::validateProfileGDSExist(pathProfile=dataDir,
#'     profile="ex1_demo")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateProfileGDSExist <- function(pathProfile, profile) {

    pathFile <- file.path(pathProfile, paste0(profile, ".gds"))

    ## The Profile GDS file must exists
    if (!(file.exists(pathFile))) {
        stop("The Profile GDS file \'", pathFile, " does not exist.")
    }

    return(pathFile)
}


#' @title Validate the input parameters for computePCAMultiSynthetic()
#' function
#'
#' @description The function validates the input parameters for the
#' \code{\link{computePCAMultiSynthetic}} function.
#' When a parameter is not as expected, an error message is generated.
#'
#' @param gdsProfile an object of class \link[gdsfmt]{gds.class} (a GDS file),
#' an opened Profile GDS file.
#'
#' @param listPCA a \code{list} containing the PCA \code{object} generated
#' with the 1KG reference profiles (excluding the ones used to generate the
#' synthetic data set) in an entry called \code{"pca.unrel"}.
#'
#' @param sampleRef a \code{vector} of \code{character} strings representing
#' the identifiers of the 1KG reference profiles that should not be used to
#' create the reference PCA.
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier. The study identifier must be present in the Profile GDS file.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function.
#'
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Loading demo PCA on subset of 1KG reference dataset
#' data(demoPCA1KG)
#'
#' ## Path to the demo GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#' fileProfileGDS <- file.path(dataDir, "ex1.gds")
#'
#' ## Open GDS files
#' gdsProfile <- openfn.gds(fileProfileGDS)
#'
#' ## The function returns 0L when all parameters are valid
#' RAIDS:::validateComputePCAMultiSynthetic(gdsProfile=gdsProfile,
#'     listPCA=demoPCA1KG, sampleRef=c("HG00246", "HG00325"),
#'     studyIDSyn="MyStudy", verbose=FALSE)
#'
#' ## Close GDS file (it is important to always close the GDS files)
#' closefn.gds(gdsProfile)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateComputePCAMultiSynthetic <- function(gdsProfile, listPCA, sampleRef,
                                                    studyIDSyn, verbose) {

    validateGDSClass(gds=gdsProfile, "gdsProfile")

    if(!(is.list(listPCA) &&
            "pca.unrel" %in% names(listPCA))) {
        stop("The \'listPCA\' parameter must be a list with the ",
                    "entries \'pca.unrel\'.")
    }

    if(!(is.character(sampleRef) && is.vector(sampleRef))) {
        stop("The \'sampleRef\' parameter must be a vector of ",
                    "character strings.")
    }

    ## The studyID must be a character string
    if (!(is.character(studyIDSyn) && length(studyIDSyn) == 1)) {
        stop("The \'studyIDSyn\' parameter must be a character string.")
    }

    ## Validate that verbose is a logical
    validateLogical(logical=verbose, name="verbose")

    return(0L)
}


#' @title Validate the input parameters for computeKNNRefSynthetic()
#' function
#'
#' @description The function validates the input parameters for the
#' \code{\link{computeKNNRefSynthetic}} function.
#' When a parameter is not as expected, an error message is generated.
#'
#' @param gdsProfile an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, the
#' opened Profile GDS file.
#'
#' @param listEigenvector a \code{list} with 3 entries:
#' 'sample.id', 'eigenvector.ref' and 'eigenvector'. The \code{list} represents
#' the PCA done on the 1KG reference profiles and the synthetic profiles
#' projected onto it.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations.
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier.
#' The study identifier must be present in the Profile GDS file.
#'
#' @param spRef \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset.
#'
#' @param kList  a \code{vector} of \code{integer} representing  the list of
#' values tested for the  K parameter. The K parameter represents the
#' number of neighbors used in the K-nearest neighbors analysis.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  D parameter. The D parameter represents the
#' number of dimensions used in the PCA analysis.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#' fileProfileGDS <- file.path(dataDir, "ex1.gds")
#'
#' ## Open GDS files
#' gdsProfile <- openfn.gds(fileProfileGDS)
#'
#' ## The function returns 0L when all parameters are valid
#' RAIDS:::validateComputeKNNRefSynthetic(gdsProfile=gdsProfile,
#'     listEigenvector=demoPCASyntheticProfiles,
#'     listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
#'     studyIDSyn="MyStudy", spRef=demoKnownSuperPop1KG,
#'     fieldPopInfAnc="Superpop", kList=c(10, 11, 12),
#'     pcaList=c(13, 14, 15))
#'
#' ## Close GDS file (it is important to always close the GDS files)
#' closefn.gds(gdsProfile)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateComputeKNNRefSynthetic <- function(gdsProfile, listEigenvector,
    listCatPop, studyIDSyn, spRef, fieldPopInfAnc, kList, pcaList) {

    validateGDSClass(gds=gdsProfile, "gdsProfile")

    if(!(is.list(listEigenvector) &&
            all(c("sample.id", "eigenvector.ref", "eigenvector") %in%
                names(listEigenvector)))) {
        stop("The \'listEigenvector\' parameter must be a list with 3 ",
            "entries: \'sample.id\', \'eigenvector.ref\' and \'eigenvector\'.")
    }

    if(!(is.character(listCatPop) && is.vector(listCatPop))) {
        stop("The \'listCatPop\' parameter must be a vector of ",
                "character strings.")
    }

    ## The studyIDSyn must be a character string
    if (!(is.character(studyIDSyn) && length(studyIDSyn) == 1)) {
        stop("The \'studyIDSyn\' parameter must be a character string.")
    }

    if(!(is.character(spRef) && is.vector(spRef))) {
        stop("The \'spRef\' parameter must be a vector of ",
                "character strings.")
    }

    ## The fieldPopInfAnc must be a character string
    if (!(is.character(fieldPopInfAnc) && length(fieldPopInfAnc) == 1)) {
        stop("The \'fieldPopInfAnc\' parameter must be a character string.")
    }

    ## The kList must be a vector of positive numerics
    if (!(is.numeric(kList) && is.vector(kList))) {
        stop("The \'kList\' parameter must be a vector of positive numerics",
                " representing the K-neighbors values tested.")
    }

    ## The pcaList must be a vector of positive numerics
    if (!(is.numeric(pcaList) && is.vector(pcaList))) {
        stop("The \'pcaList\' parameter must be a vector of positive ",
                "numerics representing the PCA dimensions that are tested.")
    }

    return(0L)
}


#' @title Validate the input parameters for computeKNNRefSample()
#' function
#'
#' @description The function validates the input parameters for the
#' \code{\link{computeKNNRefSample}} function.
#' When a parameter is not as expected, an error message is generated.
#'
#' @param listEigenvector a \code{list} with 3 entries:
#' 'sample.id', 'eigenvector.ref' and 'eigenvector'. The \code{list} represents
#' the PCA done on the 1KG reference profiles and the specific profile
#' projected onto it. The 'sample.id' entry must contain only one identifier
#' (one profile).
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations.
#'
#' @param spRef \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset.
#'
#' @param kList  a \code{vector} of \code{integer} representing  the list of
#' values tested for the  K parameter. The K parameter represents the
#' number of neighbors used in the K-nearest neighbors analysis.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  D parameter. The D parameter represents the
#' number of dimensions used in the PCA analysis.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' pcaSynthetic <- demoPCASyntheticProfiles
#' pcaSynthetic$sample.id <- pcaSynthetic$sample.id[1]
#'
#' ## The function returns 0L when all parameters are valid
#' RAIDS:::validateComputeKNNRefSample(listEigenvector=pcaSynthetic,
#'     listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
#'     spRef=demoKnownSuperPop1KG, fieldPopInfAnc="Superpop",
#'     kList=c(10, 11, 12), pcaList=c(13, 14, 15))
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateComputeKNNRefSample <- function(listEigenvector, listCatPop, spRef,
                                    fieldPopInfAnc, kList, pcaList) {

    if(!(is.list(listEigenvector) &&
            all(c("sample.id", "eigenvector.ref", "eigenvector") %in%
                    names(listEigenvector)))) {
        stop("The \'listEigenvector\' parameter must be a list with 3 ",
            "entries: \'sample.id\', \'eigenvector.ref\' and \'eigenvector\'.")
    }

    if(length(listEigenvector$sample.id) != 1) {
        stop("Only one profile can be present in the \'sample.id\'",
                " entry from the \'listEigenvector\' parameter.\n")
    }

    if(!(is.character(listCatPop) && is.vector(listCatPop))) {
        stop("The \'listCatPop\' parameter must be a vector of ",
                    "character strings.")
    }

    if(!(is.character(spRef) && is.vector(spRef))) {
        stop("The \'spRef\' parameter must be a vector of ",
                    "character strings.")
    }

    ## The fieldPopInfAnc must be a character string
    if (!(is.character(fieldPopInfAnc) && length(fieldPopInfAnc) == 1)) {
        stop("The \'fieldPopInfAnc\' parameter must be a character string.")
    }

    ## The kList must be a vector of positive numerics
    if (!(is.numeric(kList) && is.vector(kList))) {
        stop("The \'kList\' parameter must be a vector of positive numerics",
                " representing the K-neighbors values tested.")
    }

    ## The pcaList must be a vector of positive numerics
    if (!(is.numeric(pcaList) && is.vector(pcaList))) {
        stop("The \'pcaList\' parameter must be a vector of positive ",
                "numerics representing the PCA dimensions that are tested.")
    }

    return(0L)
}


#' @title Calculate Principal Component Analysis (PCA) on SNV genotype data set
#'
#' @description The functions calculates the principal component analysis (PCA)
#' for a list of pruned SNVs present in a Profile GDS file. The
#' \link[SNPRelate]{snpgdsPCA} function is used to do the calculation.
#'
#' @param gdsProfile an object of class \link[SNPRelate]{SNPGDSFileClass},
#' the opened Profile GDS file.
#'
#' @param refProfileIDs a \code{vector} of reference 1KG profile identifiers
#' that are present in the Profile GDS file.
#' Those profiles minus the one present in the \code{listRM} vector will be
#' used to run the PCA analysis.
#'
#' @param listRM a \code{vector} of \code{character} strings containing the
#' identifiers for the reference samples that need to be removed for the
#' PCA analysis.
#'
#' @param np a single positive \code{integer} representing the number of CPU
#' that will be used. Default: \code{1L}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missingRate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA} function
#' with "<= missingRate" only; if \code{NaN}, no missing threshold.
#' Default: \code{0.025}.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
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
#' ## Required library
#' library(SNPRelate)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' ## Profiles that should be removed from the PCA analysis
#' ## Those profiles has been used to generate the synthetic data set
#' samplesRM <- c("HG00246", "HG00325", "HG00611", "HG01173", "HG02165",
#'     "HG01112", "HG01615", "HG01968", "HG02658", "HG01850", "HG02013",
#'     "HG02465", "HG02974", "HG03814", "HG03445", "HG03689", "HG03789",
#'     "NA12751", "NA19107", "NA18548", "NA19075", "NA19475", "NA19712",
#'     "NA19731", "NA20528", "NA20908")
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#'
#' ## Compute PCA for the 1KG reference profiles excluding
#' ## the profiles used to generate the synthetic profiles
#' results <- RAIDS:::computePCARefRMMulti(gdsProfile=gdsProfile,
#'     refProfileIDs=names(demoKnownSuperPop1KG), listRM=samplesRM, np=1L,
#'     algorithm="exact", eigenCount=32L, missingRate=0.025, verbose=FALSE)
#'
#' ## The PCA on the pruned SNVs data set for selected profiles
#' head(results$pca.unrel$eigenvect)
#'
#' ## Close Profile GDS file (important)
#' closefn.gds(gdsProfile)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA
#' @encoding UTF-8
#' @keywords internal
computePCARefRMMulti <- function(gdsProfile, refProfileIDs, listRM, np=1L,
                                    algorithm="exact", eigenCount=32L,
                                    missingRate=0.025, verbose) {

    unrelatedSamples <- refProfileIDs[which(!(refProfileIDs %in% listRM))]

    listPCA <- list()

    listPCA[["pruned"]] <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))

    ## Calculate Principal Component Analysis (PCA) on SNV genotype dataset
    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsobj=gdsProfile,
                                            sample.id=unrelatedSamples,
                                            snp.id=listPCA[["pruned"]],
                                            num.thread=np,
                                            missing.rate=missingRate,
                                            algorithm=algorithm,
                                            eigen.cnt=eigenCount,
                                            verbose=verbose)

    return(listPCA)
}


#' @title Compile all the inferred ancestry results done on the
#' synthetic profiles for different D and K values in the objective of
#' selecting the optimal D and K values for a specific profile
#'
#' @description The function calculates the accuracy of the inferred ancestry
#' called done on the synthetic profiles for different D and K values. The
#' accuracy is also calculated for each super-population used to generate
#' the synthetic profiles. The known ancestry from the reference profiles
#' used to generate the synthetic profiles is required to calculate the
#' accuracy.
#'
#' @param matKNN a \code{data.frame} containing the inferred ancestry for the
#' synthetic profiles for different _K_ and _D_ values. The \code{data.frame}
#' must contained those columns: "sample.id", "D", "K" and the fourth column
#' name must correspond to the \code{predCall} argument.
#'
#' @param pedCall a \code{data.frame} containing the information about
#' the super-population information from the 1KG GDS file
#' for profiles used to generate the synthetic profiles. The \code{data.frame}
#' must contained a column named as the \code{refCall} argument.
#'
#' @param refCall a \code{character} string representing the name of the
#' column that contains the known ancestry for the reference profiles in
#' the Reference GDS file.
#'
#' @param predCall a \code{character} string representing the name of
#' the column that contains the inferred ancestry for the specified
#' profiles. The column must be present in the \code{matKNN} \code{data.frame}
#' argument.
#'
#' @param listCall a \code{vector} of \code{character} strings representing
#' the list of possible ancestry assignations.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis.
#' Default: \code{seq(3,15,1)}.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.
#' Default: \code{seq(2,15,1)}.
#'
#' @return a \code{list} containing 5 entries:
#' \describe{
#' \item{\code{dfPCA}}{ a \code{data.frame} containing statistical results
#' on all combined synthetic results done with a fixed value of \code{D} (the
#' number of dimensions). The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{median}}{ a \code{numeric} representing the median of the
#' minimum AUROC obtained (within super populations) for all combination of
#' the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{mad}}{ a \code{numeric} representing the MAD of the minimum
#' AUROC obtained (within super populations) for all combination of the fixed
#' \code{D} value and all tested \code{K} values. }
#' \item{\code{upQuartile}}{ a \code{numeric} representing the upper quartile
#' of the minimum AUROC obtained (within super populations) for all
#' combination of the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{k}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for a fixed \code{D} value. }
#' }
#' }
#' \item{\code{dfPop}}{ a \code{data.frame} containing statistical results on
#' all combined synthetic results done with different values of \code{D} (the
#' number of dimensions) and \code{K} (the number of neighbors).
#' The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{AUROC.min}}{ a \code{numeric} representing the minimum accuracy
#' obtained by grouping all the synthetic results by super-populations, for
#' the specified values of \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing the accuracy obtained
#' by grouping all the synthetic results for the specified values of \code{D}
#' and \code{K}.}
#' \item{\code{Accu.CM}}{ a \code{numeric} representing the value of accuracy
#' of the confusion matrix obtained by grouping all the synthetic results for
#' the specified values of \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{D}}{ a \code{numeric} representing the optimal \code{D} value
#' (the number of dimensions) for the specific profile.}
#' \item{\code{K}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for the specific profile.}
#' \item{\code{listD}}{ a \code{numeric} representing the optimal \code{D}
#' values (the number of dimensions) for the specific profile. More than one
#' \code{D} is possible.}
#' }
#'
#' @examples
#'
#' ## Loading demo dataset containing pedigree information for synthetic
#' ## profiles and known ancestry of the profiles used to generate the
#' ## synthetic profiles
#' data(pedSynthetic)
#'
#' ## Loading demo dataset containing the inferred ancestry results
#' ## for the synthetic data
#' data(matKNNSynthetic)
#'
#' ## Compile all the results for ancestry inference done on the
#' ## synthetic profiles for different D and K values
#' ## Select the optimal D and K values
#' results <- RAIDS:::selParaPCAUpQuartile(matKNN=matKNNSynthetic,
#'     pedCall=pedSynthetic, refCall="superPop", predCall="SuperPop",
#'     listCall=c("EAS", "EUR", "AFR", "AMR", "SAS"), kList=seq(3,15,1),
#'     pcaList=seq(2,15,1))
#' results$D
#' results$K
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom stats mad median quantile
#' @encoding UTF-8
#' @keywords internal
selParaPCAUpQuartile <- function(matKNN, pedCall, refCall,
                                    predCall, listCall, kList=seq(3,15,1),
                                    pcaList=seq(2,15,1)) {

    if (min(kList) < 3) {
        warning("A K smaller than 3 could not give robust results.\n")
    }

    tableSyn <- list()
    tableCall <- list()
    tableAUROC <- list()
    i <- 1

    ## Loop on all PCA dimension values
    for (D in pcaList) {
        matKNNCurD <- matKNN[which(matKNN$D == D), ]
        listTMP <- list()
        listTMP.AUROC <- list()
        j <- 1
        ## Loop on all k neighbor values
        for (K in kList) {
            matKNNCur <- matKNNCurD[which(matKNNCurD$K == K), ]
            ## Calculate accuracy for fixed D and K values
            res <- computeSyntheticConfMat(matKNN=matKNNCur,
                    matKNNAncestryColumn=predCall, pedCall=pedCall,
                    pedCallAncestryColumn=refCall, listCall=listCall)
            resROC <- computeSyntheticROC(matKNN=matKNNCur,
                    matKNNAncestryColumn=predCall, pedCall=pedCall,
                    pedCallAncestryColumn=refCall, listCall=listCall)

            df <- data.frame(D=D, K=K, AUROC.min=min(resROC$matAUROC.Call$AUC),
                                AUROC=resROC$matAUROC.All$ROC.AUC,
                                Accu.CM=res$matAccuracy$Accu.CM)

            listTMP[[j]] <- df
            listTMP.AUROC[[j]] <- resROC$matAUROC.Call
            j <- j + 1
        }
        df <- do.call(rbind, listTMP)

        tableCall[[i]] <- df
        tableAUROC[[i]] <- do.call(rbind, listTMP.AUROC)
        maxAUROC <- max(df[df$K %in% kList, "AUROC.min"])
        kMax <- df[df$K %in% kList & abs(df$AUROC.min-maxAUROC) < 1e-3, "K"]
        kV <- kMax[(length(kMax) + length(kMax)%%2)/2]
        dfPCA <- data.frame(D=D,
            median=median(df[df$K %in% kList, "AUROC.min"]),
            mad=mad(df[df$K %in% kList, "AUROC.min"]),
            upQuartile=quantile(df[df$K %in% kList, "AUROC.min"], 0.75), K=kV)
        tableSyn[[i]] <- dfPCA
        i <- i + 1
    }

    dfPCA <- do.call(rbind, tableSyn)
    dfCall <- do.call(rbind, tableCall)
    dfAUROC <- do.call(rbind, tableAUROC)
    selD <- dfPCA$D[which.max(dfPCA$upQuartile)]
    selK <- dfPCA$K[which.max(dfPCA$upQuartile)]
    tmp <- max(dfPCA$upQuartile)
    listD <- dfPCA$D[which(abs(dfPCA$upQuartile - tmp) < 1e-3)]

    res <- list(dfPCA=dfPCA, dfPop=dfCall, dfAUROC=dfAUROC,
                D=selD, K=selK, listD=listD)
    return(res)
}

#' @title Run most steps leading to the ancestry inference call on a
#' specific profile (RNA or DNA)
#'
#' @description This function runs most steps leading to the ancestry inference
#' call on a specific profile. First, the function creates the Profile GDS file
#' for the specific profile using the information from a RDS Sample
#' description file and the Population reference GDS file.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Population Reference GDS file.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Population Reference SNV Annotation GDS file.
#' This parameter is RNA specific.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param currentProfile a \code{character} string representing the profile
#' identifier.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Profile files will be created.
#' Default: \code{NULL}.
#'
#' @param pathOut a \code{character} string representing the path to
#' the directory where the output files are created.
#'
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section.
#'
#' @param syntheticRefDF a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The \code{data.frame} must have those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @param studyDFSyn a \code{data.frame} containing the information about the
#' synthetic data to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param listProfileRef a \code{vector} of \code{character} string
#' representing the
#' identifiers of the selected 1KG profiles that will be used as reference to
#' generate the synthetic profiles.
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "DNA" and "RNA". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"DNA"}.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param blockTypeID a \code{character} string corresponding to the block
#' type used to extract the block identifiers. The block type must be
#' present in the GDS Reference Annotation file.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{FALSE}.
#'
#' @return The integer \code{0L} when successful. See details section for
#' more information about the generated output files.
#'
#' @details
#'
#' The runWrapperAncestry() function generates 3 types of files
#' in the \code{pathOut} directory:
#' \describe{
#' \item{Ancestry Inference}{ The ancestry inference CSV file
#' (".Ancestry.csv" file)}
#' \item{Inference Informaton}{ The inference information RDS file
#' (".infoCall.rds" file)}
#' \item{Synthetic Information}{ The parameter information RDS files
#' from the synthetic inference ("KNN.synt.*.rds" files in a sub-directory)}
#' }
#'
#' In addition, a sub-directory (named using the profile ID) is
#' also created.
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
#' ## Required library for GDS
#' library(SNPRelate)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' #################################################################
#' ## Load the information about the profile
#' #################################################################
#' data(demoPedigreeEx1)
#' head(demoPedigreeEx1)
#'
#' #################################################################
#' ## The 1KG GDS file and the 1KG SNV Annotation GDS file
#' ## need to be located in the same directory
#' ## Note that the 1KG GDS file used for this example is a
#' ## simplified version and CANNOT be used for any real analysis
#' #################################################################
#' path1KG <- file.path(dataDir, "tests")
#'
#' fileReferenceGDS  <- file.path(path1KG, "ex1_good_small_1KG.gds")
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#'
#' #################################################################
#' ## The Sample SNP pileup files (one per sample) need
#' ## to be located in the same directory.
#' #################################################################
#' pathGeno <- file.path(dataDir, "example", "snpPileup")
#'
#' #################################################################
#' ## The path where the Profile GDS Files (one per sample)
#' ## will be created need to be specified.
#' #################################################################
#'
#' pathProfileGDS <- file.path(tempdir(), "outTest.tmp")
#'
#' pathOut <- file.path(tempdir(), "resTest.out")
#'
#' #################################################################
#' ## A data frame containing general information about the study
#' ## is also required. The data frame must have
#' ## those 3 columns: "studyID", "study.desc", "study.platform"
#' #################################################################
#' studyDF <- data.frame(study.id="MYDATA",
#'                         study.desc="Description",
#'                         study.platform="PLATFORM",
#'                         stringsAsFactors=FALSE)
#'
#' ####################################################################
#' ## Fix seed to ensure reproducible results
#' ####################################################################
#' set.seed(3043)
#'
#' gdsReference <- snpgdsOpen(fileReferenceGDS)
#' dataRef <- select1KGPop(gdsReference, nbProfiles=2L)
#' closefn.gds(gdsReference)
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#'
#'     studyDFSyn <- data.frame(study.id=paste0(studyDF$study.id, ".Synthetic"),
#'         study.desc=paste0(studyDF$study.id, " synthetic data"),
#'         study.platform=studyDF$study.platform, stringsAsFactors=FALSE)
#'
#'     listProfileRef <- dataRef$sample.id
#'     profileFile <- file.path(pathProfileGDS, "ex1.gds")
#'
#'     \dontrun{
#'
#'         dir.create(pathProfileGDS)
#'         dir.create(pathOut)
#'         file.copy(file.path(dataDir, "tests", "ex1_demo.gds"), profileFile)
#'
#'         gdsReference <- snpgdsOpen(fileReferenceGDS)
#'         gdsRefAnnot <- openfn.gds(fileAnnotGDS)
#'
#'         RAIDS:::runProfileAncestry(gdsReference=gdsReference,
#'              gdsRefAnnot=gdsRefAnnot,
#'              studyDF=studyDF, currentProfile=ped[1,"Name.ID"],
#'              pathProfileGDS=pathProfileGDS,
#'              pathOut=pathOut,
#'              chrInfo=chrInfo,
#'              syntheticRefDF=dataRef,
#'              studyDFSyn=studyDFSyn,
#'              listProfileRef=listProfileRef,
#'              studyType="DNA")
#'
#'         closefn.gds(gdsReference)
#'         closefn.gds(gdsRefAnnot)
#'
#'         unlink(pathProfileGDS, recursive=TRUE, force=TRUE)
#'         unlink(pathOut, recursive=TRUE, force=TRUE)
#'
#'     }
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @keywords internal
runProfileAncestry <- function(gdsReference, gdsRefAnnot, studyDF,
    currentProfile, pathProfileGDS, pathOut, chrInfo, syntheticRefDF,
    studyDFSyn, listProfileRef, studyType=c("DNA", "RNA"),
    np=1L, blockTypeID=NULL, verbose=FALSE) {

    studyType <- arg_match(studyType)

    pruningSample(gdsReference=gdsReference, currentProfile=currentProfile,
            studyID=studyDF$study.id, pathProfileGDS=pathProfileGDS, np=np)

    fileGDSProfile <- file.path(pathProfileGDS,
                            paste0(currentProfile, ".gds"))
    add1KG2SampleGDS(gdsReference=gdsReference, fileProfileGDS=fileGDSProfile,
                    currentProfile=currentProfile, studyID=studyDF$study.id)

    addStudy1Kg(gdsReference, fileGDSProfile)

    gdsProfile <- openfn.gds(fileGDSProfile, readonly=FALSE)

    estimateAllelicFraction(gdsReference=gdsReference, gdsProfile=gdsProfile,
                currentProfile=currentProfile, studyID=studyDF$study.id,
                chrInfo=chrInfo, studyType=studyType, gdsRefAnnot=gdsRefAnnot,
                blockID=blockTypeID, verbose=verbose)
    closefn.gds(gdsProfile)

    ## Add information related to the synthetic profiles in Profile GDS file
    prepSynthetic(fileProfileGDS=fileGDSProfile,
            listSampleRef=listProfileRef,  profileID=currentProfile,
            studyDF=studyDFSyn, prefix="1", verbose=verbose)

    resG <- syntheticGeno(gdsReference=gdsReference, gdsRefAnnot=gdsRefAnnot,
                    fileProfileGDS=fileGDSProfile, profileID=currentProfile,
                    listSampleRef=listProfileRef, prefix="1")

    if(! file.exists(pathOut)) {
        dir.create(pathOut)
    }
    spRef <- getRef1KGPop(gdsReference, "superPop")
    sampleRM <- splitSelectByPop(syntheticRefDF)

    pathOutProfile <- file.path(pathOut, currentProfile)
    if(! file.exists(pathOutProfile)) {
        dir.create(pathOutProfile)
    }

    ## Open the Profile GDS file
    gdsProfile <- snpgdsOpen(fileGDSProfile)


    ## This variable will contain the results from the PCA analyses
    ## For each row of the sampleRM matrix
    apply(t(t(seq_len(nrow(sampleRM)))), 1, FUN=function(x, sampleRM,
                                        gdsProfile, studyDFSyn, spRef,
                                        pathOutProfile, currentProfile) {
            synthKNN <- computePoolSyntheticAncestryGr(gdsProfile=gdsProfile,
                            sampleRM=sampleRM[x,],
                            studyIDSyn=studyDFSyn$study.id,
                            np=np, spRef=spRef, eigenCount=15L,
                            verbose=verbose)

            ## Results are saved
            saveRDS(synthKNN$matKNN, file.path(pathOutProfile,
                        paste0("KNN.synt.", currentProfile, ".", x, ".rds")))
            return(NULL)
        }, sampleRM=sampleRM, gdsProfile=gdsProfile,
        studyDFSyn=studyDFSyn, spRef=spRef, pathOutProfile=pathOutProfile,
        currentProfile=currentProfile)


    ## Directory where the KNN results have been saved
    pathKNN <- file.path(pathOut, currentProfile)
    listFilesName <- dir(file.path(pathKNN), ".rds")
    ## List of the KNN result files from PCA on synthetic data
    listFiles <- file.path(file.path(pathKNN) , listFilesName)

    resCall <- computeAncestryFromSyntheticFile(gdsReference=gdsReference,
                        gdsProfile=gdsProfile, listFiles=listFiles,
                        currentProfile=currentProfile, spRef=spRef,
                        studyIDSyn=studyDFSyn$study.id, np=np)

    saveRDS(resCall, file.path(pathOut,
                            paste0(currentProfile, ".infoCall", ".rds")))

    write.csv(x=resCall$Ancestry, file=file.path(pathOut,
        paste0(currentProfile, ".Ancestry",".csv")), quote=FALSE,
        row.names=FALSE)

    ## Close Profile GDS file (important)
    closefn.gds(gdsProfile)

    return(0L)
}

#' @title Run most steps leading to the ancestry inference call on a
#' specific profile (LD or geneAware)
#'
#' @description This function runs most steps leading to the ancestry inference
#' call on a specific profile. First, the function creates the Profile GDS file
#' for the specific profile using the information from a RDS Sample
#' description file and the Population reference GDS file.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Population Reference GDS file.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Population Reference SNV Annotation GDS file.
#' This parameter is RNA specific.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param currentProfile a \code{character} string representing the profile
#' identifier.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Profile files will be created.
#' Default: \code{NULL}.
#'
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section.
#'
#' @param syntheticRefDF a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The \code{data.frame} must have those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @param studyDFSyn a \code{data.frame} containing the information about the
#' synthetic data to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param listProfileRef a \code{vector} of \code{character} string
#' representing the
#' identifiers of the selected 1KG profiles that will be used as reference to
#' generate the synthetic profiles.
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "LD" and "GeneAware". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"LD"}.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param blockTypeID a \code{character} string corresponding to the block
#' type used to extract the block identifiers. The block type must be
#' present in the GDS Reference Annotation file.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{FALSE}.
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{pcaSample}}{ a \code{list} containing the information related
#' to the eigenvectors. The \code{list} contains those 3 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{eigenvector.ref}}{ a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{\code{eigenvector}}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current profile projected on the PCA from the
#' reference profiles.}
#' }
#' }
#' \item{\code{paraSample}}{ a \code{list} containing the results with
#' different \code{D} and \code{K} values that lead to optimal parameter
#' selection. The \code{list} contains those entries:
#' \describe{
#' \item{\code{dfPCA}}{ a \code{data.frame} containing statistical results
#' on all combined synthetic results done with a fixed value of \code{D} (the
#' number of dimensions). The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{median}}{ a \code{numeric} representing the median of the
#' minimum AUROC obtained (within super populations) for all combination of
#' the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{mad}}{ a \code{numeric} representing the MAD of the minimum
#' AUROC obtained (within super populations) for all combination of the fixed
#' \code{D} value and all tested \code{K} values. }
#' \item{\code{upQuartile}}{ a \code{numeric} representing the upper quartile
#' of the minimum AUROC obtained (within super populations) for all
#' combination of the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{k}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for a fixed \code{D} value. }
#' }
#' }
#' \item{\code{dfPop}}{ a \code{data.frame} containing statistical results on
#' all combined synthetic results done with different values of \code{D} (the
#' number of dimensions) and \code{K} (the number of neighbors).
#' The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{AUROC.min}}{ a \code{numeric} representing the minimum accuracy
#' obtained by grouping all the synthetic results by super-populations, for
#' the specified values of \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing the accuracy obtained
#' by grouping all the synthetic results for the specified values of \code{D}
#' and \code{K}.}
#' \item{\code{Accu.CM}}{ a \code{numeric} representing the value of accuracy
#' of the confusion matrix obtained by grouping all the synthetic results for
#' the specified values of \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{dfAUROC}}{ a \code{data.frame} the summary of the results by
#' super-population. The \code{data.frame} contains
#' those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{Call}}{ a \code{character} string representing the
#' super-population.}
#' \item{\code{L}}{ a \code{numeric} representing the lower value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing  the AUROC obtained for the
#' fixed values of super-population, \code{D} and \code{K}.}
#' \item{\code{H}}{ a \code{numeric} representing the higher value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{D}}{ a \code{numeric} representing the optimal \code{D} value
#' (the number of dimensions) for the specific profile.}
#' \item{\code{K}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for the specific profile.}
#' \item{\code{listD}}{ a \code{numeric} representing the optimal \code{D}
#' values (the number of dimensions) for the specific profile. More than one
#' \code{D} is possible.}
#' }
#' }
#' \item{\code{KNNSample}}{  a \code{data.frame} containing the inferred ancestry
#' for different values of \code{K} and \code{D}. The \code{data.frame}
#' contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' }
#' }
#' \item{\code{KNNSynthetic}}{  a \code{data.frame} containing the inferred ancestry
#' for each synthetic data for different values of \code{K} and \code{D}.
#' The \code{data.frame}
#' contains those columns: "sample.id", "D", "K", "infer.superPop", "ref.superPop"
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current synthetic data.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{infer.superPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' \item{\code{ref.superPop}}{ a \code{character} string representing the known
#' ancestry from the reference}
#' }
#' }
#' \item{\code{Ancestry}}{ a \code{data.frame} containing the inferred
#' ancestry for the current profile. The \code{data.frame} contains those
#' columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry.}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry.}
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry.}
#' }
#' }
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
#' ## Required library for GDS
#' library(SNPRelate)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' #################################################################
#' ## Load the information about the profile
#' #################################################################
#' data(demoPedigreeEx1)
#' head(demoPedigreeEx1)
#'
#' #################################################################
#' ## The 1KG GDS file and the 1KG SNV Annotation GDS file
#' ## need to be located in the same directory
#' ## Note that the 1KG GDS file used for this example is a
#' ## simplified version and CANNOT be used for any real analysis
#' #################################################################
#' path1KG <- file.path(dataDir, "tests")
#'
#' fileReferenceGDS  <- file.path(path1KG, "ex1_good_small_1KG.gds")
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#'
#' #################################################################
#' ## The Sample SNP pileup files (one per sample) need
#' ## to be located in the same directory.
#' #################################################################
#' pathGeno <- file.path(dataDir, "example", "snpPileup")
#'
#' #################################################################
#' ## The path where the Profile GDS Files (one per sample)
#' ## will be created need to be specified.
#' #################################################################
#'
#' pathProfileGDS <- file.path(tempdir(), "outTest.tmp")
#'
#'
#' #################################################################
#' ## A data frame containing general information about the study
#' ## is also required. The data frame must have
#' ## those 3 columns: "studyID", "study.desc", "study.platform"
#' #################################################################
#' studyDF <- data.frame(study.id="MYDATA",
#'                         study.desc="Description",
#'                         study.platform="PLATFORM",
#'                         stringsAsFactors=FALSE)
#'
#' ####################################################################
#' ## Fix seed to ensure reproducible results
#' ####################################################################
#' set.seed(3043)
#'
#' dataRef <- select1KGPopForSynthetic(fileReferenceGDS, nbProfiles=2L)
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#'
#'     studyDFSyn <- data.frame(study.id=paste0(studyDF$study.id, ".Synthetic"),
#'         study.desc=paste0(studyDF$study.id, " synthetic data"),
#'         study.platform=studyDF$study.platform, stringsAsFactors=FALSE)
#'
#'     listProfileRef <- dataRef$sample.id
#'     profileFile <- file.path(pathProfileGDS, "ex1.gds")
#'
#'     \dontrun{
#'
#'         dir.create(pathProfileGDS)
#'         file.copy(file.path(dataDir, "tests", "ex1_demo.gds"), profileFile)
#'
#'         gdsReference <- snpgdsOpen(fileReferenceGDS)
#'         gdsRefAnnot <- openfn.gds(fileAnnotGDS)
#'
#'         res <- RAIDS:::profileAncestry(gdsReference=gdsReference,
#'              gdsRefAnnot=gdsRefAnnot,
#'              studyDF=studyDF, currentProfile=demoPedigreeEx1[1,"Name.ID"],
#'              pathProfileGDS=pathProfileGDS,
#'              chrInfo=chrInfo,
#'              syntheticRefDF=dataRef,
#'              studyDFSyn=studyDFSyn,
#'              listProfileRef=listProfileRef,
#'              studyType="LD")
#'
#'         closefn.gds(gdsReference)
#'         closefn.gds(gdsRefAnnot)
#'
#'
#'     }
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @keywords internal
profileAncestry <- function(gdsReference, gdsRefAnnot, studyDF,
            currentProfile, pathProfileGDS, chrInfo, syntheticRefDF,
            studyDFSyn, listProfileRef, studyType=c("LD", "GeneAware"),
            np=1L, blockTypeID=NULL, verbose=FALSE) {
    # This part can be share with runProfileAncestry
    studyType <- arg_match(studyType)
    if(verbose){
        message("pruningSample start ", Sys.time())
    }
    pruningSample(gdsReference=gdsReference, currentProfile=currentProfile,
        studyID=studyDF$study.id, pathProfileGDS=pathProfileGDS, np=np)
    if(verbose){
        message("pruningSample end ", Sys.time())
    }
    fileGDSProfile <- file.path(pathProfileGDS,
                                    paste0(currentProfile, ".gds"))

    add1KG2SampleGDS(gdsReference=gdsReference, fileProfileGDS=fileGDSProfile,
            currentProfile=currentProfile, studyID=studyDF$study.id)
    if(verbose){
        message("add1KG start ", Sys.time())
    }
    addStudy1Kg(gdsReference, fileGDSProfile)

    gdsProfile <- openfn.gds(fileGDSProfile, readonly=FALSE)
    # Change for the old studyType
    studyTypeLeg <- ifelse(studyType=="LD", "DNA", "RNA")
    estimateAllelicFraction(gdsReference=gdsReference, gdsProfile=gdsProfile,
            currentProfile=currentProfile, studyID=studyDF$study.id,
            chrInfo=chrInfo, studyType=studyTypeLeg, gdsRefAnnot=gdsRefAnnot,
            blockID=blockTypeID, verbose=verbose)
    closefn.gds(gdsProfile)

    ## Add information related to the synthetic profiles in Profile GDS file
    prepSynthetic(fileProfileGDS=fileGDSProfile,
        listSampleRef=listProfileRef,  profileID=currentProfile,
        studyDF=studyDFSyn, prefix="1", verbose=verbose)
    if(verbose){
        message("syntheticGeno start ", Sys.time())
    }
    resG <- syntheticGeno(gdsReference=gdsReference, gdsRefAnnot=gdsRefAnnot,
                fileProfileGDS=fileGDSProfile, profileID=currentProfile,
                listSampleRef=listProfileRef, prefix="1")

    # if(! file.exists(pathOut)) {
    #     dir.create(pathOut)
    # }
    #
    spRef <- getRef1KGPop(gdsReference, "superPop")
    sampleRM <- splitSelectByPop(syntheticRefDF)
    #
    # pathOutProfile <- file.path(pathOut, currentProfile)
    # if(! file.exists(pathOutProfile)) {
    #     dir.create(pathOutProfile)
    # }

    ##### End share with runProfileAncestry

    ## Open the Profile GDS file
    gdsProfile <- snpgdsOpen(fileGDSProfile)

    if(verbose){
        message("SyntheticAncestry start ", Sys.time())
    }
    ## This variable will contain the results from the PCA analyses
    ## For each row of the sampleRM matrix
    resSyn <- lapply(seq_len(nrow(sampleRM)), FUN=function(x, sampleRM,
                        gdsProfile, studyDFSyn, spRef, currentProfile) {
        synthKNN <- computePoolSyntheticAncestryGr(gdsProfile=gdsProfile,
                            sampleRM=sampleRM[x,],
                            studyIDSyn=studyDFSyn$study.id,
                            np=np, spRef=spRef, eigenCount=15L,
                            verbose=verbose)
        ## Results are saved
        # saveRDS(synthKNN$matKNN, file.path(pathOutProfile,
        #    paste0("KNN.synt.", currentProfile, ".", x, ".rds")))
        return(synthKNN$matKNN)
    }, sampleRM=sampleRM, gdsProfile=gdsProfile, studyDFSyn=studyDFSyn,
            spRef=spRef, currentProfile=currentProfile)

    if(verbose){
        message("SyntheticAncestry end ", Sys.time())
    }
    resSyn <- do.call(rbind, resSyn)
    ## Extract the super-population information from the 1KG GDS file
    ## for profiles associated to the synthetic study
    pedSyn <- prepPedSynthetic1KG(gdsReference=gdsReference,
        gdsSample=gdsProfile, studyID=studyDFSyn$study.id, popName="superPop")


    # idCur <- matrix(unlist(strsplit(resSyn$sample.id, "\\.")), nr=4)
    # resKNN$superPopObs <- df[idCur[3,], 3]
    # ## Directory where the KNN results have been saved
    # pathKNN <- file.path(pathOut, currentProfile)
    # listFilesName <- dir(file.path(pathKNN), ".rds")
    # ## List of the KNN result files from PCA on synthetic data
    # listFiles <- file.path(file.path(pathKNN) , listFilesName)

    resCall <- computeAncestryFromSynthetic(gdsReference=gdsReference,
                    gdsProfile=gdsProfile, syntheticKNN=resSyn,
                    pedSyn=pedSyn, currentProfile=currentProfile, spRef=spRef,
                    studyIDSyn=studyDFSyn$study.id, np=np)
    if(verbose){
        message("Ancestry end ", Sys.time())
    }
    # saveRDS(resCall, file.path(pathOut,
    #                            paste0(currentProfile, ".infoCall", ".rds")))
    #
    # write.csv(x=resCall$Ancestry, file=file.path(pathOut,
    #       paste0(currentProfile, ".Ancestry",".csv")), quote=FALSE,
    #           row.names=FALSE)

    ## Close Profile GDS file (important)
    closefn.gds(gdsProfile)
    resSyn[[paste0("ref.superPop")]] <- pedSyn[resSyn$sample.id, "superPop"]

    colnames(resSyn) <- c("sample.id", "D", "K", "infer.superPop",
                            "ref.superPop")

    res <- list(pcaSample=resCall$pcaSample, # PCA of the profile + 1KG
            paraSample=resCall$paraSample, # Result of the parameter selection
            KNNSample=resCall$KNNSample$matKNN, # KNN for the profile
            KNNSynthetic=resSyn, # KNN results for synthetic data
            Ancestry=resCall$Ancestry) # the ancestry call fo the profile

    return(res)
}

#' @title Select the optimal K and D parameters using the synthetic data and
#' infer the ancestry of a specific profile
#'
#' @description The function select the optimal K and D parameters for a
#' specific profile. The results on the synthetic data are used for the
#' parameter selection. Once the optimal parameters are selected, the
#' ancestry is inferred for the specific profile.
#'
#' @param gdsReference an object of class \link[gdsfmt]{gds.class} (a GDS
#' file), the opened 1KG GDS file.
#'
#' @param gdsProfile an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Profile GDS file.
#'
#' @param syntheticKNN a \code{vector} of \code{character} strings representing
#' the name of files that contain the results of ancestry inference done on
#' the synthetic profiles for multiple values of _D_ and _K_. The files must
#' exist.
#'
#' @param pedSyn a \code{data.frame} containing the columns extracted from the
#' GDS Sample 'study.annot' node with a extra column named as the 'popName'
#' parameter that has been extracted from the 1KG GDS 'sample.annot' node.
#'
#' @param currentProfile a \code{character} string representing the profile
#' identifier of the current profile on which ancestry will be inferred.
#'
#' @param spRef a \code{vector} of \code{character} strings representing the
#' known super population ancestry for the 1KG profiles. The 1KG profile
#' identifiers are used as names for the \code{vector}.
#'
#' @param studyIDSyn a \code{character} string corresponding to the study
#' identifier. The study identifier must be present in the GDS Sample file.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @param listCatPop a \code{vector} of \code{character} string
#' representing the list of possible ancestry assignations. Default:
#' \code{("EAS", "EUR", "AFR", "AMR", "SAS")}.
#'
#' @param fieldPopIn1KG a \code{character} string representing the name of the
#' column that contains the known ancestry for the reference profiles in
#' the Reference GDS file.
#'
#' @param fieldPopInfAnc a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' profiles. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.  If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#' Default: \code{seq(2,15,1)}.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missingRate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA}
#' with "<= missingRate" only; if \code{NaN}, no missing threshold.
#' Default: \code{NaN}.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{FALSE}.
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{pcaSample}}{ a \code{list} containing the information related
#' to the eigenvectors. The \code{list} contains those 3 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{eigenvector.ref}}{ a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{\code{eigenvector}}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current profile projected on the PCA from the
#' reference profiles.}
#' }
#' }
#' \item{\code{paraSample}}{ a \code{list} containing the results with
#' different \code{D} and \code{K} values that lead to optimal parameter
#' selection. The \code{list} contains those entries:
#' \describe{
#' \item{\code{dfPCA}}{ a \code{data.frame} containing statistical results
#' on all combined synthetic results done with a fixed value of \code{D} (the
#' number of dimensions). The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{median}}{ a \code{numeric} representing the median of the
#' minimum AUROC obtained (within super populations) for all combination of
#' the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{mad}}{ a \code{numeric} representing the MAD of the minimum
#' AUROC obtained (within super populations) for all combination of the fixed
#' \code{D} value and all tested \code{K} values. }
#' \item{\code{upQuartile}}{ a \code{numeric} representing the upper quartile
#' of the minimum AUROC obtained (within super populations) for all
#' combination of the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{k}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for a fixed \code{D} value. }
#' }
#' }
#' \item{\code{dfPop}}{ a \code{data.frame} containing statistical results on
#' all combined synthetic results done with different values of \code{D} (the
#' number of dimensions) and \code{K} (the number of neighbors).
#' The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{AUROC.min}}{ a \code{numeric} representing the minimum accuracy
#' obtained by grouping all the synthetic results by super-populations, for
#' the specified values of \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing the accuracy obtained
#' by grouping all the synthetic results for the specified values of \code{D}
#' and \code{K}.}
#' \item{\code{Accu.CM}}{ a \code{numeric} representing the value of accuracy
#' of the confusion matrix obtained by grouping all the synthetic results for
#' the specified values of \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{dfAUROC}}{ a \code{data.frame} the summary of the results by
#' super-population. The \code{data.frame} contains
#' those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{Call}}{ a \code{character} string representing the
#' super-population.}
#' \item{\code{L}}{ a \code{numeric} representing the lower value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing  the AUROC obtained for the
#' fixed values of super-population, \code{D} and \code{K}.}
#' \item{\code{H}}{ a \code{numeric} representing the higher value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{D}}{ a \code{numeric} representing the optimal \code{D} value
#' (the number of dimensions) for the specific profile.}
#' \item{\code{K}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for the specific profile.}
#' \item{\code{listD}}{ a \code{numeric} representing the optimal \code{D}
#' values (the number of dimensions) for the specific profile. More than one
#' \code{D} is possible.}
#' }
#' }
#' \item{\code{KNNSample}}{ a \code{list} containing the inferred ancestry
#' using different \code{D} and \code{K} values. The \code{list} contains
#' those entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{matKNN}}{ a \code{data.frame} containing the inferred ancestry
#' for different values of \code{K} and \code{D}. The \code{data.frame}
#' contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' }
#' }
#' }
#' }
#' \item{\code{Ancestry}}{ a \code{data.frame} containing the inferred
#' ancestry for the current profile. The \code{data.frame} contains those
#' columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry.}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry.}
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry.}
#' }
#' }
#' }
#'
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
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## The Reference GDS file
#' path1KG <- system.file("extdata/tests", package="RAIDS")
#'
#' ## Open the Reference GDS file
#' gdsRef <- snpgdsOpen(file.path(path1KG, "ex1_good_small_1KG.gds"))
#'
#' ## Path to the demo synthetic results files
#' ## List of the KNN result files from PCA run on synthetic data
#' dataDirRes <- system.file("extdata/demoAncestryCall/ex1", package="RAIDS")
#' listFilesName <- dir(file.path(dataDirRes), ".rds")
#' listFiles <- file.path(file.path(dataDirRes) , listFilesName)
#' syntheticKNN <- lapply(listFiles, FUN=function(x){return(readRDS(x))})
#' syntheticKNN <- do.call(rbind, syntheticKNN)
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#' \dontrun{
#'     pedSyn <- RAIDS:::prepPedSynthetic1KG(gdsReference=gdsRef,
#'               gdsSample=gdsProfile, studyID=studyID, popName="superPop")
#'
#'     ## Run the ancestry inference on one profile called 'ex1'
#'     ## The values of K and D used for the inference are selected using the
#'     ## synthetic results listFiles=listFiles,
#'     resCall <- RAIDS:::computeAncestryFromSynthetic(gdsReference=gdsRef,
#'                                 gdsProfile=gdsProfile,
#'                                 syntheticKNN = syntheticKNN,
#'                                 pedSyn = pedSyn,
#'                                 currentProfile=c("ex1"),
#'                                 spRef=demoKnownSuperPop1KG,
#'                                 studyIDSyn=studyID, np=1L)
#'
#'     ## The ancestry called with the optimal D and K values
#'     resCall$Ancestry
#' }
#' ## Close the GDS files (important)
#' closefn.gds(gdsProfile)
#' closefn.gds(gdsRef)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @export
computeAncestryFromSynthetic <- function(gdsReference, gdsProfile,
                                            syntheticKNN,
                                            pedSyn,
                                            currentProfile,
                                            spRef,
                                            studyIDSyn,
                                            np=1L,
                                            listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                                            fieldPopIn1KG="superPop",
                                            fieldPopInfAnc="SuperPop",
                                            kList=seq(2, 15, 1),
                                            pcaList=seq(2, 15, 1),
                                            algorithm=c("exact", "randomized"),
                                            eigenCount=32L,
                                            missingRate=NaN, verbose=FALSE) {

    if(is.null(pcaList)) {
        pcaList <- seq(2, 15, 1)
    }

    if (is.null(kList)) {
        kList <- seq(2, 15, 1)
    }

    ## Validate input parameters
    # validateComputeAncestryFromSynthetic(gdsReference=gdsReference,
    #                                          gdsProfile=gdsProfile, syntheticKNN=syntheticKNN,
    #                                          pedSyn=pedSyn,
    #                                          currentProfile=currentProfile, spRef=spRef, studyIDSyn=studyIDSyn,
    #                                          np=np, listCatPop=listCatPop, fieldPopIn1KG=fieldPopIn1KG,
    #                                          fieldPopInfAnc=fieldPopInfAnc, kList=kList, pcaList=pcaList,
    #                                          algorithm=algorithm, eigenCount=eigenCount, missingRate=missingRate,
    #                                          verbose=verbose)

    ## Matches a character method against a table of candidate values
    algorithm <- arg_match(algorithm)



    ## Compile all the inferred ancestry results for different values of
    ## D and K to select the optimal parameters
    listParaSample <- selParaPCAUpQuartile(matKNN=syntheticKNN,
                                           pedCall=pedSyn, refCall=fieldPopIn1KG, predCall=fieldPopInfAnc,
                                           listCall=listCatPop)

    ## Project profile on the PCA created with the reference profiles
    listPCAProfile <- computePCARefSample(gdsProfile=gdsProfile,
                                          currentProfile=currentProfile, studyIDRef="Ref.1KG", np=np,
                                          algorithm=algorithm, eigenCount=eigenCount, missingRate=missingRate,
                                          verbose=verbose)

    ## Run a k-nearest neighbors analysis on one specific profile
    listKNNSample <- computeKNNRefSample(listEigenvector=listPCAProfile,
                                         listCatPop=listCatPop, spRef=spRef, fieldPopInfAnc=fieldPopInfAnc,
                                         kList=kList, pcaList=pcaList)

    ## The ancestry call for the current profile
    resCall <- listKNNSample$matKNN[
        which(listKNNSample$matKNN$D == listParaSample$D &
                  listKNNSample$matKNN$K == listParaSample$K ),]
    colnames(listParaSample$dfAUROC) <- c("D", "K", "Call", "L", "AUROC", "H")
    res <- list(pcaSample=listPCAProfile, # PCA of the profile + 1KG
                paraSample=listParaSample, # Result of the parameter selection
                KNNSample=listKNNSample, # KNN for the profile
                Ancestry=resCall) # the ancestry call fo the profile

    return(res)
}


#' @title Run most steps leading to the ancestry inference call
#' on a specific profile (RNA or DNA)
#'
#' @description This function runs most steps leading to the ancestry inference
#' call on a specific profile. First, the function creates the Profile GDS file
#' for the specific profile using the information from a RDS Sample
#' description file and the Population reference GDS file.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Profile files will be created.
#' Default: \code{NULL}.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the VCF output of SNP-pileup for each sample. The
#' SNP-pileup files must be compressed (gz files) and have the name identifiers
#' of the samples. A sample with "Name.ID" identifier would have an
#' associated file called
#' if genoSource is "VCF", then "Name.ID.vcf.gz",
#' if genoSource is "generic", then "Name.ID.generic.txt.gz"
#' if genoSource is "snp-pileup", then "Name.ID.txt.gz".
#'
#' @param pathOut a \code{character} string representing the path to
#' the directory where the output files are created.
#'
#' @param fileReferenceGDS  a \code{character} string representing the file
#' name of the Reference GDS file. The file must exist.
#'
#' @param fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Reference GDS Annotation file. The file must exist.
#'
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section.
#'
#' @param syntheticRefDF a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The \code{data.frame} must have those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "DNA" and "RNA". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"DNA"}.
#'
#' @param genoSource a \code{character} string with two possible values:
#' 'snp-pileup', 'generic' or 'VCF'. It specifies if the genotype files
#' are generated by snp-pileup (Facets) or are a generic format CSV file
#' with at least those columns:
#' 'Chromosome', 'Position', 'Ref', 'Alt', 'Count', 'File1R' and 'File1A'.
#' The 'Count' is the depth at the specified position;
#' 'FileR' is the depth of the reference allele and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally the file can be a VCF file with at least those genotype
#' fields: GT, AD, DP.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param blockTypeID a \code{character} string corresponding to the block
#' type used to extract the block identifiers. The block type must be
#' present in the GDS Reference Annotation file.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{FALSE}.
#'
#' @return The integer \code{0L} when successful. See details section for
#' more information about the generated output files.
#'
#' @details
#'
#' The runWrapperAncestry() function generates 3 types of files
#' in the \code{pathOut} directory.
#' \describe{
#' \item{Ancestry Inference}{ The ancestry inference CSV file
#' (".Ancestry.csv" file)}
#' \item{Inference Informaton}{ The inference information RDS file
#' (".infoCall.rds" file)}
#' \item{Synthetic Information}{ The parameter information RDS files
#' from the synthetic inference ("KNN.synt.*.rds" files in a sub-directory)}
#' }
#'
#' In addition, a sub-directory (named using the profile ID) is
#' also created.
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
#' ## Required library for GDS
#' library(SNPRelate)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' #################################################################
#' ## Load the information about the profile
#' #################################################################
#' data(demoPedigreeEx1)
#' head(demoPedigreeEx1)
#'
#' #################################################################
#' ## The 1KG GDS file and the 1KG SNV Annotation GDS file
#' ## need to be located in the same directory
#' ## Note that the 1KG GDS file used for this example is a
#' ## simplified version and CANNOT be used for any real analysis
#' #################################################################
#' path1KG <- file.path(dataDir, "tests")
#'
#' fileReferenceGDS  <- file.path(path1KG, "ex1_good_small_1KG.gds")
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#'
#' #################################################################
#' ## The Sample SNP pileup files (one per sample) need
#' ## to be located in the same directory.
#' #################################################################
#' pathGeno <- file.path(dataDir, "example", "snpPileup")
#'
#' #################################################################
#' ## The path where the Profile GDS Files (one per sample)
#' ## will be created need to be specified.
#' #################################################################
#' pathProfileGDS <- file.path(tempdir(), "out.tmp")
#'
#' pathOut <- file.path(tempdir(), "res.out")
#'
#' #################################################################
#' ## A data frame containing general information about the study
#' ## is also required. The data frame must have
#' ## those 3 columns: "studyID", "study.desc", "study.platform"
#' #################################################################
#' studyDF <- data.frame(study.id="MYDATA",
#'                         study.desc="Description",
#'                         study.platform="PLATFORM",
#'                         stringsAsFactors=FALSE)
#'
#' ####################################################################
#' ## Fix seed to ensure reproducible results
#' ####################################################################
#' set.seed(3043)
#'
#' gds1KG <- snpgdsOpen(fileReferenceGDS)
#' dataRef <- select1KGPop(gds1KG, nbProfiles=2L)
#' closefn.gds(gds1KG)
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#'
#'     \dontrun{
#'
#'         RAIDS:::runWrapperAncestry(pedStudy=demoPedigreeEx1, studyDF=studyDF,
#'             pathProfileGDS=pathProfileGDS,
#'             pathGeno=pathGeno, pathOut=pathOut,
#'             fileReferenceGDS=fileReferenceGDS,
#'             fileReferenceAnnotGDS=fileAnnotGDS,
#'             chrInfo=chrInfo, syntheticRefDF=dataRef,
#'             studyType="DNA", genoSource="snp-pileup")
#'
#'         unlink(pathProfileGDS, recursive=TRUE, force=TRUE)
#'         unlink(pathOut, recursive=TRUE, force=TRUE)
#'
#'     }
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @keywords internal
runWrapperAncestry <- function(pedStudy, studyDF, pathProfileGDS,
                pathGeno, pathOut, fileReferenceGDS, fileReferenceAnnotGDS,
                chrInfo, syntheticRefDF,
                genoSource=c("snp-pileup", "generic", "VCF"),
                studyType=c("DNA", "RNA"), np=1L, blockTypeID=NULL,
                verbose=FALSE) {

    genoSource <- arg_match(genoSource)

    listProfiles <- pedStudy[, "Name.ID"]

    createStudy2GDS1KG(pathGeno=pathGeno, pedStudy=pedStudy,
        fileNameGDS=fileReferenceGDS, listProfiles=listProfiles,
        studyDF=studyDF, pathProfileGDS=pathProfileGDS, genoSource=genoSource,
        verbose=verbose)

    ## Open the 1KG GDS file (demo version)
    gdsReference <- snpgdsOpen(fileReferenceGDS)
    ## Open the 1KG GDS file and 1KG SNV Annotation file
    gdsRefAnnot <- openfn.gds(fileReferenceAnnotGDS)

    listProfileRef <- syntheticRefDF$sample.id
    studyDFSyn <- data.frame(study.id=paste0(studyDF$study.id, ".Synthetic"),
        study.desc=paste0(studyDF$study.id, " synthetic data"),
        study.platform=studyDF$study.platform, stringsAsFactors=FALSE)

    apply(pedStudy[,"Name.ID", drop=FALSE],1,FUN=function(x,gdsReference,
                        gdsRefAnnot, studyDF,pathProfileGDS,
                        pathOut, chrInfo, syntheticRefDF, studyDFSyn,
                        listProfileRef, studyType, verbose) {
        runProfileAncestry(gdsReference, gdsRefAnnot, studyDF,
            currentProfile=x, pathProfileGDS, pathOut, chrInfo,
            syntheticRefDF, studyDFSyn, listProfileRef,
            studyType, np=np, blockTypeID=blockTypeID, verbose=verbose)
        return(NULL)
    }, gdsReference=gdsReference, gdsRefAnnot=gdsRefAnnot,
        studyDF=studyDF, pathProfileGDS=pathProfileGDS, pathOut=pathOut,
        chrInfo=chrInfo, syntheticRefDF=syntheticRefDF,
        listProfileRef=listProfileRef,
        studyDFSyn=studyDFSyn, studyType=studyType, verbose=verbose)

    ## Close all GDS files
    closefn.gds(gdsReference)
    closefn.gds(gdsRefAnnot)

    ## Successful
    return(0L)
}

#' @title Create the Profile GDS file(s) for one or multiple specific profiles
#' using the information from a RDS Sample description file and the 1KG
#' GDS file
#'
#' @description The function uses the information for the Reference GDS file
#' and the RDS Sample Description file to create the Profile GDS file. One
#' Profile GDS file is created per profile. One Profile GDS file will be
#' created for each entry present in the \code{listProfiles} parameter.
#'
#' @param profileFile a \code{character} string representing the path to the
#' file: with genotype and the allele information of the profile A profile would have an
#' associated file called
#' if genoSource is "VCF", then "*vcf.gz",
#' if genoSource is "generic", then "*.txt.gz"
#' if genoSource is "snp-pileup", then "*.txt.gz".
#' if genoSource is "bam", then "*.bam" and "*.bai".
#'
#' @param profileName a \code{character} string representing the the profile Name.ID
#'
#' @param filePedRDS a \code{character} string representing the path to the
#' RDS file that contains the information about the sample to analyse.
#' The RDS file must
#' include a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings. The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined.
#'
#' @param fileNameGDS a \code{character} string representing the file name of
#' the Reference GDS file. The file must exist.
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
#' Default: \code{1}.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param listProfiles a \code{vector} of \code{character} string corresponding
#' to the profile identifiers that will have a Profile GDS file created. The
#' profile identifiers must be present in the "Name.ID" column of the Profile
#' RDS file passed to the \code{filePedRDS} parameter.
#' If \code{NULL}, all profiles present in the \code{filePedRDS} are selected.
#' Default: \code{NULL}.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the Profile GDS files will be created.
#' Default: \code{NULL}.
#'
#' @param genoSource a \code{character} string with two possible values:
#' 'snp-pileup', 'generic' or 'VCF'. It specifies if the genotype files
#' are generated by snp-pileup (Facets) or are a generic format CSV file
#' with at least those columns:
#' 'Chromosome', 'Position', 'Ref', 'Alt', 'Count', 'File1R' and 'File1A'.
#' The 'Count' is the depth at the specified position;
#' 'FileR' is the depth of the reference allele and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally the file can be a VCF file with at least those genotype
#' fields: GT, AD, DP.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{FALSE}.
#'
#' @return The function returns \code{0L} when successful.
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
#' ## Create the Profile GDS File for samples in 'listSamples' vector
#' ## (in this case, samples "ex1")
#' ## The Profile GDS file is created in the pathProfileGDS directory
#' result <- RAIDS:::createProfile(profileFile=file.path(dataDir, "ex1.txt.gz"),
#'             profileName="ex1",
#'             pedStudy=samplePED, fileNameGDS=fileGDS,
#'             studyDF=studyDF, listProfiles=c("ex1"),
#'             pathProfileGDS=tempdir(),
#'             genoSource="snp-pileup",
#'             verbose=FALSE)
#'
#' ## The function returns OL when successful
#' result
#'
#' ## The Profile GDS file 'ex1.gds' has been created in the
#' ## specified directory
#' list.files(tempdir())
#'
#' ## Remove Profile GDS file (created for demo purpose)
#' unlink(file.path(tempdir(), "ex1.gds"), force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt createfn.gds put.attr.gdsn closefn.gds read.gdsn
#' @importFrom S4Vectors isSingleNumber
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @keywords internal

createProfile <- function(profileFile, profileName,
                          filePedRDS=NULL, pedStudy=NULL, fileNameGDS,
                          batch=1, studyDF, listProfiles=NULL,
                          pathProfileGDS=NULL,
                          genoSource=c("snp-pileup", "generic", "VCF", "bam"),
                          paramProfile=list(ScanBamParam=NULL,
                                            PileupParam=NULL,
                                            yieldSize=10000000),
                          verbose=FALSE) {

    ## When filePedRDS is defined and pedStudy is null
    if (!(is.null(filePedRDS)) && is.null(pedStudy)) {
        ## The filePedRDS must be a character string and the file must exists
        if (!(is.character(filePedRDS) && (file.exists(filePedRDS)))) {
            stop("The \'filePedRDS\' must be a character string representing",
                 " the RDS Sample information file. The file must exist.")
        }
        ## Open the RDS Sample information file
        pedStudy <- readRDS(file=filePedRDS)
    } else if (!(is.null(filePedRDS) || is.null(pedStudy))) {
        stop("Both \'filePedRDS\' and \'pedStudy\' parameters cannot be ",
             "defined at the same time.")
    } else if (is.null(filePedRDS) && is.null(pedStudy)) {
        stop("One of the parameter \'fineNamePED\' of \'pedStudy\' must ",
             "be defined.")
    }

    ## Validate input parameters
    validatecreateProfile( pedStudy=pedStudy,
                               fileNameGDS=fileNameGDS, batch=batch, studyDF=studyDF,
                               listProfiles=listProfiles, pathProfileGDS=pathProfileGDS,
                               genoSource=genoSource, verbose=verbose)

    genoSource <- arg_match(genoSource)

    ## Read the Reference GDS file
    gdsReference <- snpgdsOpen(filename=fileNameGDS)

    ## Extract the chromosome and position information for all SNPs in 1KG GDS

    listPos <- NULL
    if(genoSource == "bam"){
        alDf <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))
        alDf <- matrix(unlist(strsplit(alDf,"\\/")),nrow=2)
        listPos <- data.frame(chr = read.gdsn(index.gdsn(gdsReference, "snp.chromosome")),
                            start = read.gdsn(index.gdsn(gdsReference, "snp.position")),
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
        listPos <- data.frame(snp.chromosome=read.gdsn(index.gdsn(node=gdsReference, "snp.chromosome")),
                              snp.position=read.gdsn(index.gdsn(node=gdsReference, "snp.position")))
    }
    ## Create a data.frame containing the information

    # Need to reformat for bam in varDf

    if(verbose) {
        message("Start ", Sys.time())
        message("Sample info DONE ", Sys.time())
    }

    generateProfileGDS(profileFile=profileFile, profileName=profileName,
                       listPos=listPos, offset=-1, minCov=10,
                       minProb=0.999, seqError=0.001, dfPedProfile=pedStudy, batch=batch,
                       studyDF=studyDF, pathProfileGDS=pathProfileGDS,
                       genoSource=genoSource, paramProfileGDS=paramProfile, verbose=verbose)

    if(verbose) {
        message("Genotype DONE ", Sys.time())
    }

    ## Close 1KG GDS file
    closefn.gds(gdsReference)

    ## Return successful code
    return(0L)
}


#' @title Run most steps leading to the ancestry inference call
#' on a specific profile (RNA or DNA)
#'
#' @description This function runs most steps leading to the ancestry inference
#' call on a specific profile. First, the function creates the Profile GDS file
#' for the specific profile using the information from a RDS Sample
#' description file and the Population reference GDS file.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{filePedRDS} or \code{pedStudy}
#' can be defined.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factor).
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Profile files will be created.
#' Default: \code{NULL}.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the VCF output of SNP-pileup for each sample. The
#' SNP-pileup files must be compressed (gz files) and have the name identifiers
#' of the samples. A sample with "Name.ID" identifier would have an
#' associated file called
#' if genoSource is "VCF", then "Name.ID.vcf.gz",
#' if genoSource is "generic", then "Name.ID.generic.txt.gz"
#' if genoSource is "snp-pileup", then "Name.ID.txt.gz".
#'
#' @param fileReferenceGDS  a \code{character} string representing the file
#' name of the Reference GDS file. The file must exist.
#'
#' @param fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Reference GDS Annotation file. The file must exist.
#'
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section.
#'
#' @param syntheticRefDF a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The \code{data.frame} must have those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#'
#' @param studyType a \code{character} string representing the type of study.
#' The possible choices are: "DNA" and "RNA". The type of study affects the
#' way the estimation of the allelic fraction is done. Default: \code{"DNA"}.
#'
#' @param genoSource a \code{character} string with two possible values:
#' 'snp-pileup', 'generic' or 'VCF', "bam". It specifies if the genotype files
#' are generated by snp-pileup (Facets) or are a generic format CSV file
#' with at least those columns:
#' 'Chromosome', 'Position', 'Ref', 'Alt', 'Count', 'File1R' and 'File1A'.
#' The 'Count' is the depth at the specified position;
#' 'FileR' is the depth of the reference allele and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally the file can be a VCF file with at least those genotype
#' fields: GT, AD, DP.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param blockTypeID a \code{character} string corresponding to the block
#' type used to extract the block identifiers. The block type must be
#' present in the GDS Reference Annotation file.
#'
#' @param paramAncestry a \code{list} parameters ...
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{FALSE}.
#'
#' @return a \code{list} containing 4 entries:
#' \describe{
#' \item{\code{pcaSample}}{ a \code{list} containing the information related
#' to the eigenvectors. The \code{list} contains those 3 entries:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{eigenvector.ref}}{ a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{\code{eigenvector}}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current profile projected on the PCA from the
#' reference profiles.}
#' }
#' }
#' \item{\code{paraSample}}{ a \code{list} containing the results with
#' different \code{D} and \code{K} values that lead to optimal parameter
#' selection. The \code{list} contains those entries:
#' \describe{
#' \item{\code{dfPCA}}{ a \code{data.frame} containing statistical results
#' on all combined synthetic results done with a fixed value of \code{D} (the
#' number of dimensions). The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{median}}{ a \code{numeric} representing the median of the
#' minimum AUROC obtained (within super populations) for all combination of
#' the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{mad}}{ a \code{numeric} representing the MAD of the minimum
#' AUROC obtained (within super populations) for all combination of the fixed
#' \code{D} value and all tested \code{K} values. }
#' \item{\code{upQuartile}}{ a \code{numeric} representing the upper quartile
#' of the minimum AUROC obtained (within super populations) for all
#' combination of the fixed \code{D} value and all tested \code{K} values. }
#' \item{\code{k}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for a fixed \code{D} value. }
#' }
#' }
#' \item{\code{dfPop}}{ a \code{data.frame} containing statistical results on
#' all combined synthetic results done with different values of \code{D} (the
#' number of dimensions) and \code{K} (the number of neighbors).
#' The \code{data.frame} contains those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{AUROC.min}}{ a \code{numeric} representing the minimum accuracy
#' obtained by grouping all the synthetic results by super-populations, for
#' the specified values of \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing the accuracy obtained
#' by grouping all the synthetic results for the specified values of \code{D}
#' and \code{K}.}
#' \item{\code{Accu.CM}}{ a \code{numeric} representing the value of accuracy
#' of the confusion matrix obtained by grouping all the synthetic results for
#' the specified values of \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{dfAUROC}}{ a \code{data.frame} the summary of the results by
#' super-population. The \code{data.frame} contains
#' those columns:
#' \describe{
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions).}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors).}
#' \item{\code{Call}}{ a \code{character} string representing the
#' super-population.}
#' \item{\code{L}}{ a \code{numeric} representing the lower value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' \item{\code{AUROC}}{ a \code{numeric} representing  the AUROC obtained for the
#' fixed values of super-population, \code{D} and \code{K}.}
#' \item{\code{H}}{ a \code{numeric} representing the higher value of the 95%
#' confidence interval for the AUROC obtained for the fixed values of
#' super-population, \code{D} and \code{K}.}
#' }
#' }
#' \item{\code{D}}{ a \code{numeric} representing the optimal \code{D} value
#' (the number of dimensions) for the specific profile.}
#' \item{\code{K}}{ a \code{numeric} representing the optimal \code{K} value
#' (the number of neighbors) for the specific profile.}
#' \item{\code{listD}}{ a \code{numeric} representing the optimal \code{D}
#' values (the number of dimensions) for the specific profile. More than one
#' \code{D} is possible.}
#' }
#' }
#' \item{\code{KNNSample}}{  a \code{data.frame} containing the inferred ancestry
#' for different values of \code{K} and \code{D}. The \code{data.frame}
#' contains those columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' }
#' }
#' \item{\code{KNNSynthetic}}{  a \code{data.frame} containing the inferred ancestry
#' for each synthetic data for different values of \code{K} and \code{D}.
#' The \code{data.frame}
#' contains those columns: "sample.id", "D", "K", "infer.superPop", "ref.superPop"
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current synthetic data.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry. }
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry. }
#' \item{\code{infer.superPop}}{ a \code{character} string representing the inferred
#' ancestry for the specified \code{D} and \code{K} values.}
#' \item{\code{ref.superPop}}{ a \code{character} string representing the known
#' ancestry from the reference}
#' }
#' }
#' \item{\code{Ancestry}}{ a \code{data.frame} containing the inferred
#' ancestry for the current profile. The \code{data.frame} contains those
#' columns:
#' \describe{
#' \item{\code{sample.id}}{ a \code{character} string representing the unique
#' identifier of the current profile.}
#' \item{\code{D}}{ a \code{numeric} representing the value of \code{D} (the
#' number of dimensions) used to infer the ancestry.}
#' \item{\code{K}}{ a \code{numeric} representing the value of \code{K} (the
#' number of neighbors) used to infer the ancestry.}
#' \item{\code{SuperPop}}{ a \code{character} string representing the inferred
#' ancestry.}
#' }
#' }
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
#' ## Required library for GDS
#' library(SNPRelate)
#'
#' ## Path to the demo 1KG GDS file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' #################################################################
#' ## Load the information about the profile
#' #################################################################
#' data(demoPedigreeEx1)
#' head(demoPedigreeEx1)
#'
#' #################################################################
#' ## The 1KG GDS file and the 1KG SNV Annotation GDS file
#' ## need to be located in the same directory
#' ## Note that the 1KG GDS file used for this example is a
#' ## simplified version and CANNOT be used for any real analysis
#' #################################################################
#' path1KG <- file.path(dataDir, "tests")
#'
#' fileReferenceGDS  <- file.path(path1KG, "ex1_good_small_1KG.gds")
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#'
#' #################################################################
#' ## The Sample SNP pileup files (one per sample) need
#' ## to be located in the same directory.
#' #################################################################
#' pathGeno <- file.path(dataDir, "example", "snpPileup")
#'
#' #################################################################
#' ## The path where the Profile GDS Files (one per sample)
#' ## will be created need to be specified.
#' #################################################################
#' pathProfileGDS <- file.path(tempdir(), "out.tmp")
#'
#'
#' #################################################################
#' ## A data frame containing general information about the study
#' ## is also required. The data frame must have
#' ## those 3 columns: "studyID", "study.desc", "study.platform"
#' #################################################################
#' studyDF <- data.frame(study.id="MYDATA",
#'                         study.desc="Description",
#'                         study.platform="PLATFORM",
#'                         stringsAsFactors=FALSE)
#'
#' ####################################################################
#' ## Fix seed to ensure reproducible results
#' ####################################################################
#' set.seed(3043)
#'
#' gds1KG <- snpgdsOpen(fileReferenceGDS)
#' dataRef <- select1KGPop(gds1KG, nbProfiles=2L)
#' closefn.gds(gds1KG)
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'
#'     ## Chromosome length information
#'     ## chr23 is chrX, chr24 is chrY and chrM is 25
#'     chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#'
#'     \dontrun{
#'
#'         res <- RAIDS:::wrapperAncestry(pedStudy=demoPedigreeEx1, studyDF=studyDF,
#'             pathProfileGDS=pathProfileGDS,
#'             pathGeno=pathGeno,
#'             fileReferenceGDS=fileReferenceGDS,
#'             fileReferenceAnnotGDS=fileAnnotGDS,
#'             chrInfo=chrInfo, syntheticRefDF=dataRef,
#'             studyType="LD", genoSource="snp-pileup")
#'
#'         unlink(pathProfileGDS, recursive=TRUE, force=TRUE)
#'
#'     }
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv
#' @importFrom rlang arg_match
#' @encoding UTF-8
#' @keywords internal
wrapperAncestry <- function(pedStudy, studyDF, pathProfileGDS,
                            profileFile, fileReferenceGDS, fileReferenceAnnotGDS,
                            chrInfo, syntheticRefDF,
                            genoSource=c("snp-pileup", "generic", "VCF", "bam"),
                            studyType=c("LD", "GeneAware"), np=1L, blockTypeID=NULL,
                            paramAncestry=list(ScanBamParam=NULL,
                                               PileupParam=NULL,
                                               yieldSize=10000000),
                            verbose=FALSE) {

    if(genoSource == "bam") {
        message("Process from bam is a new feature;",
                    " if you have an issue, please let us know")
    }
    genoSource <- arg_match(genoSource)

    listProfiles <- pedStudy[, "Name.ID"]
    # createProfile <- function(profileFile, profileName,
    #                           filePedRDS=NULL, pedStudy=NULL, fileNameGDS,
    #                           batch=1, studyDF, listProfiles=NULL,
    #                           pathProfileGDS=NULL,
    #                           genoSource=c("snp-pileup", "generic", "VCF", "bam"),
    #                           paramProfile=list(ScanBamParam=NULL,
    #                                             PileupParam=NULL,
    #                                             yieldSize=5000000),
    #                           verbose=FALSE)
    if(is.character(listProfiles)){
        for(profileCur in listProfiles){
            if(file.exists(file.path(pathProfileGDS, paste0(profileCur, ".gds")))){
                stop(paste0("The gds file for ", profileCur, " already exist."))
            }
        }
    }
    if(file.exists(file.path(pathProfileGDS, paste0(pedStudy$Name.ID[1], ".gds")))){
        stop(paste0("The gds file for ", pedStudy$Name.ID[1], " already exist."))
    }
    createProfile(profileFile=profileFile, profileName=pedStudy$Name.ID[1],
                  pedStudy=pedStudy, fileNameGDS=fileReferenceGDS,
                  studyDF=studyDF, pathProfileGDS=pathProfileGDS,
                  genoSource=genoSource, paramProfile=paramAncestry,
                  verbose=verbose)
    # createStudy2GDS1KG(pathGeno=pathGeno, pedStudy=pedStudy,
    #                    fileNameGDS=fileReferenceGDS, listProfiles=listProfiles,
    #                    studyDF=studyDF, pathProfileGDS=pathProfileGDS, genoSource=genoSource,
    #                    verbose=verbose)



    ## Open the 1KG GDS file (demo version)
    gdsReference <- snpgdsOpen(fileReferenceGDS)
    ## Open the 1KG GDS file and 1KG SNV Annotation file
    gdsRefAnnot <- openfn.gds(fileReferenceAnnotGDS)

    listProfileRef <- syntheticRefDF$sample.id
    studyDFSyn <- data.frame(study.id=paste0(studyDF$study.id, ".Synthetic"),
                             study.desc=paste0(studyDF$study.id, " synthetic data"),
                             study.platform=studyDF$study.platform, stringsAsFactors=FALSE)

    res <- profileAncestry(gdsReference, gdsRefAnnot, studyDF,
                           currentProfile=pedStudy[1,"Name.ID"], pathProfileGDS, chrInfo,
                           syntheticRefDF, studyDFSyn, listProfileRef,
                           studyType, np=np, blockTypeID=blockTypeID, verbose=verbose)


    ## Close all GDS files
    closefn.gds(gdsReference)
    closefn.gds(gdsRefAnnot)

    ## Successful
    return(res)
}

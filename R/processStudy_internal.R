#' @title Validate input parameters for pruningSample() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{pruningSample}} function.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1 KG GDS file.
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
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validatePruningSample(gds=gds1KG, method="corr",
#'      currentProfile="TGCA_01", studyID="TCGA",
#'      listSNP=c("sr10103", "sr10202"), slideWindowMaxBP=1000L,
#'      thresholdLD=0.008, np=1L, verbose=TRUE, chr=1,
#'      superPopMinAF=0.002, keepPrunedGDS=TRUE, pathProfileGDS=data.dir,
#'      keepFile=FALSE, pathPrunedGDS=data.dir, outPrefix="test")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validatePruningSample <- function(gds, method, currentProfile, studyID,
                        listSNP, slideWindowMaxBP, thresholdLD, np, verbose,
                        chr, superPopMinAF, keepPrunedGDS, pathProfileGDS,
                        keepFile, pathPrunedGDS, outPrefix) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gds, name="gds")

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


#' @title Validate input parameters for pruningSample() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{pruningSample}} function.
#'
#' @param gdsSample an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, the
#' GDS Sample file.
#'
#' @param sampleRM a \code{vector} of \code{character} strings representing
#' the identifiers of the 1KG reference samples that should not be used to
#' create the reference PCA.
#'
#' @param spRef TODO
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
#' @param fieldPopIn1KG TODO. Default: \code{"superPop"}.
#'
#' @param fieldPopInfAnc TODO. Default: \code{"SuperPop"}.
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#'
#' @param pcaList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis.  If \code{NULL},
#' the value \code{seq(2,15,1)} is assigned.
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016).
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA} only
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
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
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validateComputePoolSyntheticAncestryGr(gdsSample=gdsSample,
#'      sampleRM="TGCA_01", spRef="TCGA",
#'      studyIDSyn="TCGA", np=1L, listCatPop=c("AFR", "EAS", "SAS"),
#'      fieldPopIn1KG="SuperPop",  fieldPopInfAnc="Pop", kList=seq_len(3),
#'      pcaList=seq_len(10), algorithm="exact", eigen.cnt=12L,
#'      missing.rate=0.02)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateComputePoolSyntheticAncestryGr <- function(gdsSample, sampleRM,
        spRef, studyIDSyn, np, listCatPop, fieldPopIn1KG,
        fieldPopInfAnc, kList, pcaList, algorithm, eigen.cnt, missing.rate) {

    ## The gdsSample must be objects of class "gds.class"
    validateGDSClass(gdsSample, "gdsSample")

    ## The parameter sampleRM must be a single positive integer
    if(!(is.character(sampleRM))) {
        stop("The \'sampleRM\' parameter must be a vector of character ",
                "strings.")
    }

    ## The parameter studyIDSyn must be a character string
    if(!(is.character(studyIDSyn))) {
        stop("The \'studyIDSyn\' parameter must be a character string.")
    }

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

    ## The parameter eigen.cnt must be a single integer
    if(!(isSingleNumber(eigen.cnt))) {
        stop("The \'eigen.cnt\' parameter must be a single integer.")
    }

    ## The parameter missing.rate must be a single numer [0,1] or NaN
    if(!(((isSingleNumber(missing.rate) && missing.rate >= 0.0 &&
                    missing.rate <= 1.0)) || is.nan(missing.rate)))  {
        stop("The \'missing.rate\' parameter must be a single positive ",
                "numeric between zero and one or NaN.")
    }

    return(0L)
}


#' @title Validate input parameters for estimateAllelicFraction() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{estimateAllelicFraction}} function.
#'
#' @param gds an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the GDS Sample file.
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
#' @param minProb a single \code{numeric} between 0 and 1 representing TODO.
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
#' (a GDS file), the1 1KG SNV Annotation GDS file.
#'  **This parameter is RNA specific.**
#'
#' @param block.id a \code{character} string corresponding to the block
#' identifier in \code{gdsRefAnnot}. **This parameter is RNA specific.**
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## Get chromosome length information
#' ## Information from BSgenome.Hsapiens.UCSC.hg38 package version 1.4.4
#' ## Order by chromosomes 1 to 25
#' ## chr23 is chrX, chr24 is chrY and chrM is 25
#' chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
#'     181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
#'     135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
#'     83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
#'     156040895L, 57227415L,  16569L)
#'
#' ## The validatiion should be successful
#' RAIDS:::validateEstimateAllelicFraction(gds=gds1KG, gdsSample=gdsSample,
#'     currentProfile="Sample01", studyID="Synthetic", chrInfo=chrInfo,
#'     studyType="DNA", minCov=10L, minProb=0.03, eProb=0.002, cutOffLOH=10,
#'     cutOffHomoScore=11, wAR=2, cutOffAR=10, gdsRefAnnot=gds1KG,
#'     block.id="1")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateEstimateAllelicFraction <- function(gds, gdsSample, currentProfile,
        studyID, chrInfo, studyType, minCov, minProb, eProb, cutOffLOH,
        cutOffHomoScore, wAR, cutOffAR, gdsRefAnnot, block.id) {

    ## The gds and gdsSample must be objects of class "gds.class"
    validateGDSClass(gds, "gds")
    validateGDSClass(gdsSample, "gdsSample")

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

    return(0L)
}


#' @title Validate input parameters for createStudy2GDS1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{createStudy2GDS1KG}} function.
#'
#' @param pedStudy a \code{data.frame} with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings (no factor). The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter. Only \code{fileNamePED} or \code{pedStudy}
#' can be defined.
#'
#' @param fileNameGDS a \code{character} string representing the file name of
#' the 1KG GDS file. The file must exist.
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
#' passed to the \code{fileNamePED} parameter.
#' If \code{NULL}, all profiles in the \code{fileNamePED} are selected.
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
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' gds1KG <- file.path(data.dir, "gds1KG.gds")
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
#' ## The validatiion should be successful
#' RAIDS:::validateCreateStudy2GDS1KG(pedStudy=ped, fileNameGDS=gds1KG,
#'             batch=1, studyDF=studyInfo,
#'             listProfiles=c("Sample_01", "Sample_02"),
#'             pathProfileGDS=data.dir, verbose=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateCreateStudy2GDS1KG <- function(pedStudy, fileNameGDS, batch, studyDF,
                                       listProfiles, pathProfileGDS, verbose) {

    ## The PED study must have the mandatory columns
    if (!(all(c("Name.ID", "Case.ID", "Sample.Type", "Diagnosis", "Source")
                        %in% colnames(pedStudy)))) {
        stop("The PED study data frame is incomplete. ",
                "One or more mandatory columns are missing.")
    }

    ## The fileNameGDS must be a character string and the file must exists
    if (!(is.character(fileNameGDS) && (file.exists(fileNameGDS)))) {
        stop("The \'fileNameGDS\' must be a character string representing ",
                "the GDS 1KG file. The file must exist.")
    }

    ## The batch must be a single numeric
    if(!(isSingleNumber(batch))) {
        stop("The \'batch\' must be a single integer.")
    }

    ## The listProfiles must be a vector of character string
    if (!(is.character(listProfiles) || is.null(listProfiles))) {
        stop("The \'listProfiles\' must be a vector ",
                "of character strings (1 entry or more) or NULL.")
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, "verbose")

    return(0L)
}


#' @title Validate input parameters for createStudy2GDS1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{createStudy2GDS1KG}} function.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1KG GDS file.
#'
#' @param gdsSample an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the GDS Sample file.
#'
#' @param listFiles TODO.
#'
#' @param sample.ana.id TODO
#'
#' @param spRef TODO
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
#' @param fieldPopIn1KG a \code{character} string representing the name of TODO
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
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA}
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
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
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' listFiles <- file.path(data.dir,  "listSNPIndexes_Demo.rds")
#'
#' ## The validatiion should be successful
#' RAIDS:::validateComputeAncestryFromSyntheticFile(gds=gds1KG,
#'     gdsSample=gdsSample, listFiles=listFiles, sample.ana.id="sample01",
#'     spRef=NULL, studyIDSyn="Synthetic", np=1L, listCatPop=c("AFR", "EUR"),
#'     fieldPopIn1KG="superpop", fieldPopInfAnc="Superpop", kList=c(2, 3, 4),
#'     pcaList=c(3, 4, 5), algorithm="exact", eigen.cnt=32L, missing.rate=0.2)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateComputeAncestryFromSyntheticFile <- function(gds, gdsSample,
                listFiles, sample.ana.id, spRef, studyIDSyn, np, listCatPop,
                fieldPopIn1KG, fieldPopInfAnc, kList, pcaList,
                algorithm, eigen.cnt, missing.rate) {

    ## The gds and gdsSample must be objects of class "gds.class"
    validateGDSClass(gds, "gds")
    validateGDSClass(gdsSample, "gdsSample")

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

    ## The missing.rate must be a positive numeric between zero and one or NaN
    if (!is.nan(missing.rate)) {
        if (!(isSingleNumber(missing.rate) && (missing.rate >= 0.0) &&
                    (missing.rate <= 1.0))) {
            stop("The \'missing.rate\' must be a single numeric positive ",
                        "value between 0 and 1 or NaN.")
        }
    }

    return(0L)
}


#' @title Validate input parameters for computePCARefSample() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{computePCARefSample}} function.
#'
#' @param gdsSample an object of class \link[gdsfmt]{gds.class},
#' a GDS Sample file.
#'
#' @param name.id a single \code{character} string representing the sample
#' identifier.
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
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; if \code{NaN}, no missing
#' threshold.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The GDS Sample (opened)
#' gdsSample <- openfn.gds(file.path(data.dir,
#'                     "GDS_Sample_with_study_demo.gds"), readonly=TRUE)
#'
#' ## The validatiion should be successful
#' RAIDS:::validateComputePCARefSample(gdsSample=gdsSample, name.id="HCC01",
#'     studyIDRef="1KG", np=1L, algorithm="exact", eigen.cnt=32L,
#'     missing.rate=0.02)
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gdsSample)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateComputePCARefSample <- function(gdsSample, name.id, studyIDRef,
                                            np, algorithm,
                                            eigen.cnt, missing.rate) {

    ## The gdsSample must be object of class "gds.class"
    validateGDSClass(gdsSample, "gdsSample")

    ## Validate that name.id is a string
    if(!(is.character(name.id) && length(name.id) == 1)) {
        stop("The \'name.id\' parameter must be a single character string.")
    }

    ## Validate that studyIDRef is a string
    if(!(is.character(studyIDRef) && length(name.id) == 1)) {
        stop("The \'studyIDRef\' parameter must be a character string.")
    }

    ## The parameter np must be a single positive integer
    if(!(isSingleNumber(np) && (np > 0))) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    ## Validate that algorithm is a string
    if(!(is.character(algorithm))) {
        stop("The \'algorithm\' parameter must be a character string.")
    }

    ## The parameter eigen.cnt must be a single integer
    if(!(isSingleNumber(eigen.cnt))) {
        stop("The \'eigen.cnt\' parameter must be a single integer.")
    }

    ## The missing.rate must be a positive numeric between zero and one or NaN
    if (!is.nan(missing.rate)) {
        if (!(isSingleNumber(missing.rate) && (missing.rate >= 0.0) &&
                (missing.rate <= 1.0))) {
            stop("The \'missing.rate\' must be a single numeric positive ",
                    "value between 0 and 1 or NaN.")
        }
    }

    return(0L)
}

#' @title Validate input parameters for appendStudy2GDS1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{appendStudy2GDS1KG}} function.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the output of SNP-pileup, a VCF Sample file, for
#' each sample.
#'
#' @param fileNamePED a \code{character} string representing the path to the
#' RDS file that contains the information about the sample to analyse.
#'
#' @param fileNameGDS a \code{character} string representing the file name of
#' the 1KG GDS file. The file must exist.
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "studyID", "study.desc", "study.platform". All columns
#' must be in \code{character} strings.
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to the sample identifiers that will have a GDS Sample file created. The
#' sample identifiers must be present in the "Name.ID" column of the RDS file
#' passed to the \code{fileNamePED} parameter.
#' If \code{NULL}, all samples in the \code{fileNamePED} are selected.
#'
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Sample files will be created.
#'
#' @param verbose a \code{logical} indicating if message information should be
#' printed.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' gds1KG <- file.path(data.dir, "1KG_Demo.gds")
#' ped <- file.path(data.dir, "unrelatedPatientsInfo_Demo.rds")
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "studyID", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyInfo <- data.frame(study.id="Pancreatic.WES",
#'                 study.desc="Pancreatic study",
#'                 study.platform="WES",
#'                 stringsAsFactors=FALSE)
#'
#' ## The validatiion should be successful
#' RAIDS:::validateAppendStudy2GDS1KG(pathGeno=data.dir,
#'     fileNamePED=ped, fileNameGDS=gds1KG,
#'     batch=1L, studyDF=studyInfo, listSamples=c("HC01", "HC02"),
#'     pathProfileGDS=data.dir, verbose=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateAppendStudy2GDS1KG <- function(pathGeno, fileNamePED, fileNameGDS,
                                batch, studyDF, listSamples, pathProfileGDS,
                                verbose) {

    ## The pathGeno must be a character string and the path must exists
    if (!(is.character(pathGeno) && (dir.exists(pathGeno)))) {
        stop("The \'pathGeno\' must be a character string representing ",
                "a path. The path must exist.")
    }

    ## The fileNamePED must be a character string and the file must exists
    if (!(is.character(fileNamePED) && (file.exists(fileNamePED)))) {
        stop("The \'fileNamePED\' must be a character string representing ",
                "the RDS Sample information file. The file must exist.")
    }

    ## The fileNameGDS must be a character string and the file must exists
    if (!(is.character(fileNameGDS) && (file.exists(fileNameGDS)))) {
        stop("The \'fileNameGDS\' must be a character string representing ",
                "the GDS 1KG file. The file must exist.")
    }

    ## The batch must be a single numeric
    if (!(isSingleNumber(batch))) {
        stop("The \'batch\' must be a single integer.")
    }

    if (!(is.data.frame(studyDF) && all(c("study.id", "study.desc",
                        "study.platform") %in% colnames(studyDF)))) {
        stop("The \'studyDF\' must be a data.frame and contain those 3 ",
            "columns: \'study.id\', \'study.desc\' and \'study.platform\'.")
    }

    ## The listSamples must be a vector of character string
    if (!(is.character(listSamples) || is.null(listSamples))) {
        stop("The \'listSamples\' must be a vector ",
                "of character strings (1 entry or more) or NULL.")
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, name="verbose")

    return(0L)
}

#' @title Validate input parameters for add1KG2SampleGDS() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{add1KG2SampleGDS}} function.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
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
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(data.dir, "gds1KG.gds"), readonly=TRUE)
#'
#' ## The validatiion should be successful
#' RAIDS:::validateAdd1KG2SampleGDS(gds=gds1KG,
#'     gdsProfileFile=file.path(data.dir, "GDS_Sample_with_study_demo.gds"),
#'     currentProfile="Sample01", studyID="Synthetic")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateAdd1KG2SampleGDS <- function(gds, gdsProfileFile, currentProfile,
                                        studyID) {

    ## The gds must be an object of class "gds.class"
    validateGDSClass(gds=gds, name="gds")

    ## The gdsProfileFile must be a character string and the file must exists
    if(!(is.character(gdsProfileFile) && (file.exists(gdsProfileFile)))) {
        stop("The \'gdsSampleFile\' must be a character string representing ",
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


#' @title Calculate Principal Component Analysis (PCA) on SNV genotype dataset
#'
#' @description The functions calculates the principal component analysis (PCA)
#' for a list of pruned SNVs present in a GDS Sample file. The
#' \link[SNPRelate]{snpgdsPCA} function is used to do the calculation.
#'
#' @param gdsSample an object of class \link[SNPRelate]{SNPGDSFileClass},
#' the GDS Sample file.
#'
#' @param sample.ref TODO
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
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @param missing.rate a \code{numeric} value representing the threshold
#' missing rate at with the SNVs are discarded; the SNVs are retained in the
#' \link[SNPRelate]{snpgdsPCA} function
#' with "<= missing.rate" only; if \code{NaN}, no missing threshold.
#' Default: \code{0.025}.
#'
#' @return a \code{list} containing 2 entries:
#' \itemize{
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
#' # TODO
#' gds <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA
#' @encoding UTF-8
#' @keywords internal
computePCARefRMMulti <- function(gdsSample, sample.ref, listRM, np=1L,
                                    algorithm="exact", eigen.cnt=32L,
                                    missing.rate=0.025) {

    sample.Unrel <- sample.ref[which(!(sample.ref %in% listRM))]

    listPCA <- list()

    listPCA[["pruned"]] <- read.gdsn(index.gdsn(gdsSample, "pruned.study"))

    ## Calculate Principal Component Analysis (PCA) on SNV genotype dataset
    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsobj=gdsSample,
                                            sample.id=sample.Unrel,
                                            snp.id=listPCA[["pruned"]],
                                            num.thread=np,
                                            missing.rate=missing.rate,
                                            algorithm=algorithm,
                                            eigen.cnt=eigen.cnt,
                                            verbose=TRUE)

    return(listPCA)
}


#' @title Deprecated
#'
#' @description Deprecated
#'
#' @param listEigenvector TODO see return of computePCAsynthetic
#'
#' @param sample.ref TODO
#'
#' @param study.annot a  \code{data.frame} with one entry from study.annot in
#' the gds
#'
#' @param spRef TODO
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  K parameter. The K parameter represents the
#' number of neighbors used in the K-nearest neighbors analysis. If
#' \code{NULL}, the value \code{seq_len(15)} is assigned.
#' Default: \code{seq_len(15)}.
#'
#' @param pcaList TODO array of the pca dimension possible values
#'
#' @return A \code{list} TODO with the sample.id and eigenvectors
#' and a table with KNN callfor different K and pca dimension.
#'
#' @examples
#'
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASampLoading snpgdsPCASampLoading
#' @importFrom class knn
#' @encoding UTF-8
#' @keywords internal
computeKNNSuperPoprSynthetic <- function(listEigenvector, sample.ref,
                                    study.annot, spRef, kList=seq_len(15),
                                    pcaList=seq(2, 15, 1)) {

    ## The number of rows in study.annot must be one.
    if(nrow(study.annot) != 1) {
        stop("Number of samples in study.annot not equal to 1\n")
    }

    if(is.null(kList)){
        kList <- seq_len(15) #c(seq_len(14), seq(15,100, by=5))
    }
    if(is.null(pcaList)){
        pcaList <- seq(2, 15, 1)
    }

    resMat <- data.frame(sample.id=rep(listEigenvector$sample.id,
                                            length(pcaList) * length(kList)),
                            D=rep(0,length(pcaList) * length(kList)),
                            K=rep(0,length(pcaList) * length(kList)),
                            SuperPop=character(length(pcaList) * length(kList)),
                            stringsAsFactors=FALSE)

    listSuperPop <- c("EAS", "EUR", "AFR", "AMR", "SAS")

    #curPCA <- listPCA.Samples[[sample.id[sample.pos]]]
    eigenvect <- rbind(listEigenvector$eigenvector.ref,
                            listEigenvector$eigenvector)

    rownames(eigenvect) <- c(sample.ref[which(sample.ref !=
                                            study.annot$case.id[1])],
                                        listEigenvector$sample.id)

    totR <- 1
    for(pcaD in pcaList) {
        for(kV in  seq_len(length(kList))) {
            dCur <- paste0("d", pcaD)
            kCur <- paste0("k", kList[kV])
            resMat[totR,c("D", "K")] <- c(pcaD, kList[kV])

            pcaND <- eigenvect[ ,seq_len(pcaD)]
            y_pred <- knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                        test=pcaND[rownames(eigenvect)[nrow(eigenvect)],,
                                                    drop=FALSE],
                    cl=factor(spRef[rownames(eigenvect)[-1*nrow(eigenvect)]],
                                levels=listSuperPop, labels=listSuperPop),
                        k=kList[kV],
                        prob=FALSE)

            resMat[totR, paste0("SuperPop")] <- listSuperPop[as.integer(y_pred)]

            totR <- totR + 1
        } # end k
    } # end pca Dim
    listKNN <- list(sample.id=listEigenvector$sample.id,
                        sample1Kg=study.annot$case.id[1],
                        sp=spRef[study.annot$case.id[1]], matKNN=resMat)

    return(listKNN)
}


#' @title Deprecated
#'
#' @description Deprecated
#'
#' @param listEigenvector TODO see return of computePCARefSample
#'
#' @param sample.ref TODO
#'
#' @param study.annot a  \code{data.frame} with one entry from study.annot in
#' the gds
#'
#' @param spRef TODO
#'
#' @param kList a \code{vector} of \code{integer} representing  the list of
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. If \code{NULL},
#' the value \code{seq_len(15)} is assigned.
#' Default: \code{seq_len(15)}.
#'
#' @param pcaList TODO array of the pca dimension possible values
#'
#' @return A \code{list} TODO with the sample.id and eigenvectors
#' and a table with KNN callfor different K and pca dimension.
#'
#' @examples
#'
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASampLoading snpgdsPCASampLoading
#' @importFrom class knn
#' @encoding UTF-8
#' @keywords internal
computeKNNSuperPopSample <- function(gdsSample, listEigenvector, name.id,
                                        spRef, studyIDRef="Ref.1KG",
                                        kList=seq_len(15), pcaList=2:15) {

    if(is.null(kList)) {
        kList <- seq_len(15)#c(seq_len(14), seq(15,100, by=5))
    }
    if(is.null(pcaList)) {
        pcaList <- 2:15
    }
    if(length(name.id) != 1) {
        stop("Number of sample in study.annot not equal to 1\n")
    }

    study.annot.all <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

    sample.ref <- study.annot.all[which(study.annot.all$study.id ==
                                            studyIDRef), "data.id"]

    resMat <- data.frame(sample.id=rep(listEigenvector$sample.id,
                                        length(pcaList) * length(kList)),
                    D=rep(0,length(pcaList) * length(kList)),
                    K=rep(0,length(pcaList) * length(kList)),
                    SuperPop=character(length(pcaList) * length(kList)),
                    stringsAsFactors=FALSE)

    listSuperPop <- c("EAS", "EUR", "AFR", "AMR", "SAS")

    eigenvect <- rbind(listEigenvector$eigenvector.ref,
                        listEigenvector$eigenvector)

    rownames(eigenvect) <- c(sample.ref, listEigenvector$sample.id)

    totR <- 1
    for(pcaD in pcaList) {
        for(kV in  seq_len(length(kList))) {
            dCur <- paste0("d", pcaD)
            kCur <- paste0("k", kList[kV])
            resMat[totR,c("D", "K")] <- c(pcaD, kList[kV])

            pcaND <- eigenvect[ ,seq_len(pcaD)]
            y_pred <- knn(train=pcaND[rownames(eigenvect)[-1*nrow(eigenvect)],],
                test=pcaND[rownames(eigenvect)[nrow(eigenvect)],, drop=FALSE],
                cl=factor(spRef[rownames(eigenvect)[-1*nrow(eigenvect)]],
                                levels=listSuperPop, labels=listSuperPop),
                k=kList[kV],prob=FALSE)

            resMat[totR, paste0("SuperPop")] <- listSuperPop[as.integer(y_pred)]

            totR <- totR + 1
        } # end k
    } # end pca Dim
    listKNN <- list(sample.id=listEigenvector$sample.id, matKNN=resMat)

    return(listKNN)
}


#' @title Deprecated - Project patients onto existing principal component
#' axes (PCA)
#'
#' @description This function calculates the patient eigenvectors using
#' the specified SNP loadings. Deprecated
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param pathProfileGDS the path of an object of class \code{gds} related to
#' the sample
#'
#' @param listSamples a \code{vector} of string representing the samples for
#' which compute the PCA.
#'
#' @param np a single positive \code{integer} representing the number of
#' threads. Default: \code{1L}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @details
#'
#' More information about the method used to calculate the patient eigenvectors
#' can be found at the Bioconductor SNPRelate website:
#' https://bioconductor.org/packages/SNPRelate/
#'
#' @examples
#'
#' ## TODO
#' gds <- "TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom SNPRelate snpgdsPCASampLoading
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
computePCAForSamples <- function(gds, pathProfileGDS, listSamples, np=1L) {

    ## Validate that np is a single positive integer
    if(! (isSingleNumber(np) && np > 0)) {
        stop("The \'np\' parameter must be a single positive integer.")
    }

    for(i in seq_len(length(listSamples)) ){

        gdsSample <- openfn.gds(file.path(pathProfileGDS,
                                            paste0(listSamples[i], ".gds")))
        study.annot <- read.gdsn(index.gdsn(gdsSample, "study.annot"))

        if(length(which(study.annot$study.id == "Ref.1KG")) == 0) {
            stop("The study Ref.1KG is not define you must run the ",
                    "function addStudy1Kg \n")
        }

        sample.Unrel.All <- study.annot$data.id[study.annot$study.id ==
                                                        "Ref.1KG"]
        #sample.ref <- sample.Unrel.All$data.id
        listPCA <- computePrunedPCARef(gdsSample, sample.Unrel.All, np)

        listPCA[["samp.load"]] <- projectSample2PCA(gdsSample, listPCA,
                                                        listSamples[i], np)
        closefn.gds(gdsSample)

        saveRDS(listPCA, file.path(pathProfileGDS, paste0(listSamples[i],
                                            ".pca.pruned.rds")))

    }

    return(0L)
}


#' @title Deprecated Function
#'
#' @description TODO Deprecated
#'
#' @param gdsSample an object of class \code{gds} opened related to
#' the sample
#'
#' @param pruned TODO
#'
#' @param sample.id TODO
#'
#' @param sample.ref TODO
#'
#' @param study.annot a  \code{data.frame} with one entry from study.annot in
#' the gds
#'
#' @param algorithm a \code{character} string representing the algorithm used
#' to calculate the PCA. The 2 choices are "exact" (traditional exact
#' calculation) and "randomized" (fast PCA with randomized algorithm
#' introduced in Galinsky et al. 2016). Default: \code{"exact"}.
#'
#' @param eigen.cnt a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigen.cnt' <= 0, then all eigenvectors are returned.
#' Default: \code{32L}.
#'
#' @return A \code{list} TODO with the sample.id and eigenvectors.
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
#' # TODO
#' gds <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn
#' @importFrom SNPRelate snpgdsPCA snpgdsPCASampLoading snpgdsPCASampLoading
#' @encoding UTF-8
#' @keywords internal
computePCAsynthetic <- function(gdsSample, pruned, sample.id,
                                    sample.ref, study.annot,
                                    algorithm="exact", eigen.cnt=32L) {

    if(nrow(study.annot) != 1) {
        stop("Number of sample in study.annot not equal to 1\n")
    }

    sample.pos <- which(sample.id == study.annot$data.id[1])
    sample.Unrel <- sample.ref[which(sample.ref != study.annot$case.id[1])]

    g <- read.gdsn(index.gdsn(gdsSample, "genotype"),
                        start=c(1, sample.pos), count=c(-1, 1))

    listPCA <- list()

    listPCA[["pruned"]] <- pruned[which(g != 3)]
    rm(g)

    listPCA[["pca.unrel"]] <- snpgdsPCA(gdsSample,
                                            sample.id=sample.Unrel,
                                            snp.id=listPCA[["pruned"]],
                                            num.thread=1,
                                            algorithm=algorithm,
                                            eigen.cnt=eigen.cnt,
                                            verbose=TRUE)

    listPCA[["snp.load"]] <- snpgdsPCASNPLoading(listPCA[["pca.unrel"]],
                                                    gdsobj=gdsSample,
                                                    num.thread=1,
                                                    verbose=TRUE)

    listPCA[["samp.load"]] <- snpgdsPCASampLoading(listPCA[["snp.load"]],
                                            gdsobj=gdsSample,
                                            sample.id=sample.id[sample.pos],
                                            num.thread=1, verbose=TRUE)

    listRes <- list(sample.id=sample.id[sample.pos],
                        eigenvector.ref=listPCA[["pca.unrel"]]$eigenvect,
                        eigenvector=listPCA[["samp.load"]]$eigenvect)
    return(listRes)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param matKNN.All TODO see it is rbind matKNN of the
#' computeKNNSuperPoprSynthetic return from group synthetic data
#'
#' @param pedCall TODO see return of prepPedSynthetic1KG
#'
#' @param refCall TODO column name in pedCall with the call
#'
#' @param predCall a \code{character} string representing the name of
#' the column that will contain the inferred ancestry for the specified
#' dataset.
#'
#' @param listCall TODO array of the possible call
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
#' \itemize{
#' \item{dfPCA} {TODO}
#' \item{dfPop} {TODO}
#' \item{D} {TODO}
#' \item{K} {TODO}
#' \item{listD} {TODO}
#' }
#'
#' @examples
#'
#' # TODO
#' listEigenvector <- "TOTO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom stats mad median quantile
#' @encoding UTF-8
#' @keywords internal
selParaPCAUpQuartile <- function(matKNN.All, pedCall, refCall,
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
        matKNNCurD <- matKNN.All[which(matKNN.All$D == D ), ]
        listTMP <- list()
        listTMP.AUROC <- list()
        j <- 1
        ## Loop on all k neighbor values
        for (K in kList) {
            matKNNCur <- matKNNCurD[which(matKNNCurD$K == K), ]
            res <- computeSyntheticConfMat(matKNN=matKNNCur, pedCall=pedCall,
                        refCall=refCall, predCall=predCall, listCall=listCall)
            resROC <- computeSyntheticROC(matKNNCur, pedCall, refCall,
                                          predCall, listCall)

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
                            upQuartile=quantile(df[df$K %in% kList,
                                                   "AUROC.min"], 0.75),
                            K=kV)
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

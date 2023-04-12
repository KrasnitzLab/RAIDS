### Unit tests for gdsWrapper_internal.R functions

library(RAIDS)
library(withr)
library(testthat)
library(SNPRelate)

remove_local_GDS_file <- function(path) {
    unlink(x=path, force=TRUE)
}

local_GDS_file <- function(path) {
    GDS_file_tmp  <- createfn.gds(filename=path)
    defer_parent(remove_local_GDS_file(path=path))

    add.gdsn(GDS_file_tmp, "Ref.count", rep(10L, 12))
    add.gdsn(GDS_file_tmp, "Alt.count", rep(12L, 12))
    add.gdsn(GDS_file_tmp, "Total.count", rep(22L, 12))
    add.gdsn(GDS_file_tmp, "lap", rep(0.5, 12))
    sync.gds(GDS_file_tmp)

    return(GDS_file_tmp)
}


#############################################################################
### Tests runIBDKING() results
#############################################################################

context("runIBDKING() results")


test_that("runIBDKING() must return expected results", {

    ## Open genotype GDS file
    genoFile <- snpgdsOpen(snpgdsExampleFileName())
    withr::defer((snpgdsClose(genoFile)), envir = parent.frame())

    samples <- read.gdsn(index.gdsn(genoFile, "sample.id"))
    CEU <- samples[
            read.gdsn(index.gdsn(genoFile, "sample.annot/pop.group"))=="CEU"]
    CEU <- CEU[1:3]

    result <- RAIDS:::runIBDKING(gds=genoFile, profileID=CEU, snpID=NULL,
                                        maf=0.05, verbose=FALSE)

    sampleIDs <- c("NA07034", "NA07055", "NA12814")
    ibs0 <- matrix(c(0.000000000000000, 0.079705823891870, 0.088194306191519,
        0.079705823891870, 0.000000000000000, 0.094035785288270,
        0.088194306191519, 0.094035785288270, 0.000000000000000), nrow=3,
        byrow=TRUE)
    kinship <- matrix(c(0.500000000000000, 0.017237640936687, 0.001294777729823,
        0.017237640936687, 0.500000000000000, -0.014341590612777,
        0.001294777729823, -0.014341590612777, 0.500000000000000), nrow=3,
        byrow=TRUE)

    expect_equal(result$sample.id, sampleIDs)
    expect_equal(result$IBS0, ibs0)
    expect_equal(result$kinship, kinship)
})


#############################################################################
### Tests runLDPruning() results
#############################################################################

context("runLDPruning() results")


test_that("runLDPruning() must return expected results", {

    ## Open genotype GDS file
    genoFile <- snpgdsOpen(snpgdsExampleFileName())
    withr::defer((snpgdsClose(genoFile)), envir = parent.frame())

    set.seed(121)

    result <- RAIDS:::runLDPruning(gds=genoFile, verbose=FALSE)

    expect_equal(class(result), "list")
    expect_equal(length(result), 22)
    expect_equal(names(result) , paste0("chr", 1:22))
})


#############################################################################
### Tests appendGDSSample() results
#############################################################################

context("appendGDSSample() results")


test_that("appendGDSSample() must copy the expected entry in \"sample.annot\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(data.dir, "GDS_TEMP_05.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create sample.id field
    add.gdsn(node=GDS_file_tmp, name="sample.id", val=c("sample_01",
                                                            "sample_02"))

    dataInit <- data.frame(sex=c(1,1),  pop.group=c("ACB", "ACB"),
            superPop=c("AFR", "AFR"), batch=c(1, 1), stringsAsFactors=FALSE)

    add.gdsn(node=GDS_file_tmp, name="sample.annot", val=dataInit)
    sync.gds(gdsfile=GDS_file_tmp)

    ## Vector of SNV names
    samples <- c('sample_05', 'sample_08', 'sample_11')

    sample_info <- data.frame(Name.ID=samples, sex=c(1,2,1),
                                pop.group=c("ACB", "ACB", "ACB"),
                                superPop=c("AFR", "AFR", "AFR"),
                                stringsAsFactors=FALSE)
    rownames(sample_info) <- samples

    ## Add samples to the GDS file
    RAIDS:::appendGDSSample(gds=GDS_file_tmp,  pedDF=sample_info, batch=2,
        listSamples=c("sample_05", "sample_11"), verbose=FALSE)

    ## Read sample names from GDS file
    results1 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.id"))

    results2 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.annot"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected1 <- c("sample_01", "sample_02", samples[samples != "sample_08"])

    dataF <- sample_info[c("sample_05", "sample_11"),]
    dataF$Name.ID <- NULL
    dataF$batch <- rep(2, 2)
    expected2 <- rbind(dataInit, dataF)
    rownames(expected2) <- NULL

    expect_equal(results1, expected1)
    expect_equal(results2, expected2)
})


test_that("appendGDSSample() must print the expected message", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(data.dir, "GDS_TEMP_05.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create sample.id field
    add.gdsn(node=GDS_file_tmp, name="sample.id", val=c("sample_01",
                                                            "sample_02"))

    dataInit <- data.frame(sex=c(1,1),  pop.group=c("ACB", "ACB"),
            superPop=c("AFR", "AFR"), batch=c(1, 1), stringsAsFactors=FALSE)

    add.gdsn(node=GDS_file_tmp, name="sample.annot", val=dataInit)
    sync.gds(gdsfile=GDS_file_tmp)

    ## Vector of SNV names
    samples <- c('sample_05', 'sample_08', 'sample_11')

    sample_info <- data.frame(Name.ID=samples, sex=c(1,2,1),
                                pop.group=c("ACB", "ACB", "ACB"),
                                superPop=c("AFR", "AFR", "AFR"),
                                stringsAsFactors=FALSE)
    rownames(sample_info) <- samples

    message <- "Annot"

    ## Add samples to the GDS file
    expect_message(RAIDS:::appendGDSSample(gds=GDS_file_tmp,
        pedDF=sample_info, batch=2, listSamples=c("sample_05", "sample_11"),
        verbose=TRUE), regexp=message, all=TRUE, perl=TRUE)

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)
})


#############################################################################
### Tests addStudyGDSSample() results
#############################################################################

context("addStudyGDSSample() results")


test_that("addStudyGDSSample() must copy the expected entry in \"study.annot\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_11.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create PEF information
    pedInformation <- data.frame(Name.ID=c("sample_01", "sample_02", "sample_03"),
        Case.ID=c("sample_01", "sample_02", "sample_03"),
        Sample.Type=rep("Reference", 3), Diagnosis=rep("Reference", 3),
        Source=rep("IGSR", 3), stringsAsFactors=FALSE)

    rownames(pedInformation) <- pedInformation$Name.ID

    ## Create Study information
    studyInfo <- data.frame(study.id="Ref.1KG",
                        study.desc="Unrelated samples from 1000 Genomes",
                        study.platform="GRCh38 1000 Genotypes",
                        stringsAsFactors=FALSE)

    ## Add samples to the GDS file
    results3 <- RAIDS:::addStudyGDSSample(gds=GDS_file_tmp,  pedDF=pedInformation,
        batch=2, listSamples=c("sample_01", "sample_02"),
        studyDF=studyInfo, verbose=FALSE)

    ## Read sample names from GDS file
    results1 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="study.list"))

    results2 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="study.annot"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected3 <- c("sample_01", "sample_02")

    expected1 <- studyInfo

    expected2 <- data.frame(data.id=expected3, case.id=expected3,
        sample.type=rep("Reference", 2), diagnosis=rep("Reference", 2),
        source=rep("IGSR", 2), study.id=rep("Ref.1KG", 2),
        stringsAsFactors=FALSE)

    expect_equal(results1, expected1)
    expect_equal(results2, expected2)
    expect_equal(results3, expected3)
})


test_that("addStudyGDSSample() must copy the expected entry in \"study.annot\" node when the node already exists", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_12.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create Study information initial
    studyInfoInit <- data.frame(study.id="Ref.1KG",
            study.desc="Unrelated samples from 1000 Genomes",
            study.platform="GRCh38 1000 Genotypes",
            stringsAsFactors=FALSE)

    add.gdsn(GDS_file_tmp, "study.list", studyInfoInit)

    ## Create sample information initial
    sampleInfo <- data.frame(data.id=c("sample_01", "sample_02"), case.id=c("sample_01", "sample_02"),
        sample.type=rep("Reference", 2), diagnosis=rep("Reference", 2),
        source=rep("IGSR", 2), study.id=rep("Ref.1KG", 2),
        stringsAsFactors=FALSE)

    add.gdsn(GDS_file_tmp, "study.annot", sampleInfo)

    sync.gds(GDS_file_tmp)

    ## Create sample information novel
    pedInformation <- data.frame(Name.ID=c("sample_11", "sample_12"),
        Case.ID=c("sample_11", "sample_12"),
        Sample.Type=rep("Reference 2", 2), Diagnosis=rep("Reference 2", 2),
        Source=rep("IGSR", 2), stringsAsFactors=FALSE)

    ## Create study information novel
    studyInfo <- data.frame(study.id="Ref.1KG New",
                        study.desc="Unrelated samples from 1000 Genomes New",
                        study.platform="GRCh38 1000 Genotypes",
                        stringsAsFactors=FALSE)

    ## Add samples to the GDS file
    results3 <- RAIDS:::addStudyGDSSample(gds=GDS_file_tmp,  pedDF=pedInformation,
        batch=2, listSamples=NULL, studyDF=studyInfo, verbose=FALSE)

    ## Read sample names from GDS file
    results1 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="study.list"))

    results2 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="study.annot"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected3 <- c("sample_11", "sample_12")

    expected1 <- rbind(studyInfoInit, studyInfo)

    expected2 <- data.frame(data.id=c("sample_01", "sample_02", expected3),
        case.id=c("sample_01", "sample_02", expected3),
        sample.type=c(rep("Reference", 2), rep("Reference 2",2)),
        diagnosis=c(rep("Reference", 2), rep("Reference 2",2)),
        source=rep("IGSR", 4),
        study.id=c("Ref.1KG", "Ref.1KG", "Ref.1KG New", "Ref.1KG New"),
        stringsAsFactors=FALSE)

    expect_equal(results1, expected1)
    expect_equal(results2, expected2)
    expect_equal(results3, expected3)
})


test_that("addStudyGDSSample() must generate messages when verbose is TRUE", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_13.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create Study information initial
    studyInfoInit <- data.frame(study.id="Ref.1KG",
        study.desc="Unrelated samples from 1000 Genomes",
        study.platform="GRCh38 1000 Genotypes", stringsAsFactors=FALSE)

    add.gdsn(GDS_file_tmp, "study.list", studyInfoInit)

    ## Create sample information initial
    sampleInfo <- data.frame(data.id=c("sample_01", "sample_02"),
        case.id=c("sample_01", "sample_02"),
        sample.type=rep("Reference", 2), diagnosis=rep("Reference", 2),
        source=rep("IGSR", 2), study.id=rep("Ref.1KG", 2),
        stringsAsFactors=FALSE)

    add.gdsn(GDS_file_tmp, "study.annot", sampleInfo)

    sync.gds(GDS_file_tmp)

    ## Create sample information novel
    pedInformation <- data.frame(Name.ID=c("sample_11", "sample_12"),
        Case.ID=c("sample_11", "sample_12"),
        Sample.Type=rep("Reference 2", 2), Diagnosis=rep("Reference 2", 2),
        Source=rep("IGSR", 2), stringsAsFactors=FALSE)

    ## Create study information novel
    studyInfo <- data.frame(study.id="Ref.1KG New",
        study.desc="Unrelated samples from 1000 Genomes New",
        study.platform="GRCh38 1000 Genotypes", stringsAsFactors=FALSE)

    ## Add samples to the GDS file
    results3 <- RAIDS:::addStudyGDSSample(gds=GDS_file_tmp,
        pedDF=pedInformation, batch=2, listSamples=NULL,
        studyDF=studyInfo, verbose=FALSE)

    message <- "DONE"

    ## Add samples to the GDS file
    expect_message(RAIDS:::addStudyGDSSample(gds=GDS_file_tmp, pedDF=pedInformation,
        batch=2, listSamples=NULL, studyDF=studyInfo, verbose=TRUE),
        regexp=message, all=TRUE, perl=TRUE)

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)
})


test_that("addStudyGDSSample() must generate messages when verbose is TRUE", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_13.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create Study information initial
    studyInfoInit <- data.frame(study.id="Ref.1KG",
        study.desc="Unrelated samples from 1000 Genomes",
        study.platform="GRCh38 1000 Genotypes", stringsAsFactors=FALSE)

    add.gdsn(GDS_file_tmp, "study.list", studyInfoInit)

    ## Create sample information initial
    sampleInfo <- data.frame(data.id=c("sample_01", "sample_02"),
        case.id=c("sample_01", "sample_02"),
        sample.type=rep("Reference", 2), diagnosis=rep("Reference", 2),
        source=rep("IGSR", 2), study.id=rep("Ref.1KG", 2),
        stringsAsFactors=FALSE)

    add.gdsn(GDS_file_tmp, "study.annot", sampleInfo)

    sync.gds(GDS_file_tmp)

    ## Create sample information novel
    pedInformation <- data.frame(Name.ID=c("sample_11", "sample_12"),
        Case.ID=c("sample_11", "sample_12"),
        Sample.Type=rep("Reference 2", 2), Diagnosis=rep("Reference 2", 2),
        Source=rep("IGSR", 2), stringsAsFactors=FALSE)

    ## Create study information novel
    studyInfo <- data.frame(study.id="Ref.1KG New",
        study.desc="Unrelated samples from 1000 Genomes New",
        study.platform="GRCh38 1000 Genotypes", stringsAsFactors=FALSE)

    ## Add samples to the GDS file
    results3 <- RAIDS:::addStudyGDSSample(gds=GDS_file_tmp,  pedDF=pedInformation,
        batch=2, listSamples=NULL, studyDF=studyInfo, verbose=FALSE)

    message <- "DONE"

    ## Add samples to the GDS file
    expect_message(RAIDS:::addStudyGDSSample(gds=GDS_file_tmp,
        pedDF=pedInformation, batch=2, listSamples=NULL, studyDF=studyInfo,
        verbose=TRUE), regexp=message, all=TRUE, perl=TRUE)

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)
})


#############################################################################
### Tests addGDSStudyPruning() results
#############################################################################

context("addGDSStudyPruning() results")


test_that("addGDSStudyPruning() must copy the expected entry in \"pruned.study\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_03.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Vector of SNV names
    study <- c('s19771', 's19999', 's20122', 's21222')

    ## Add name of SNVs to the GDS file
    RAIDS:::addGDSStudyPruning(gds=GDS_file_tmp, pruned=study)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="pruned.study"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, study)
})


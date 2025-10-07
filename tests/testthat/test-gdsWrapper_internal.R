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
### Tests addGDSRef() results
#############################################################################

context("addGDSRef() results")


test_that("addGDSRef() must return expected result", {

    ## Create and open a temporary GDS file
    GDS_path  <- test_path("fixtures", "GDS_addGDSRef_Temp_01.gds")
    GDS_file_tmp  <- createfn.gds(filename=GDS_path)
    defer(unlink(x=GDS_path, force=TRUE), envir=parent.frame())

    ## Create "sample.id" node (the node must be present)
    sampleIDs <- c("HG00104", "HG00105", "HG00106", "HG00109", "HG00110")
    add.gdsn(node=GDS_file_tmp, name="sample.id", val=sampleIDs)
    sync.gds(GDS_file_tmp)

    listD <- list(rels=c("HG00106", "HG00110"), unrels=c("HG00104",
                                                    "HG00105", "HG00109"))

    RDS_file_tmp <- test_path("fixtures", "RDS_addGDSRef_Temp_01.RDS")

    saveRDS(listD, RDS_file_tmp)
    defer(unlink(RDS_file_tmp), envir=parent.frame())

    ## Add samples to the GDS file
    results3 <- RAIDS:::addGDSRef(gdsReference=GDS_file_tmp,
                                    filePart=RDS_file_tmp)

    ## Read sample names from GDS file
    results1 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.ref"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected1 <- c(1, 1, 0, 1, 0)

    expect_equal(results3, 0L)
    expect_equal(results1, expected1)
})


#############################################################################
### Tests appendGDSgenotypeMat() results
#############################################################################

context("appendGDSgenotypeMat() results")


test_that("appendGDSgenotypeMat() must copy the expected entry in \"genotype\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_06.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create a "genotype" node with initial matrix
    geno_initial <- matrix(rep(0L, 24), nrow=6)

    add.gdsn(node=GDS_file_tmp, name="genotype", val=geno_initial)
    sync.gds(GDS_file_tmp)

    ## new genotype to be added
    geno_new <- matrix(rep(1L, 12), nrow=6)

    ## Add genotype to the GDS file
    results0 <- RAIDS:::appendGDSgenotypeMat(gds=GDS_file_tmp,  matG=geno_new)

    ## Read genotype names from GDS file
    results1 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="genotype"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected1 <- cbind(geno_initial, geno_new)

    expect_equal(results1, expected1)
    expect_equal(results0, 0L)
})


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
    ibs0 <- matrix(c(0.000000000000000, 0.079928243970500, 0.088100458441300,
        0.079928243970500, 0.000000000000000, 0.093880805262109,
        0.088100458441300, 0.093880805262109, 0.000000000000000), nrow=3,
        byrow=TRUE)
    kinship <- matrix(c(0.500000000000000, 0.015481901439163, 0.001080847384349,
        0.015481901439163, 0.500000000000000, -0.014173571740078,
        0.001080847384349, -0.014173571740078, 0.500000000000000), nrow=3,
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
### Tests appendGDSRefSample() results
#############################################################################

context("appendGDSRefSample() results")


test_that("appendGDSRefSample() must copy the expected entry in \"sample.annot\" node of the GDS file", {

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
    RAIDS:::appendGDSRefSample(gdsReference=GDS_file_tmp,  dfPedReference=sample_info, batch=2,
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


test_that("appendGDSRefSample() must print the expected message", {

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
    expect_message(RAIDS:::appendGDSRefSample(gdsReference=GDS_file_tmp,
        dfPedReference=sample_info, batch=2, listSamples=c("sample_05", "sample_11"),
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
    results3 <- RAIDS:::addStudyGDSSample(gdsProfile=GDS_file_tmp,
        pedProfile=pedInformation,
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
    results3 <- RAIDS:::addStudyGDSSample(gdsProfile=GDS_file_tmp,  pedProfile=pedInformation,
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
    results3 <- RAIDS:::addStudyGDSSample(gdsProfile=GDS_file_tmp,
        pedProfile=pedInformation, batch=2, listSamples=NULL,
        studyDF=studyInfo, verbose=FALSE)

    message <- "DONE"

    ## Add samples to the GDS file
    expect_message(RAIDS:::addStudyGDSSample(gdsProfile=GDS_file_tmp,
        pedProfile=pedInformation, batch=2, listSamples=NULL,
        studyDF=studyInfo, verbose=TRUE),
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
    results3 <- RAIDS:::addStudyGDSSample(gdsProfile=GDS_file_tmp,
        pedProfile=pedInformation, batch=2, listSamples=NULL,
        studyDF=studyInfo, verbose=FALSE)

    message <- "DONE"

    ## Add samples to the GDS file
    expect_message(RAIDS:::addStudyGDSSample(gdsProfile=GDS_file_tmp,
        pedProfile=pedInformation, batch=2, listSamples=NULL, studyDF=studyInfo,
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
    RAIDS:::addGDSStudyPruning(gdsProfile=GDS_file_tmp, pruned=study)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="pruned.study"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, study)
})


#############################################################################
### Tests addUpdateLap() results
#############################################################################

context("addUpdateLap() results")


test_that("addUpdateLap() must copy the expected entry in \"lap\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_01.gds")

    ## Create and open a temporary GDS file
    tmpGDS  <- local_GDS_file(fileGDS)

    ## Vector of segment identifiers
    lap <- c(0.1, 0.22, 0.12, 0.43, 0.00, 0.00, 0.15, 0.43,
             0.00, 0.00, 0.15, 0.43)

    ## Add segments to the GDS file
    RAIDS:::addUpdateLap(gdsProfile=tmpGDS, snpLap=lap)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=tmpGDS, path="lap"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=tmpGDS)

    expect_equal(results, lap)
})


#############################################################################
### Tests getBlockIDs() results
#############################################################################


context("getBlockIDs() results")


test_that("getBlockIDs() must return the expected result", {

    ## Create a temporary GDS file in an test directory
    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "ex1_good_small_1KG_Annot_GDS.gds")

    annotFile <- openfn.gds(fileGDS)
    defer(closefn.gds(annotFile), envir=parent.frame())

    ## Vector of segment identifiers
    indexes <- c(1, 3, 6, 8, 9)

    ## Block identifiers for the selected SNVs
    result <- RAIDS:::getBlockIDs(gdsRefAnnot=annotFile, snpIndex=indexes,
                                blockTypeID="GeneS.Ensembl.Hsapiens.v86")

    expected <- rep(943, 5)

    expect_equal(result, expected)
})


test_that("getBlockIDs() must return expected error", {

    ## Create a temporary GDS file in an test directory
    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "ex1_good_small_1KG_Annot_GDS.gds")

    annotFile <- openfn.gds(fileGDS)
    defer(closefn.gds(annotFile), envir=parent.frame())

    ## Vector of segment identifiers
    indexes <- c(1, 3, 6, 8, 9)

    error_message <- paste0("The following block type is not found in the ",
                                "GDS Annotation file: \'InformationTremblay\'")

    expect_error(RAIDS:::getBlockIDs(gdsRefAnnot=annotFile, snpIndex=indexes,
                                        blockTypeID="InformationTremblay"),
                    error_message, fixed=TRUE)
})


#############################################################################
### Tests addUpdateSegment() results
#############################################################################

context("addUpdateSegment() results")


test_that("addUpdateSegment() must copy the expected entry in \"segment\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_01.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Vector of segment identifiers
    segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)

    ## Add segments to the GDS file
    RAIDS:::addUpdateSegment(gdsProfile=GDS_file_tmp, snpSeg=segments)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="segment"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, segments)
})


test_that("addUpdateSegment() must copy the expected entry in \"segment\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_02.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Vector of segment identifiers
    segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
    ## Vector of segment identifiers
    segments2 <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 5L, 5L)

    ## Add segments to the GDS file
    RAIDS:::addUpdateSegment(gdsProfile=GDS_file_tmp, snpSeg=segments)

    ## Update segments to the GDS file
    RAIDS:::addUpdateSegment(gdsProfile=GDS_file_tmp, snpSeg=segments2)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="segment"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, segments2)
})


#############################################################################
### Tests appendGDSSampleOnly() results
#############################################################################

context("appendGDSSampleOnly() results")


test_that("appendGDSSampleOnly() must copy the expected entry in \"sample.id\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_04.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create sample.id field
    add.gdsn(node=GDS_file_tmp, name="sample.id", val=c("sample_01",
                                                            "sample_02"))
    sync.gds(gdsfile=GDS_file_tmp)

    ## Vector of SNV names
    samples <- c('sample_05', 'sample_08', 'sample_11')

    ## Add name of samples to the GDS file
    RAIDS:::appendGDSSampleOnly(gds=GDS_file_tmp, listSamples=samples)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.id"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected <- c("sample_01", "sample_02", samples)

    expect_equal(results, expected)
})


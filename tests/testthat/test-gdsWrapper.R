### Unit tests for gdsWrapper functions

library(RAIDS)
library(withr)
library(testthat)


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
### Tests addUpdateSegment() results
#############################################################################

context("addUpdateSegment() results")


test_that("addUpdateSegment() must copy the expected entry in \"segment\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "GDS_TEMP_01.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

    ## Vector of segment identifiers
    segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)

    ## Add segments to the GDS file
    RAIDS:::addUpdateSegment(gds=GDS_file_tmp, snp.seg=segments)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="segment"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, segments)
})


test_that("addUpdateSegment() must copy the expected entry in \"segment\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "GDS_TEMP_02.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

    ## Vector of segment identifiers
    segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
    ## Vector of segment identifiers
    segments2 <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 5L, 5L)

    ## Add segments to the GDS file
    RAIDS:::addUpdateSegment(gds=GDS_file_tmp, snp.seg=segments)

    ## Update segments to the GDS file
    RAIDS:::addUpdateSegment(gds=GDS_file_tmp, snp.seg=segments2)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="segment"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, segments2)
})


#############################################################################
### Tests addUpdateLap() results
#############################################################################

context("addUpdateLap() results")


test_that("addUpdateLap() must copy the expected entry in \"lap\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "GDS_TEMP_01.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

    ## Vector of segment identifiers
    lap <- c(0.1, 0.22, 0.12, 0.43, 0.00, 0.00, 0.15, 0.43,
                0.00, 0.00, 0.15, 0.43)

    ## Add segments to the GDS file
    RAIDS:::addUpdateLap(gds=GDS_file_tmp, snp.lap=lap)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="lap"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, lap)
})



#############################################################################
### Tests addGDSStudyPruning() results
#############################################################################

context("addGDSStudyPruning() results")


test_that("addGDSStudyPruning() must copy the expected entry in \"pruned.study\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "GDS_TEMP_03.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

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


#############################################################################
### Tests appendGDSSampleOnly() results
#############################################################################

context("appendGDSSampleOnly() results")


test_that("appendGDSSampleOnly() must copy the expected entry in \"sample.id\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "GDS_TEMP_04.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

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


#############################################################################
### Tests appendGDSSample() results
#############################################################################

context("appendGDSSample() results")


test_that("appendGDSSample() must copy the expected entry in \"sample.annot\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "GDS_TEMP_05.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

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
    gdsFile <- file.path(data.dir, "GDS_TEMP_05.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

    ## Create sample.id field
    add.gdsn(node=GDS_file_tmp, name="sample.id", val=c("sample_01", "sample_02"))

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
    expect_message(RAIDS:::appendGDSSample(gds=GDS_file_tmp,  pedDF=sample_info, batch=2,
                listSamples=c("sample_05", "sample_11"), verbose=TRUE),
                regexp=message, all=TRUE, perl=TRUE)

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)
})



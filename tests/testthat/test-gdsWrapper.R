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
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_01.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Vector of segment identifiers
    segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)

    ## Add segments to the GDS file
    RAIDS:::addUpdateSegment(gdsProfile=GDS_file_tmp, snp.seg=segments)

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
    RAIDS:::addUpdateSegment(gdsProfile=GDS_file_tmp, snp.seg=segments)

    ## Update segments to the GDS file
    RAIDS:::addUpdateSegment(gdsProfile=GDS_file_tmp, snp.seg=segments2)

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
### Tests generateGDSRefSample() results
#############################################################################

context("generateGDSRefSample() results")


test_that("generateGDSRefSample() must copy the expected entry in \"sample.annot\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_06.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create PEF information
    pedInformation <- data.frame(sample.id=c("sample_01", "sample_02", "sample_03"),
                            Name.ID=c("sample_01", "sample_02", "sample_03"),
                            sex=c(1, 1, 2),  # 1:Male  2: Female
                            pop.group=c("ACB", "ACB", "ACB"),
                            superPop=c("AFR", "AFR", "AFR"),
                            batch=c(1, 1, 1), stringsAsFactors=FALSE)
    rownames(pedInformation) <- pedInformation$Name.ID

    ## Add samples to the GDS file
    results3 <- RAIDS:::generateGDSRefSample(gdsReference=GDS_file_tmp,  pedDF=pedInformation,
                                listSamples=c("sample_01", "sample_02"))

    ## Read sample names from GDS file
    results1 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.id"))

    results2 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.annot"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected1 <- c("sample_01", "sample_02")

    expected2 <- pedInformation[c("sample_01", "sample_02"),]
    expected2$Name.ID <- NULL
    expected2$sample.id <- NULL
    rownames(expected2) <- NULL

    expect_equal(results1, expected1)
    expect_equal(results2, expected2)
    expect_equal(results3, expected1)
})


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
    results3 <- RAIDS:::addGDSRef(gdsReference=GDS_file_tmp,  filePart=RDS_file_tmp)

    ## Read sample names from GDS file
    results1 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.ref"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected1 <- c(1, 1, 0, 1, 0)

    expect_equal(results3, 0L)
    expect_equal(results1, expected1)
})



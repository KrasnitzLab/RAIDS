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



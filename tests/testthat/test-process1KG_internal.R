### Unit tests for process1KG_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests validatePrepPed1KG() results
#############################################################################

context("validatePrepPed1KG() results")


test_that("validatePrepPed1KG() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    filePed <- file.path(dataDir, "PedigreeDemo.ped")

    result1 <- RAIDS:::validatePrepPed1KG(filePed=filePed, pathGeno=dataDir,
                                                batch=1L)

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateGenerateGDS1KG() results
#############################################################################

context("validateGenerateGDS1KG() results")


test_that("validateGenerateGDS1KG() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    filePed <- file.path(dataDir, "PedigreeDemo.ped")
    snpIndexFile <- file.path(dataDir, "listSNPIndexes_Demo.rds")
    fileSNV <- file.path(dataDir, "mapSNVSelected_Demo.rds")
    gdsFile <- file.path(dataDir, "1KG_Test.gds")

    result1 <- RAIDS:::validateGenerateGDS1KG(pathGeno=dataDir,
        filePedRDS=filePed, fileSNVIndex=snpIndexFile,
        fileSNVSelected=fileSNV, fileNameGDS=gdsFile, listSamples=NULL,
        verbose=FALSE)

    expect_identical(result1, 0L)
})

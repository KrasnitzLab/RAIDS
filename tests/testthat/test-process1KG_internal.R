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




#############################################################################
### Tests pruning1KGbyChr() results
#############################################################################

context("pruning1KGbyChr() results")


test_that("pruning1KGbyChr() must return expected results when all input are valid", {
    
    set.seed(121)
    
    dataDir <- test_path("fixtures")
    
    gds1KG <- snpgdsOpen(file.path(dataDir, "1KG_Test.gds"))
    withr::defer((snpgdsClose(gds1KG)), envir = parent.frame())
    
    outPrefix <- file.path(tempdir(), "Pruned_Test")
    result1 <- RAIDS:::pruning1KGbyChr(gdsReference=gds1KG, 
                                                    outPrefix=outPrefix)
    withr::defer(if(file.exists(paste0(outPrefix, ".rds"))) 
                    {unlink(paste0(outPrefix, ".rds"), force=TRUE)}, 
                    envir=parent.frame())
    
    expect_identical(result1, 0L)
    expect_true(file.exists(paste0(outPrefix, ".rds")))
    
    test <- readRDS(paste0(outPrefix, ".rds"))
    
    expect_identical(test, c("s5", "s7"))
})
### Unit tests for visualization_internal.R functions

library(RAIDS)
library(testthat)


#############################################################################
### Tests validateCreateAccuracyGraph() results
#############################################################################

context("validateCreateAccuracyGraph() results")


test_that("validateCreateAccuracyGraph() must return expected results when all input are valid", {
    
    dataDir <- system.file("extdata", package="RAIDS")
    fileRDS <- file.path(dataDir, "TEST_01.infoCall.RDS")
    
    result1 <- RAIDS:::validateCreateAccuracyGraph(fileRDS=fileRDS, 
        title="test", selectD=c(5, 15), selectColor=c("red","darkblue"))
    
    expect_identical(result1, 0L)
})

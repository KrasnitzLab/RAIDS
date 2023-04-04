### Unit tests for process1KG_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests validatePrepPed1KG() results
#############################################################################

context("validatePrepPed1KG() results")


test_that("validatePrepPed1KG() must return expected results when all input are valid", {

    data.dir <- test_path("fixtures")
    pedFile <- file.path(data.dir, "PedigreeDemo.ped")

    result1 <- RAIDS:::validatePrepPed1KG(pedFile=pedFile, PATHGENO=data.dir,
                                                batch=1L)

    expect_identical(result1, 0L)
})

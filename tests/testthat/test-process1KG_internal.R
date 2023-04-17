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
    filePed <- file.path(data.dir, "PedigreeDemo.ped")

    result1 <- RAIDS:::validatePrepPed1KG(filePed=filePed, pathGeno=data.dir,
                                                batch=1L)

    expect_identical(result1, 0L)
})

### Unit tests for allelicFraction_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests getTableSNV() results
#############################################################################

context("validateGetTableSNV() results")


test_that("validateGetTableSNV() must return epxected results when all input are valid", {

    data.dir <- test_path("fixtures")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validateGetTableSNV(gds=gdsF, gdsSample=gdsF,
                    sampleCurrent="TEST", study.id="TEST",
                    minCov=10L, minProb=0.999, eProb=0.001,verbose=TRUE)

    expect_identical(result1, 0L)
})

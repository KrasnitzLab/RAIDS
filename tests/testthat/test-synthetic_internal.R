### Unit tests for synthetic_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests validateSyntheticGeno() results
#############################################################################

context("validateSyntheticGeno() results")


test_that("validateSyntheticGeno() must return expected results when all input are valid", {

    data.dir <- test_path("fixtures")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validateSyntheticGeno(gds=gdsF, gdsRefAnnot=gdsF,
            gdsSampleFile=gdsFIle, data.id.profile="TCGA_001",
            listSampleRef=c("Sample01", "Sample02"), nbSim=1L,
            prefId="TEST", pRecomb=0.2, minProb=0.3, seqError=0.2)

    expect_identical(result1, 0L)
})

### Unit tests for synthetic functions

library(RAIDS)
library(withr)
library(testthat)




#############################################################################
### Tests select1KGPop() results
#############################################################################

context("select1KGPop() results")


test_that("select1KGPop() must return error when gds is a character string", {

    error_message <- paste0("The \'gds\' must be an object of class ",
                                "\'gdsn.class\' or \'gds.class\'")

    expect_error(select1KGPop(gds="test.gds", nbSamples=10L),
                    error_message)
})


test_that("select1KGPop() must return error when nbSamples is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'nbSamples\' parameter must be a ",
                                "single positive integer.")

    expect_error(select1KGPop(gds=gdsF, nbSamples="CANADA"),
                    error_message)
})


#############################################################################
### Tests splitSelectByPop() results
#############################################################################

context("splitSelectByPop() results")


test_that("splitSelectByPop() must return error when dataRef is a character string", {

    error_message <- "The \'dataRef\' must be a data.frame object."
    expect_error(splitSelectByPop(dataRef="GDS"), error_message)
})


test_that("splitSelectByPop() must return error when dataRef does not have the pop.group column", {

    demo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC", "SampleD"),
                            pop=c("TSI", "TSI", "YRI", "YRI"),
                            superPop=c("EUR", "EUR", "AFR", "AFR"))

    error_message <- "The \'dataRef\' must have a column named \'pop.group\'."
    expect_error(splitSelectByPop(dataRef=demo), error_message)
})


test_that("splitSelectByPop() must return error when dataRef does not have the sample.id column", {

    demo <- data.frame(sample=c("SampleA", "SampleB", "SampleC", "SampleD"),
                       pop.group=c("TSI", "TSI", "YRI", "YRI"),
                       superPop=c("EUR", "EUR", "AFR", "AFR"))

    error_message <- "The \'dataRef\' must have a column named \'sample.id\'."
    expect_error(splitSelectByPop(dataRef=demo), error_message)
})


test_that("splitSelectByPop() must return expected results", {

    demo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC", "SampleD",
                                    "SampleE", "SampleF"),
                       pop.group=c("TSI", "TSI", "YRI", "YRI", "LWK", "LWK"),
                       superPop=c("EUR", "EUR", "AFR", "AFR", "AFR", "AFR"))

    expected <- matrix(data=c("SampleA", "SampleB", "SampleC", "SampleD",
                              "SampleE", "SampleF"), byrow = FALSE, nrow = 2)
    colnames(expected) <- c("TSI", "YRI", "LWK")

    result <- splitSelectByPop(dataRef=demo)

    expect_identical(result, expected)
})


#############################################################################
### Tests syntheticGeno() results
#############################################################################

context("syntheticGeno() results")


test_that("syntheticGeno() must return error when seqError is a character string", {

    error_message <- paste0("The \'seqError\' parameter must be a single ",
                "positive numeric value between 0 and 1.")

    expect_error(syntheticGeno(gds="test.gds",
                              gdsRefAnnot="test.annot.gds",
                              gdsSampleFile="test.gds",
                              data.id.profile="test",
                              listSampleRef=c("Sample1", "sample2"),
                              nbSim=1, prefId="", pRecomb=0.01,
                              minProb=0.999, seqError="0.001"), error_message)
})


test_that("syntheticGeno() must return error when minProb is a character string", {

    error_message <- paste0("The \'minProb\' parameter must be a single ",
                            "positive numeric value between 0 and 1.")

    expect_error(syntheticGeno(gds="test.gds",
                               gdsRefAnnot="test.annot.gds",
                               gdsSampleFile="test.gds",
                               data.id.profile="test",
                               listSampleRef=c("Sample1", "sample2"),
                               nbSim=1, prefId="", pRecomb=0.01,
                               minProb="0.999", seqError=0.001), error_message)
})



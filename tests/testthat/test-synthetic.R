### Unit tests for synthetic functions

library(RAIDS)
library(withr)




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


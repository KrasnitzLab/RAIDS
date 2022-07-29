### Unit tests for processStudy.R functions

library(RAIDS)
library(withr)




#############################################################################
### Tests projectSample2PCA() results
#############################################################################

context("projectSample2PCA() results")


test_that("projectSample2PCA() must return error when np is character string", {

    data.dir <- system.file("extdata", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Demo.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gds=gdsFIle, listPCA=list(),
                    sample.current="sample1",
                    np="test"), error_message)
})

test_that("projectSample2PCA() must return error when np is negative integer", {

    data.dir <- system.file("extdata", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Demo.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gds=gdsFIle, listPCA=list(),
                                    sample.current="sample1",
                                    np=-1L), error_message)
})

test_that("projectSample2PCA() must return error when np is zero", {

    data.dir <- system.file("extdata", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Demo.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gds=gdsFIle, listPCA=list(),
                                    sample.current="sample1",
                                    np=0L), error_message)
})

test_that("projectSample2PCA() must return error when sample.current is number", {

    data.dir <- system.file("extdata", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Demo.gds")

    error_message <- paste0("The \'sample.current\' ",
                                "parameter must be a character string.")

    expect_error(projectSample2PCA(gds=gdsFIle, listPCA=list(),
                                        sample.current=101,
                                        np=1L), error_message)
})



context("appendStudy2GDS1KG() results")


test_that("appendStudy2GDS1KG() must return error when fileNameGDS is numeric", {

    error_message <- paste0("The \'fileNameGDS\' must be a character string representing ",
                                "the GDS study file. The file must exist.")

    expect_error(appendStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
        fileNamePED="test", fileNameGDS=33, batch=1,
        studyDF="test", listSamples=NULL, PATHSAMPLEGDS=NULL), error_message)
})

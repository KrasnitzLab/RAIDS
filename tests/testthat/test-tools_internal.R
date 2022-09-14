### Unit tests for tools_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests validateGDSClass() results
#############################################################################

context("validateGDSClass() results")


test_that("validateGDSClass() must return expected results when all input are valid", {

    data.dir <- test_path("fixtures")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validateGDSClass(gds=gdsF, name="gds")

    expect_identical(result1, 0L)
})


test_that("validateGDSClass() must return error when input is not valid", {

    expected <- "The \'todo_01\' must be an object of class \'gds.class\'."

    expect_error(RAIDS:::validateGDSClass(gds="toto.gds", name="todo_01"),
                    expected, fixed=TRUE)
})


#############################################################################
### Tests validateLogical() results
#############################################################################

context("validateLogical() results")


test_that("validateLogical() must return expected results when all input are valid", {

    result1 <- RAIDS:::validateLogical(logical=TRUE, name="parameter_01")

    expect_identical(result1, 0L)
})


test_that("validateLogical() must return error when input is not valid", {

    expected <- "The \'param_01\' parameter must be a logical (TRUE or FALSE)."

    expect_error(RAIDS:::validateLogical(logical="toto", name="param_01"),
                    expected, fixed=TRUE)
})


#############################################################################
### Tests validateSingleRatio() results
#############################################################################

context("validateSingleRatio() results")


test_that("validateSingleRatio() must return expected results when all input are valid", {

    result1 <- RAIDS:::validateSingleRatio(value=0.02, name="parameter_01")

    expect_identical(result1, 0L)
})


test_that("validateSingleRatio() must return error when input is not valid", {

    expected <- paste0("The \'param_01\' must be a single numeric positive ",
                        "value between 0 and 1.")

    expect_error(RAIDS:::validateSingleRatio(value="toto", name="param_01"),
                    expected, fixed=TRUE)
})


test_that("validateSingleRatio() must return error when input is negative value", {

    expected <- paste0("The \'param_02\' must be a single numeric positive ",
                        "value between 0 and 1.")

    expect_error(RAIDS:::validateSingleRatio(value=-0.01, name="param_02"),
                 expected, fixed=TRUE)
})


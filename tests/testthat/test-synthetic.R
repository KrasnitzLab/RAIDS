### Unit tests for synthetic functions

library(aicsPaper)
library(withr)




#############################################################################
### Tests select1KGPop() results
#############################################################################

context("select1KGPop() results")


test_that("select1KGPop() must return error when nbSamples is a character string", {

    error_message <- "The \'nbSamples\' parameter must be a single positive integer."

    expect_error(select1KGPop(gds = "test.gds", nbSamples = "CANADA"),
                    error_message)
})

### Unit tests for process1KG.R functions

library(aicsPaper)




#############################################################################
### Tests prepPed1KG() results
#############################################################################

context("prepPed1KG() results")


test_that("prepPed1KG() must return error when batch.v is a character string", {

    data.dir <- system.file("extdata", package="aicsPaper")

    pedDemoFile <- file.path(data.dir, "PedigreeDemo.txt")

    error_message <- "The batch.v must be an integer."

    expect_error(prepPed1KG(pedFile=pedFile, PATHGENO=data.dir,
                                batch.v="SAVE"), error_message)
})


test_that("prepPed1KG() must return error when batch.v is a float", {

    data.dir <- system.file("extdata", package="aicsPaper")

    pedDemoFile <- file.path(data.dir, "PedigreeDemo.txt")

    error_message <- "The batch.v must be an integer."

    expect_error(prepPed1KG(pedFile=pedFile, PATHGENO=data.dir,
                                batch.v=0.111), error_message)
})

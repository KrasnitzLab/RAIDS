### Unit tests for tools.R functions

library(RAIDS)
library(testthat)
library(withr)

#############################################################################
### Tests createAccuracyGraph() results
#############################################################################

context("createAccuracyGraph() results")


test_that("createAccuracyGraph() must return error when fileRDS is not a character string", {
    
    error_message <- "The \'fileRDS\' parameter must be a character string."
    
    expect_error(createAccuracyGraph(fileRDS=44, title="", 
        selectD=c(3,7,11), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
        error_message)
})

test_that("createAccuracyGraph() must return error when fileRDS does not exist", {
    
    error_message <- "The \'fileRDS\' file does not exist."
    
    expect_error(createAccuracyGraph(fileRDS="./toto.RDS", title="", 
        selectD=c(3,7,11), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
        error_message)
})

test_that("createAccuracyGraph() must return error when fileRDS is not RDS", {
    
    
    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    
    error_message <- "The \'fileRDS\' must have a RDS (or rds) extension."
    
    expect_error(createAccuracyGraph(fileRDS=fileGDS, title="", 
        selectD=c(3,7,11), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
        error_message, fixed=TRUE)
})


test_that("createAccuracyGraph() must return error when fileRDS is not RDS", {
    
    fileGDS <- test_path("fixtures", "TEST_01.infoCall.RDS")
    
    error_message <- "The \'selectD\' parameter cannot be empty."
    
    expect_error(createAccuracyGraph(fileRDS=fileGDS, title="", 
        selectD=c(), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
        error_message)
})


test_that("createAccuracyGraph() must return error when title is number", {
    
    fileGDS <- test_path("fixtures", "TEST_01.infoCall.RDS")
    
    error_message <- "The \'title\' must be a character string."
    
    expect_error(createAccuracyGraph(fileRDS=fileGDS, title=33, 
        selectD=c(3,5,6), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
                 error_message)
})


test_that("createAccuracyGraph() must return error when selectD is empty", {
    
    fileGDS <- test_path("fixtures", "TEST_01.infoCall.RDS")
    
    error_message <- "The \'selectD\' parameter cannot be empty."
    
    expect_error(createAccuracyGraph(fileRDS=fileGDS, title="", 
        selectD=c(), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
        error_message)
})


test_that("createAccuracyGraph() must return error when selectD has 6 elements", {
    
    fileGDS <- test_path("fixtures", "TEST_01.infoCall.RDS")
    
    error_message <- paste0("The \'selectD\' parameter can contain a ", 
                                "maximum of 5 elements.")
    
    expect_error(createAccuracyGraph(fileRDS=fileGDS, title="", 
        selectD=c(1,2,3,4,5,6), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
        error_message)
})


test_that("createAccuracyGraph() must return error when selectColor shorter than selectD", {
    
    fileGDS <- test_path("fixtures", "TEST_01.infoCall.RDS")
    
    error_message <- paste0("The \'selectColor\' parameter must be the ", 
        "same length than the \'selectD\' parameter.")
    
    expect_error(createAccuracyGraph(fileRDS=fileGDS, title="", 
        selectD=c(1,2,3,4), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
                 error_message)
})


test_that("createAccuracyGraph() must return a gglot object when successful", {
    
    fileGDS <- test_path("fixtures", "TEST_01.infoCall.RDS")
    
    graphE <-createAccuracyGraph(fileRDS=fileGDS, title="", 
                selectD=c(7,8,9), selectColor=c("#5e688a", "#cd5700", "#CC79A7"))
    
    testthat::expect_is(graphE, "ggplot")
})

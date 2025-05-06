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


test_that("createAccuracyGraph() must return error when selectD value not in the file", {
    
    fileGDS <- test_path("fixtures", "TEST_01.infoCall.RDS")
    
    error_message <- paste0("Not all values in \'selectD\' are present in the RDS file.")
    
    expect_error(createAccuracyGraph(fileRDS=fileGDS, title="", 
                    selectD=c(32,32,34), selectColor=c("#5e688a", "#cd5700", "#CC79A7")), 
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


#############################################################################
### Tests createAUROCGraph() results
#############################################################################

context("createAUROCGraph() results")


test_that("createAccuracyGraph() must return a gglot object when successful", {
    
    dfAuroc <- data.frame(D=c(rep(2, 15), rep(3, 15)),
                          K=c(rep(c(2, 3, 4), 10)), 
                          Call=c("EUR", "EUR", "EUR", "AMR", "AMR", "AMR",
                                 "EAS", "EAS", "EAS", "SAS", "SAS", "SAS", 
                                 "AFR", "AFR", "AFR",
                                 "EUR", "EUR", "EUR", "AMR", "AMR", "AMR",
                                 "EAS", "EAS", "EAS", "SAS", "SAS", "SAS", 
                                 "AFR", "AFR", "AFR"),
                          AUROC=c(0.95000, 0.9628737, 0.9701246, 0.8337130,
                                0.8509514, 0.9800000, 0.9158718, 0.9267399,
                                0.9386384, 0.7484138, 0.9000000, 0.9892067,
                                0.88000, 0.8758737, 0.9021246,
                                0.99000, 0.9888737, 0.9931246, 0.8837130,
                                0.8959514, 1.0000000, 0.9788718, 0.9977399,
                                0.9886384, 0.8244138, 1.0000000, 0.9982067,
                                0.92000, 0.8998737, 0.9251246))
    dfAuroc$L <- dfAuroc$AUROC - 0.03
    dfAuroc$H <- dfAuroc$AUROC + 0.02
    dfAuroc$H[which(dfAuroc$H > 1.0000)] <- 1.0000000
    
    
    graphE <-createAUROCGraph(dfAUROC=dfAuroc, title="", 
        selectD=c(2,3), selectColor=c("#5e688a", "#CC79A7"))
    
    testthat::expect_is(graphE, "ggplot")
})


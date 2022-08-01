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

test_that("appendStudy2GDS1KG() must return error when fileNamePED is numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study", study.platform="WES",
                        stringsAsFactors=FALSE)

    error_message <- paste0("The \'fileNamePED\' must be a character string representing ",
                        "the RDS Sample information file. The file must exist.")

    expect_error(appendStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
        fileNamePED=21, fileNameGDS=gdsFIle, batch=1,
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when fileNameGDS is numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- paste0("The \'fileNameGDS\' must be a character string representing ",
                                "the GDS 1KG file. The file must exist.")

    expect_error(appendStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
        fileNamePED=sampleRDS, fileNameGDS=33, batch=1,
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when batch is a vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES", study.desc="Pancreatic",
                                study.platform="WES", stringsAsFactors=FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(appendStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
        fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch=c(1,2),
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})

test_that("appendStudy2GDS1KG() must return error when batch is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(appendStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
        fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch="2",
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when listSamples is a numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study", study.platform="WES",
                        stringsAsFactors=FALSE)

    error_message <- paste0("The \'listSamples\' must be a vector ",
        "of character strings (1 entry or more) or NULL.")

    expect_error(appendStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
        fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch=2,
        studyDF=studyInfo, listSamples=33, PATHSAMPLEGDS=NULL, verbose=FALSE), error_message, fixed=TRUE)
})


test_that("appendStudy2GDS1KG() must return error when verbose is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- 'The \'verbose\' parameter must be a logical (TRUE or FALSE).'

    expect_error(appendStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
        fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch=2,
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL, verbose="TRUE"), error_message, fixed=TRUE)
})



context("pruningSample() results")


test_that("pruningSample() must return error when keepGDSpruned is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study", study.platform="WES",
                        stringsAsFactors=FALSE)

    error_message <- 'The \'keepGDSpruned\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gds=gdsFile, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned="YES",
        PATHSAMPLEGDS=NULL, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when keepFile is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                            study.desc="Pancreatic study", study.platform="WES",
                            stringsAsFactors=FALSE)

    error_message <- 'The \'keepFile\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gds=gdsFile, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=TRUE,
        PATHSAMPLEGDS=NULL, keepFile="NO", PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


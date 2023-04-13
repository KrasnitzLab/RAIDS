### Unit tests for synthetic functions

library(RAIDS)
library(withr)
library(testthat)



#############################################################################
### Tests select1KGPop() results
#############################################################################

context("select1KGPop() results")


test_that("select1KGPop() must return error when gdsReference is a character string", {

    error_message <- paste0("The \'gdsReference\' must be an object of class ",
                                "\'gds.class\'")

    expect_error(select1KGPop(gdsReference="test.gds", nbSamples=10L),
                    error_message)
})


test_that("select1KGPop() must return error when nbSamples is a character string", {

    dataDir <- system.file("extdata/tests", package="RAIDS")

    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'nbSamples\' parameter must be a ",
                                "single positive integer.")

    expect_error(select1KGPop(gdsReference=gdsF, nbSamples="CANADA"),
                    error_message)
})


test_that("select1KGPop() must return expected result", {

    dataDir <- system.file("extdata/tests", package="RAIDS")

    fileGDS <- file.path(dataDir, "1KG_Test_02.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    expected <- data.frame(sample.id=c("HG00101", "HG00097",
                                        "HG00096", "HG00100"),
                            pop.group=c("GBR", "GBR", "GBR", "GBR"),
                            superPop=rep("EUR", 4))

    set.seed(1212)
    results <- select1KGPop(gdsReference=gdsF, nbSamples=4L)

    expect_equal(results, expected)
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


test_that("splitSelectByPop() must return error when different numbers by population", {

    demo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC"),
                       pop.group=c("TSI", "TSI", "YRI"),
                       superPop=c("EUR", "EUR", "AFR"))

    error_message <- paste0("The number of samples in each subcontinental ",
                        "population has to be the same.\n")

    expect_error(splitSelectByPop(dataRef=demo), error_message)
})


test_that("splitSelectByPop() must return expected results", {

    demo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC", "SampleD",
                                    "SampleE", "SampleF"),
                       pop.group=c("TSI", "TSI", "YRI", "YRI", "LWK", "LWK"),
                       superPop=c("EUR", "EUR", "AFR", "AFR", "AFR", "AFR"))

    expected <- matrix(data=c("SampleA", "SampleB", "SampleC", "SampleD",
                              "SampleE", "SampleF"), byrow=FALSE, nrow=2)
    colnames(expected) <- c("TSI", "YRI", "LWK")

    result <- splitSelectByPop(dataRef=demo)

    expect_identical(result, expected)
})


#############################################################################
### Tests syntheticGeno() results
#############################################################################

context("syntheticGeno() results")


test_that("syntheticGeno() must return error when gds is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'gdsReference\' must be an object of class \'gds.class\'."

    expect_error(syntheticGeno(gdsReference=33, gdsRefAnnot=fileGDS,
            fileProfileGDS=2, data.id.profile="test",
            listSampleRef=c("Sample1", "sample2"),
            nbSim=1, prefId="", pRecomb=0.01,
            minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when gdsRefAnnot is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'gdsRefAnnot\' must be an object of class \'gds.class\'."

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=33,
                              fileProfileGDS=2, data.id.profile="test",
                              listSampleRef=c("Sample1", "sample2"),
                              nbSim=1, prefId="", pRecomb=0.01,
                              minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when fileProfileGDS is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'fileProfileGDS\' must be a character ",
                                "string and the file must exist.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                        fileProfileGDS=2, data.id.profile="test",
                        listSampleRef=c("Sample1", "sample2"),
                        nbSim=1, prefId="", pRecomb=0.01,
                        minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when data.id.profile is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'data.id.profile\' must be a character string."

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
            fileProfileGDS=fileGDS, data.id.profile=22,
            listSampleRef=c("Sample1", "sample2"),
            nbSim=1, prefId="", pRecomb=0.01,
            minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when listSampleRef is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'listSampleRef\' must be a vector of character strings."

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
        fileProfileGDS=fileGDS, data.id.profile="sample01", listSampleRef=33,
        nbSim=1, prefId="", pRecomb=0.01,
        minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when data.id.profile is a vector of strings", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'data.id.profile\' must be a character string."

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                    fileProfileGDS=fileGDS, data.id.profile=c("01", "02"),
                    listSampleRef=c("Sample1", "sample2"),
                    nbSim=1, prefId="", pRecomb=0.01,
                    minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when nbSim is a character string", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'nbSim\' parameter must be a single ",
                                "positive numeric value.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                        fileProfileGDS=fileGDS,
                        data.id.profile="test",
                        listSampleRef=c("Sample1", "sample2"),
                        nbSim="zero", prefId="", pRecomb=0.01,
                        minProb=0.999, seqError=0.001), error_message)
})


test_that(paste0("syntheticGeno() must return error when prefId is a vector ",
            "of character strings"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'prefId\' parameter must be a single ",
                                "character string.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                        fileProfileGDS=fileGDS,
                        data.id.profile="test",
                        listSampleRef=c("Sample1", "sample2"),
                        nbSim=2, prefId=c("", "eee"),
                        pRecomb=0.01, minProb=0.999, seqError=0.001),
                        error_message)
})


test_that("syntheticGeno() must return error when prefId is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'prefId\' parameter must be a single ",
                                "character string.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                        fileProfileGDS=fileGDS,
                        data.id.profile="test",
                        listSampleRef=c("Sample1", "sample2"),
                        nbSim=2, prefId=2,
                        pRecomb=0.01, minProb=0.999, seqError=0.001),
                 error_message)
})


test_that(paste0("syntheticGeno() must return error when pRecomb is a character string"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'pRecomb\' parameter must be a single positive ",
                                "numeric value between 0 and 1.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
            fileProfileGDS=fileGDS, data.id.profile="test",
            listSampleRef=c("Sample1", "sample2"),
            nbSim=1, prefId="", pRecomb="0.01",
            minProb=0.999, seqError=0.001), error_message)
})


test_that(paste0("syntheticGeno() must return error when seqError is a ",
                    "character string"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'seqError\' parameter must be a single ",
                "positive numeric value between 0 and 1.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                              fileProfileGDS=fileGDS,
                              data.id.profile="test",
                              listSampleRef=c("Sample1", "sample2"),
                              nbSim=1, prefId="", pRecomb=0.01,
                              minProb=0.999, seqError="0.001"), error_message)
})


test_that(paste0("syntheticGeno() must return error when minProb is a ",
                    "character string"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'minProb\' parameter must be a single ",
                            "positive numeric value between 0 and 1.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                              fileProfileGDS=fileGDS,
                              data.id.profile="test",
                              listSampleRef=c("Sample1", "sample2"),
                              nbSim=1, prefId="", pRecomb=0.01,
                              minProb="0.999", seqError=0.001), error_message)
})


#############################################################################
### Tests prepSynthetic() results
#############################################################################

context("prepSynthetic() results")


test_that(paste0("prepSynthetic() must return error when fileProfileGDS is",
                    " a numeric value"), {

    error_message <- paste0("The \'fileProfileGDS\' must be a character ",
                    "string representing the GDS Sample information file. ",
                    "The file must exist.")

    expect_error(prepSynthetic(fileProfileGDS=22, listSampleRef=c("S_1", "S_2"),
                    data.id.profile="S_1", studyDF="TCGA", nbSim=1L,
                    prefId=""), error_message)
})


test_that(paste0("prepSynthetic() must return error when nbSim is ",
                    "a character string"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    studyDF <- data.frame(study.id="Id of the study",
                            study.desc="Description",
                            study.platform="Whole-Exome",
                            stringsAsFactors=FALSE)

    error_message <- "The \'nbSim\' must be a single positive integer."

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
                    listSampleRef=c("S_1", "S_2"),
                    data.id.profile="S_1", studyDF=studyDF, nbSim="1L",
                    prefId=""), error_message)
})


test_that(paste0("prepSynthetic() must return error when listSampleRef is vector of numerics"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    error_message <- paste0("The \'listSampleRef\' must be a vector of ",
                                "character strings.")

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
                listSampleRef=c(1, 2),
                data.id.profile="S_1", studyDF=NULL, nbSim=1L,
                prefId=""), error_message)
})


test_that(paste0("prepSynthetic() must return error when studyDF is missing mandatory column"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    studyDF <- data.frame(study.id="Id of the study",
                                study.platform="Whole-Exome",
                                stringsAsFactors=FALSE)

    error_message <- paste0("The \'studyDF\' data frame is incomplete. ",
                                "One or more mandatory column is missing.\n")

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
                        listSampleRef=c("S_1", "S_2"),
                            data.id.profile="S_1", studyDF=studyDF, nbSim=1L,
                            prefId=""), error_message)
})


test_that(paste0("prepSynthetic() must return error when prefId is numeric"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    studyDF <- data.frame(study.id="Id of the study",
                          study.desc="Description",
                          study.platform="Whole-Exome",
                          stringsAsFactors=FALSE)

    error_message <- "The \'prefId\' must be a single character string."

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
                        listSampleRef=c("S_1", "S_2"),
                        data.id.profile="S_1", studyDF=studyDF, nbSim=1L,
                        prefId=33), error_message)
})

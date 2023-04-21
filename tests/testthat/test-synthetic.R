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

    expect_error(select1KGPop(gdsReference="test.gds", nbProfiles=10L),
                    error_message)
})


test_that("select1KGPop() must return error when nbProfiles is a character string", {

    dataDir <- system.file("extdata/tests", package="RAIDS")

    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'nbProfiles\' parameter must be a ",
                                "single positive integer.")

    expect_error(select1KGPop(gdsReference=gdsF, nbProfiles="CANADA"),
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
    results <- select1KGPop(gdsReference=gdsF, nbProfiles=4L)

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
            fileProfileGDS=2, profileID="test",
            listSampleRef=c("Sample1", "sample2"),
            nbSim=1, prefix="", pRecomb=0.01,
            minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when gdsRefAnnot is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'gdsRefAnnot\' must be an object of class \'gds.class\'."

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=33,
                              fileProfileGDS=2, profileID="test",
                              listSampleRef=c("Sample1", "sample2"),
                              nbSim=1, prefix="", pRecomb=0.01,
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
                        fileProfileGDS=2, profileID="test",
                        listSampleRef=c("Sample1", "sample2"),
                        nbSim=1, prefix="", pRecomb=0.01,
                        minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when profileID is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'profileID\' must be a character string."

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
            fileProfileGDS=fileGDS, profileID=22,
            listSampleRef=c("Sample1", "sample2"),
            nbSim=1, prefix="", pRecomb=0.01,
            minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when listSampleRef is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'listSampleRef\' must be a vector of character strings."

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
        fileProfileGDS=fileGDS, profileID="sample01", listSampleRef=33,
        nbSim=1, prefix="", pRecomb=0.01,
        minProb=0.999, seqError=0.001), error_message)
})


test_that("syntheticGeno() must return error when profileID is a vector of strings", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'profileID\' must be a character string."

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                    fileProfileGDS=fileGDS, profileID=c("01", "02"),
                    listSampleRef=c("Sample1", "sample2"),
                    nbSim=1, prefix="", pRecomb=0.01,
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
                        fileProfileGDS=fileGDS, profileID="test",
                        listSampleRef=c("Sample1", "sample2"),
                        nbSim="zero", prefix="", pRecomb=0.01,
                        minProb=0.999, seqError=0.001), error_message)
})


test_that(paste0("syntheticGeno() must return error when prefId is a vector ",
            "of character strings"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'prefix\' parameter must be a single ",
                                "character string.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                        fileProfileGDS=fileGDS, profileID="test",
                        listSampleRef=c("Sample1", "sample2"),
                        nbSim=2, prefix=c("", "eee"),
                        pRecomb=0.01, minProb=0.999, seqError=0.001),
                        error_message)
})


test_that("syntheticGeno() must return error when prefix is a numeric value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'prefix\' parameter must be a single ",
                                "character string.")

    expect_error(syntheticGeno(gdsReference=gdsF, gdsRefAnnot=gdsF,
                        fileProfileGDS=fileGDS, profileID="test",
                        listSampleRef=c("Sample1", "sample2"),
                        nbSim=2, prefix=2,
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
            fileProfileGDS=fileGDS, profileID="test",
            listSampleRef=c("Sample1", "sample2"),
            nbSim=1, prefix="", pRecomb="0.01",
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
                              profileID="test",
                              listSampleRef=c("Sample1", "sample2"),
                              nbSim=1, prefix="", pRecomb=0.01,
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
        fileProfileGDS=fileGDS, profileID="test",
        listSampleRef=c("Sample1", "sample2"), nbSim=1, prefix="",
        pRecomb=0.01, minProb="0.999", seqError=0.001), error_message)
})


test_that(paste0("syntheticGeno() must return error expected results"), {

    set.seed(121)

    dataDirSample <- test_path("fixtures/sampleGDSforEstimAlleFraction")
    file.copy(file.path(dataDirSample, "ex1_demoForEstimAllFrac.gds"),
              file.path(dataDirSample, "ex1.gds"))
    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                 envir=parent.frame())

    dataDirRef <- test_path("fixtures")

    ## Open  1KG Annotation GDS file
    fileRefAnnot <- file.path(dataDirRef, "ex1_good_small_1KG_Annot_GDS.gds")
    gdsRefAnnot <- openfn.gds(fileRefAnnot)
    withr::defer((closefn.gds(gdsRefAnnot)), envir=parent.frame())

    ## Open 1KG GDS file
    gds1KG <- snpgdsOpen(file.path(dataDirRef, "ex1_good_small_1KG_GDS.gds"))
    withr::defer((closefn.gds(gds1KG)), envir=parent.frame())

    synthStudyDF <- data.frame(study.id="MYDATA.Synthetic",
                                study.desc="MYDATA synthetic data",
                                study.platform="PLATFORM",
                                stringsAsFactors=FALSE)

    result1 <- prepSynthetic(fileProfileGDS=file.path(dataDirSample, "ex1.gds"),
        listSampleRef=c("HG00243", "HG00149"), profileID="ex1",
        studyDF=synthStudyDF, nbSim=1L, prefix="test2", verbose=FALSE)

    result2 <- syntheticGeno(gdsReference=gds1KG, gdsRefAnnot=gdsRefAnnot,
        fileProfileGDS=file.path(dataDirSample, "ex1.gds"), profileID="ex1",
        listSampleRef=c("HG00243", "HG00149"), nbSim=1, prefix="test2",
        pRecomb=0.01, minProb=0.999, seqError=0.001)

    expect_equal(result2, 0L)

    profileGDS <- openfn.gds(file.path(dataDirSample, "ex1.gds"))
    withr::defer((closefn.gds(profileGDS)), envir=parent.frame())

    sampleList <- read.gdsn(index.gdsn(profileGDS, "sample.id"))

    expect_equal(sampleList[158:159],
                 c("test2.ex1.HG00243.1", "test2.ex1.HG00149.1"))

    genotype <- read.gdsn(index.gdsn(profileGDS, "genotype"))

    expect_equal(genotype[,159], c(rep(0, 9), 2, rep(0, 5), 1, 2, rep(0, 13),
        1, 0, 0, 2, rep(0, 6), 3,  rep(0, 6), 1, rep(0, 3), 3, rep(0, 9)))

    expect_equal(genotype[,158], c(3, rep(0, 15), 2, 3, 0, 1,  rep(0, 13), 3,
        rep(0, 27)))
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
        profileID="S_1", studyDF="TCGA", nbSim=1L, prefix=""), error_message)
})


test_that(paste0("prepSynthetic() must return error when profileID is numeric"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    studyDF <- data.frame(study.id="Id of the study",
                study.desc="Description", study.platform="Whole-Exome",
                stringsAsFactors=FALSE)

    error_message <- "The \'profileID\' must be a single character string."

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
            listSampleRef=c("S_1", "S_2"), profileID=33, studyDF=studyDF,
            nbSim=1L, prefix="test"), error_message)
})


test_that(paste0("prepSynthetic() must return error when nbSim is ",
                    "a character string"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    studyDF <- data.frame(study.id="Id of the study",
            study.desc="Description", study.platform="Whole-Exome",
            stringsAsFactors=FALSE)

    error_message <- "The \'nbSim\' must be a single positive integer."

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
        listSampleRef=c("S_1", "S_2"), profileID="S_1",
        studyDF=studyDF, nbSim="1L", prefix=""), error_message)
})


test_that(paste0("prepSynthetic() must return error when listSampleRef is vector of numerics"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    error_message <- paste0("The \'listSampleRef\' must be a vector of ",
                                "character strings.")

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
        listSampleRef=c(1, 2),  profileID="S_1", studyDF=NULL, nbSim=1L,
        prefix=""), error_message)
})


test_that(paste0("prepSynthetic() must return error when studyDF is missing mandatory column"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    studyDF <- data.frame(study.id="Id of the study",
        study.platform="Whole-Exome", stringsAsFactors=FALSE)

    error_message <- paste0("The \'studyDF\' data frame is incomplete. ",
                                "One or more mandatory column is missing.\n")

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
        listSampleRef=c("S_1", "S_2"), profileID="S_1", studyDF=studyDF,
        nbSim=1L, prefix=""), error_message)
})


test_that(paste0("prepSynthetic() must return error when prefix is numeric"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    studyDF <- data.frame(study.id="Id of the study",
        study.desc="Description", study.platform="Whole-Exome",
        stringsAsFactors=FALSE)

    error_message <- "The \'prefix\' must be a single character string."

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
        listSampleRef=c("S_1", "S_2"), profileID="S_1", studyDF=studyDF,
        nbSim=1L, prefix=33), error_message)
})


test_that(paste0("prepSynthetic() must return error when verbose is numeric"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    studyDF <- data.frame(study.id="Id of the study",
        study.desc="Description", study.platform="Whole-Exome",
        stringsAsFactors=FALSE)

    error_message <- "The \'verbose\' parameter must be a logical (TRUE or FALSE)."

    expect_error(prepSynthetic(fileProfileGDS=fileGDS,
        listSampleRef=c("S_1", "S_2"), profileID="S_1", studyDF=studyDF,
        nbSim=1L, prefix="test", verbose=33), error_message, fixed=TRUE)
})


test_that(paste0("prepSynthetic() must return expected results"), {

    dataDirSample <- test_path("fixtures/sampleGDSforEstimAlleFraction")
    file.copy(file.path(dataDirSample, "ex1_demoForEstimAllFrac.gds"),
              file.path(dataDirSample, "ex1.gds"))
    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                 envir=parent.frame())

    synthStudyDF <- data.frame(study.id="MYDATA.Synthetic",
                                    study.desc="MYDATA synthetic data",
                                    study.platform="PLATFORM",
                                    stringsAsFactors=FALSE)

    result <- prepSynthetic(fileProfileGDS=file.path(dataDirSample, "ex1.gds"),
        listSampleRef=c("HG00243", "HG00149"), profileID="ex1",
        studyDF=synthStudyDF,
        nbSim=1L, prefix="test", verbose=FALSE)

    profileGDS <- openfn.gds(file.path(dataDirSample, "ex1.gds"))
    withr::defer((closefn.gds(profileGDS)), envir=parent.frame())

    studyAnnot <- read.gdsn(index.gdsn(profileGDS, "study.annot"))

    studyList <- read.gdsn(index.gdsn(profileGDS, "study.list"))

    expect_equal(studyAnnot$data.id[158:159],
                   c("test.ex1.HG00243.1", "test.ex1.HG00149.1"))
    expect_equal(studyAnnot$case.id[158:159],  c("HG00243", "HG00149"))
    expect_equal(studyAnnot$sample.type[158:159],  c("Synthetic", "Synthetic"))
    expect_equal(studyAnnot$diagnosis[158:159],  c("Cancer", "Cancer"))
    expect_equal(studyAnnot$source[158:159],  c("Synthetic", "Synthetic"))
    expect_equal(studyAnnot$study.id[158:159],  rep("MYDATA.Synthetic", 2))

    expect_equal(studyList$study.id[3],  "MYDATA.Synthetic")
    expect_equal(studyList$study.desc[3],  "MYDATA synthetic data")
    expect_equal(studyList$study.platform[3],  "Synthetic")

    expect_equal(result, 0L)
})


test_that(paste0("prepSynthetic() must return error when profileID not present"), {

    dataDirSample <- test_path("fixtures/sampleGDSforEstimAlleFraction")
    file.copy(file.path(dataDirSample, "ex1_demoForEstimAllFrac.gds"),
              file.path(dataDirSample, "ex1.gds"))
    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                 envir=parent.frame())

    synthStudyDF <- data.frame(study.id="MYDATA.Synthetic",
                                study.desc="MYDATA synthetic data",
                                study.platform="PLATFORM",
                                stringsAsFactors=FALSE)

    error_message <- paste0("The profile \'toto\' is not found in the annoted ",
            "study present in the Profile GDS file.\n")

    expect_error(prepSynthetic(fileProfileGDS=file.path(dataDirSample, "ex1.gds"),
        listSampleRef=c("HG00243", "HG00149"), profileID="toto",
        studyDF=synthStudyDF, nbSim=1L, prefix="test", verbose=FALSE),
        error_message, fixed=TRUE)
})


test_that(paste0("prepSynthetic() must return error when prefix already used"), {

    dataDirSample <- test_path("fixtures/sampleGDSforEstimAlleFraction")
    file.copy(file.path(dataDirSample, "ex1_demoForEstimAllFrac.gds"),
              file.path(dataDirSample, "ex1.gds"))
    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                 envir=parent.frame())

    synthStudyDF <- data.frame(study.id="MYDATA.Synthetic",
                                study.desc="MYDATA synthetic data",
                                study.platform="PLATFORM",
                                stringsAsFactors=FALSE)

    result <- prepSynthetic(fileProfileGDS=file.path(dataDirSample, "ex1.gds"),
        listSampleRef=c("HG00243", "HG00149"), profileID="ex1",
        studyDF=synthStudyDF, nbSim=1L, prefix="test", verbose=FALSE)

    expect_equal(result, 0L)

    error_message <- "Error data.id of the simulation exists change prefix\n"

    expect_error(prepSynthetic(fileProfileGDS=file.path(dataDirSample, "ex1.gds"),
        listSampleRef=c("HG00243", "HG00149"), profileID="ex1",
        studyDF=synthStudyDF, nbSim=1L, prefix="test", verbose=FALSE),
        error_message, fixed=TRUE)
})

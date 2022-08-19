### Unit tests for processStudy.R functions

library(RAIDS)
library(withr)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L, 181538259L,
                170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
                135086622L, 133275309L, 114364328L, 107043718L, 101991189L,
                90338345L, 83257441L, 80373285L, 58617616L, 64444167L,
                46709983L,  50818468L, 156040895L, 57227415L,  16569L)

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


#############################################################################
### Tests pruningSample() results
#############################################################################


context("pruningSample() results")


test_that("pruningSample() must return error when gds is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    error_message <- "The \'gds\' must be an object of class \'gds.class\'."

    expect_error(pruningSample(gds=gdsFile, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when keepGDSpruned is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- 'The \'keepGDSpruned\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned="YES",
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when keepFile is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- 'The \'keepFile\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="1KG_Test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=TRUE,
        PATHSAMPLEGDS=data.dir, keepFile="NO", PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when np is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'np\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np="1", verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when slide.max.bp.v is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'slide.max.bp.v\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v="4", ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when ld.threshold.v is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'ld.threshold.v\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=4, ld.threshold.v="3",
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when ld.threshold.v is a vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'ld.threshold.v\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=4, ld.threshold.v=c(3,3),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when method is a numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'method\' parameter must be a character string."

    expect_error(pruningSample(gds=gdsF, method=3, sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when method is not in the list of choices", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    expect_error(pruningSample(gds=gdsF, method="test", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"))
})


test_that("pruningSample() must return error when sampleCurrent is a numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'sampleCurrent\' parameter must be a character string."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent=2,
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})



test_that("pruningSample() must return error when PATHPRUNED is a numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'PATHPRUNED\' parameter must be ",
        "a character string representing an existing directory.")

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="Sample2",
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
        np=1L, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=2, outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when PATHSAMPLEGDS is a numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'PATHSAMPLEGDS\' parameter must be a character string ",
                                "representing an existing directory.")

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="Sample2",
            study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
            np=1L, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
            PATHSAMPLEGDS=33, keepFile=FALSE, PATHPRUNED=data.dir, outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when PATHSAMPLEGDS is a non existing directory", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'PATHSAMPLEGDS\' parameter must be a character string ",
                                "representing an existing directory.")

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="Sample2",
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
        np=1L, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=paste0(data.dir, "_NOT_EXISTING_DIRECTORY"),
        keepFile=FALSE, PATHPRUNED=data.dir, outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when GDS Sample file does not exist", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFIle <- file.path(data.dir, "1KG_Test.gds")
    sampleRDS <- file.path(data.dir, "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    gdsSampleFile <- file.path(data.dir, "A_File_That_DOES_NOT_EXIST.gds")

    error_message <- paste0("The GDS Sample file \'", gdsSampleFile, " does not exist.")

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="A_File_That_DOES_NOT_EXIST",
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L,
        ld.threshold.v=sqrt(0.1), np=1L, verbose.v=FALSE, chr=NULL,
        minAF.SuperPop=NULL, keepGDSpruned=FALSE, PATHSAMPLEGDS=data.dir,
        keepFile=FALSE, PATHPRUNED=data.dir, outPref="pruned"), error_message,
        fixed=TRUE)
})



#############################################################################
### Tests add1KG2SampleGDS() results
#############################################################################


context("add1KG2SampleGDS() results")


test_that("add1KG2SampleGDS() must return error when gds is a character string", {

    error_message <- "The \'gds\' must be an object of class \'gds.class\'."

    expect_error(add1KG2SampleGDS(gds="toto.gds", gdsSampleFile="sample.gds",
        sampleCurrent="sample", study.id="TCGA"), error_message, fixed=TRUE)
})


test_that("add1KG2SampleGDS() must return error when gdsSampleFile is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'gdsSampleFile\' must be a character string representing ",
                                "the GDS Sample file. The file must exist.")

    expect_error(add1KG2SampleGDS(gds=gdsF, gdsSampleFile=33,
        sampleCurrent="sample", study.id="TCGA"), error_message, fixed=TRUE)
})

test_that("add1KG2SampleGDS() must return error when sampleCurrent is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'sampleCurrent\' must be a character string."

    expect_error(add1KG2SampleGDS(gds=gdsF, gdsSampleFile=gdsFile,
            sampleCurrent=33, study.id="TCGA"), error_message, fixed=TRUE)
})


test_that("add1KG2SampleGDS() must return error when study.id is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'study.id\' must be a character string."

    expect_error(add1KG2SampleGDS(gds=gdsF, gdsSampleFile=gdsFile,
        sampleCurrent="Test", study.id=22), error_message, fixed=TRUE)
})




#############################################################################
### Tests computePCARefSample() results
#############################################################################


context("computePCARefSample() results")


test_that("computePCARefSample() must return error when np is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
                        study.id.ref="Ref.1KG",
                        np="1", algorithm="exact",
                        eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when np is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
                                study.id.ref="Ref.1KG",
                                np=1L, algorithm=33,
                                eigen.cnt=32L), error_message, fixed=TRUE)
})



test_that("computePCARefSample() must return error when algorithm is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
                                     study.id.ref="Ref.1KG",
                                     np=1L, algorithm=33,
                                     eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when algorithm is not a valid choice", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
                                        study.id.ref="Ref.1KG",
                                        np=1L, algorithm="sun",
                                        eigen.cnt=32L))
})


test_that("computePCARefSample() must return error when name.id isnumeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'name.id\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id=22,
                                    study.id.ref="Ref.1KG",
                                    np=1L, algorithm="exact",
                                    eigen.cnt=32L), error_message, fixed=TRUE)
})



#############################################################################
### Tests addStudy1Kg() results
#############################################################################


context("addStudy1Kg() results")


test_that("addStudy1Kg() must return error when gds is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    error_message <- "The \'gds\' must be an object of class \'gds.class\'."

    expect_error(addStudy1Kg(gds=33, gdsSampleFile=gdsFile),
                                error_message, fixed=TRUE)
})


#############################################################################
### Tests estimateAllelicFraction() results
#############################################################################


context("estimateAllelicFraction() results")


test_that("estimateAllelicFraction() must return error when gds is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'gds\' must be an object of class \'gds.class\'."

    expect_error(estimateAllelicFraction(gds="test.gds",
                    gdsSample=gdsF, sampleCurrent="gds", study.id= "TCGA",
                    chrInfo=chrInfo, studyType="DNA", minCov=10, minProb=0.999,
                    eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                    wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when gdsSample is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'gdsSample\' must be an object of class \'gds.class\'."

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample="test.gds", sampleCurrent="gds", study.id= "TCGA",
                chrInfo=chrInfo, studyType="DNA", minCov=10, minProb=0.999,
                eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when sampleCurrent is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'sampleCurrent\' must be a character string."

    expect_error(estimateAllelicFraction(gds=gdsF,
                    gdsSample=gdsF, sampleCurrent=33, study.id= "TCGA",
                    chrInfo=chrInfo, studyType="DNA", minCov=10, minProb=0.999,
                    eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                    wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when study.id is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'study.id\' must be a character string."

    expect_error(estimateAllelicFraction(gds=gdsF,
                    gdsSample=gdsF, sampleCurrent="test", study.id= 33,
                    chrInfo=chrInfo, studyType="DNA", minCov=10, minProb=0.999,
                    eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                    wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when studyType is a numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'studyType\' must be a character string."

    expect_error(estimateAllelicFraction(gds=gdsF,
            gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
            chrInfo=chrInfo, studyType=33, minCov=10, minProb=0.999,
            eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
            wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minCov is vector of numeric values", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'minCov\' must be a single numeric positive value"

    expect_error(estimateAllelicFraction(gds=gdsF,
            gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
            chrInfo=chrInfo, studyType="DNA", minCov=c(10, 11), minProb=0.999,
            eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
            wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minCov is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'minCov\' must be a single numeric positive value"

    expect_error(estimateAllelicFraction(gds=gdsF,
            gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
            chrInfo=chrInfo, studyType="DNA", minCov="3", minProb=0.999,
            eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
            wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minProb is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                      "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb="0.1",
                eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minProb is negative value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gds=gdsF,
            gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
            chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=-0.1,
            eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
            wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minProb is above 1", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=1.01,
                eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minProb is vector of numeric values", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gds=gdsF,
            gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
            chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=c(0.01, 0.02),
            eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
            wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
                eProb="0.001", cutOffLOH=-5, cutOffHomoScore=-3,
                wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is negative value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                    "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gds=gdsF,
            gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
            chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
            eProb=-0.01, cutOffLOH=-5, cutOffHomoScore=-3,
            wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is above 1", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
                eProb=1.01, cutOffLOH=-5, cutOffHomoScore=-3,
                wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is vector of numeric values", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gds=gdsF,
            gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
            chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
            eProb=c(0.01, 0.001), cutOffLOH=-5, cutOffHomoScore=-3,
            wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is vector of numeric values", {

    data.dir <- system.file("extdata/tests", package="RAIDS")

    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'wAR\' must be a single numeric positive value."

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
                eProb=0.01, cutOffLOH=-5, cutOffHomoScore=-3,
                wAR="9"), error_message, fixed=TRUE)
})


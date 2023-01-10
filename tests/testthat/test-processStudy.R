### Unit tests for processStudy.R functions

library(RAIDS)
library(withr)
library(gdsfmt)

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

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gds=gdsFIle, listPCA=list(),
                    sample.current="sample1",
                    np="test"), error_message)
})

test_that("projectSample2PCA() must return error when np is negative integer", {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gds=gdsFIle, listPCA=list(),
                                    sample.current="sample1",
                                    np=-1L), error_message)
})

test_that("projectSample2PCA() must return error when np is zero", {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gds=gdsFIle, listPCA=list(),
                                    sample.current="sample1",
                                    np=0L), error_message)
})

test_that("projectSample2PCA() must return error when sample.current is number", {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")

    error_message <- paste0("The \'sample.current\' ",
                                "parameter must be a character string.")

    expect_error(projectSample2PCA(gds=gdsFIle, listPCA=list(),
                                        sample.current=101,
                                        np=1L), error_message)
})


#############################################################################
### Tests appendStudy2GDS1KG() results
#############################################################################


context("appendStudy2GDS1KG() results")

test_that("appendStudy2GDS1KG() must return error when PATHGENO is a numeric", {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                            study.desc="Pancreatic study", study.platform="WES",
                            stringsAsFactors=FALSE)

    error_message <- paste0("The \'PATHGENO\' must be a character string ",
                            "representing a path. The path must exist.")

    expect_error(appendStudy2GDS1KG(PATHGENO=22, fileNamePED=sampleRDS,
                fileNameGDS=gdsFIle, batch=2,
                studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL,
                verbose="TRUE"), error_message, fixed=TRUE)
})


test_that("appendStudy2GDS1KG() must return error when fileNamePED is numeric", {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study", study.platform="WES",
                        stringsAsFactors=FALSE)

    error_message <- paste0("The \'fileNamePED\' must be a character string ",
        "representing the RDS Sample information file. The file must exist.")

    expect_error(appendStudy2GDS1KG(PATHGENO=test_path("fixtures"),
        fileNamePED=21, fileNameGDS=gdsFIle, batch=1,
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL,
        verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when fileNameGDS is numeric", {

    gdsFIle <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- paste0("The \'fileNameGDS\' must be a character string ",
                        "representing the GDS 1KG file. The file must exist.")

    expect_error(appendStudy2GDS1KG(PATHGENO=test_path("fixtures"),
        fileNamePED=sampleRDS, fileNameGDS=33, batch=1,
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL,
        verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when batch is a vector of numerics", {

    gdsFIle <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES", study.desc="Pancreatic",
                                study.platform="WES", stringsAsFactors=FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(appendStudy2GDS1KG(PATHGENO=test_path("fixtures"),
        fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch=c(1,2),
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL,
        verbose=TRUE), error_message)
})

test_that("appendStudy2GDS1KG() must return error when batch is a character string", {

    gdsFIle <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures",  "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(appendStudy2GDS1KG(PATHGENO=test_path("fixtures"),
        fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch="2",
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL,
        verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when studyDF is missing mandatory column", {

    gdsFIle <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures",  "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.descption="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- paste0("The \'studyDF\' must be a data.frame and contain ",
        "those 3 columns: \'study.id\', \'study.desc\' and \'study.platform\'.")

    expect_error(appendStudy2GDS1KG(PATHGENO=test_path("fixtures"),
            fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch=1,
            studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL,
            verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when listSamples is a numeric", {

    gdsFIle <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures",  "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study", study.platform="WES",
                        stringsAsFactors=FALSE)

    error_message <- paste0("The \'listSamples\' must be a vector ",
        "of character strings (1 entry or more) or NULL.")

    expect_error(appendStudy2GDS1KG(PATHGENO=test_path("fixtures"),
        fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch=2,
        studyDF=studyInfo, listSamples=33, PATHSAMPLEGDS=NULL,
        verbose=FALSE), error_message, fixed=TRUE)
})


test_that("appendStudy2GDS1KG() must return error when verbose is a character string", {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- 'The \'verbose\' parameter must be a logical (TRUE or FALSE).'

    expect_error(appendStudy2GDS1KG(PATHGENO=test_path("fixtures"),
        fileNamePED=sampleRDS, fileNameGDS=gdsFIle, batch=2,
        studyDF=studyInfo, listSamples=NULL, PATHSAMPLEGDS=NULL,
        verbose="TRUE"), error_message, fixed=TRUE)
})


#############################################################################
### Tests pruningSample() results
#############################################################################


context("pruningSample() results")


test_that("pruningSample() must return error when gds is a character string", {

    data.dir <- test_path("fixtures")

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    error_message <- "The \'gds\' must be an object of class \'gds.class\'."

    expect_error(pruningSample(gds=gdsFile, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when keepGDSpruned is a character string", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures",  "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- 'The \'keepGDSpruned\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned="YES",
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when keepFile is a character string", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- 'The \'keepFile\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="1KG_Test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=TRUE,
        PATHSAMPLEGDS=data.dir, keepFile="NO", PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when np is a character string", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'np\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=5e5, ld.threshold.v=sqrt(0.1),
        np="1", verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when slide.max.bp.v is a character string", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'slide.max.bp.v\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v="4", ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when ld.threshold.v is a character string", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'ld.threshold.v\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=4, ld.threshold.v="3",
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when ld.threshold.v is a vector of numerics", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'ld.threshold.v\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=4, ld.threshold.v=c(3,3),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when method is a numeric", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'method\' parameter must be a character string."

    expect_error(pruningSample(gds=gdsF, method=3, sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when method is not in the list of choices", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    expect_error(pruningSample(gds=gdsF, method="test", sampleCurrent="test",
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"))
})


test_that("pruningSample() must return error when sampleCurrent is a numeric", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'sampleCurrent\' parameter must be a character string."

    expect_error(pruningSample(gds=gdsF, method="corr", sampleCurrent=2,
        study.id="test", listSNP=NULL, slide.max.bp.v=50000L, ld.threshold.v=sqrt(0.1),
        np=1, verbose.v=FALSE, chr=NULL, minAF.SuperPop=NULL, keepGDSpruned=FALSE,
        PATHSAMPLEGDS=data.dir, keepFile=FALSE, PATHPRUNED=".", outPref="pruned"), error_message, fixed=TRUE)
})



test_that("pruningSample() must return error when PATHPRUNED is a numeric", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

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

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

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

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

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

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'gdsSampleFile\' must be a character string representing ",
                                "the GDS Sample file. The file must exist.")

    expect_error(add1KG2SampleGDS(gds=gdsF, gdsSampleFile=33,
        sampleCurrent="sample", study.id="TCGA"), error_message, fixed=TRUE)
})

test_that("add1KG2SampleGDS() must return error when sampleCurrent is a numeric value", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'sampleCurrent\' must be a character string."

    expect_error(add1KG2SampleGDS(gds=gdsF, gdsSampleFile=gdsFile,
            sampleCurrent=33, study.id="TCGA"), error_message, fixed=TRUE)
})


test_that("add1KG2SampleGDS() must return error when study.id is a numeric value", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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


test_that("computePCARefSample() must return error when gdsSample isnumeric value", {

    error_message <- "The \'gdsSample\' must be an object of class \'gds.class\'."

    expect_error(computePCARefSample(gdsSample="test.gds", name.id="test",
                            study.id.ref="Ref.1KG", np=1L, algorithm="exact",
                            eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when name.id is numeric value", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'name.id\' parameter must be a single character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id=22,
                            study.id.ref="Ref.1KG", np=1L, algorithm="exact",
                            eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when study.id.ref is numeric value", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'study.id.ref\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="Synthetic",
                            study.id.ref=33, np=1L, algorithm="exact",
                            eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when np is a character string", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
                        study.id.ref="Ref.1KG",
                        np="1", algorithm="exact",
                        eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when np is a numeric value", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
                                study.id.ref="Ref.1KG", np=1L, algorithm=33,
                                eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when algorithm is a numeric value", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
            study.id.ref="Ref.1KG", np=1L, algorithm=33,
            eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when algorithm is not a valid choice", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
            study.id.ref="Ref.1KG", np=1L, algorithm="sun", eigen.cnt=32L))
})


test_that("computePCARefSample() must return error when eigen.cnt is a string", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'eigen.cnt\' parameter must be a single integer."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
        study.id.ref="Ref.1KG", np=1L, algorithm="sun", eigen.cnt="32L"),
        error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when missing.rate is negative value", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'missing.rate\' must be a single numeric ",
                                "positive value between 0 and 1 or NaN.")

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
        study.id.ref="Ref.1KG", np=1L, algorithm="sun", eigen.cnt=32L,
        missing.rate=-0.02), error_message, fixed=TRUE)
})



#############################################################################
### Tests addStudy1Kg() results
#############################################################################


context("addStudy1Kg() results")


test_that("addStudy1Kg() must return error when gds is a numeric value", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'gds\' must be an object of class \'gds.class\'."

    expect_error(addStudy1Kg(gds=33, gdsSampleFile=gdsFile),
                                error_message, fixed=TRUE)
})


test_that("addStudy1Kg() must return error when gdsSampleFile is a numeric value", {

    ## Create a temporary GDS file in an test directory
    gdsFile <- test_path("fixtures", "GDS_TEMP_processStudy_101.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_Sample_file(gdsFile)

    error_message <- paste0("The \'gdsSampleFile\' must be a character ",
        "string representing the GDS Sample file. The file must exist.")

    expect_error(addStudy1Kg(gds=GDS_file_tmp, gdsSampleFile=33),
                    error_message, fixed=TRUE)

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)
})


test_that("addStudy1Kg() must return expected results", {

    ## Create a temporary GDS file in an test directory
    gdsFile1KG <- test_path("fixtures", "GDS_TEMP_processStudy_1KG_102.gds")

    ## Create and open a temporary GDS file 1KG
    GDS_file_tmp_1KG  <- local_GDS_1KG_file(gdsFile1KG, env=parent.frame())

    ## Create and open a temporary GDS Sample file
    gdsFileSample <- test_path("fixtures",
                                "GDS_TEMP_processStudy_Sample_102.gds")
    GDS_file_Sample <- createfn.gds(gdsFileSample)

    study.list <- data.frame(study.id=c("HTT Study"),
                        study.desc=c("Important Study"),
                        study.platform=c("Panel"), stringsAsFactors = FALSE)

    add.gdsn(GDS_file_Sample, "study.list", study.list)

    study.annot <- data.frame(data.id=c("TOTO1"), case.id=c("TOTO1"),
                sample.type=c("Study"), diagnosis=c("Study"),
                source=rep("IGSR"), study.id=c("Study"),
                stringsAsFactors=FALSE)

    add.gdsn(GDS_file_Sample, "study.annot", study.annot)

    sync.gds(GDS_file_Sample)

    closefn.gds(GDS_file_Sample)
    withr::defer((unlink(gdsFileSample, force = TRUE)), envir = parent.frame())

    result0 <- addStudy1Kg(gds=GDS_file_tmp_1KG, gdsSampleFile=gdsFileSample)

    gds_sample_file <- openfn.gds(gdsFileSample, readonly = TRUE)

    result1 <- read.gdsn(index.gdsn(node=gds_sample_file, path="study.list"))

    result2 <- read.gdsn(index.gdsn(node=gds_sample_file, path="study.annot"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp_1KG)
    closefn.gds(gdsfile=gds_sample_file)

    expected1 <- data.frame(study.id=c("HTT Study", "Ref.1KG"),
        study.desc=c("Important Study", "Unrelated samples from 1000 Genomes"),
        study.platform=c("Panel", "GRCh38 1000 genotypes"),
        stringsAsFactors=FALSE)

    expected2 <- data.frame(data.id=c("TOTO1", "HTT101", "HTT103"),
                case.id=c("TOTO1", "HTT101", "HTT103"),
                sample.type=c("Study", rep("Reference", 2)),
                diagnosis=c("Study", rep("Reference", 2)),
                source=rep("IGSR", 3), study.id=c("Study", "Ref.1KG", "Ref.1KG"),
                stringsAsFactors=FALSE)

    expect_equal(result0, 0L)
    expect_equal(result1, expected1)
    expect_equal(result2, expected2)
})



test_that("addStudy1Kg() must return expected results when 1KG already present", {

    ## Create a temporary GDS file in an test directory
    gdsFile1KG <- test_path("fixtures", "GDS_TEMP_processStudy_1KG_02.gds")

    ## Create and open a temporary GDS file 1KG
    GDS_file_tmp_1KG  <- local_GDS_1KG_file(gdsFile1KG, env=parent.frame())

    ## Create and open a temporary GDS Sample file
    gdsFileSample <- test_path("fixtures",
                                "GDS_TEMP_processStudy_Sample_02.gds")
    GDS_file_Sample <- createfn.gds(gdsFileSample)

    study.list <- data.frame(study.id=c("Ref.1KG"),
                        study.desc=c("Important Study"),
                        study.platform=c("Panel"), stringsAsFactors = FALSE)

    add.gdsn(GDS_file_Sample, "study.list", study.list)

    study.annot <- data.frame(data.id=c("TOTO1"), case.id=c("TOTO1"),
                             sample.type=c("Study"), diagnosis=c("Study"),
                             source=rep("IGSR"), study.id=c("Ref.1KG"),
                             stringsAsFactors=FALSE)

    add.gdsn(GDS_file_Sample, "study.annot", study.annot)

    sync.gds(GDS_file_Sample)

    closefn.gds(GDS_file_Sample)
    withr::defer((unlink(gdsFileSample)), envir = parent.frame())

    result0 <- addStudy1Kg(gds=GDS_file_tmp_1KG, gdsSampleFile=gdsFileSample)

    gds_sample_file <- openfn.gds(gdsFileSample, readonly = TRUE)

    result1 <- read.gdsn(index.gdsn(node=gds_sample_file, path="study.list"))

    result2 <- read.gdsn(index.gdsn(node=gds_sample_file, path="study.annot"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp_1KG)
    closefn.gds(gdsfile=gds_sample_file)

    expected1 <- study.list

    expected2 <- study.annot

    expect_equal(result0, 0L)
    expect_equal(result1, expected1)
    expect_equal(result2, expected2)
})


#############################################################################
### Tests estimateAllelicFraction() results
#############################################################################


context("estimateAllelicFraction() results")


test_that("estimateAllelicFraction() must return error when gds is a character string", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

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


test_that("estimateAllelicFraction() must return error when wAR is vector of numeric values", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'wAR\' parameter must be a single numeric positive value."

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
                eProb=0.01, cutOffLOH=-5, cutOffHomoScore=-3,
                wAR="9"), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when cutOffLOH is vector of numeric values", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'cutOffLOH\' parameter must be a single numeric value."

    expect_error(estimateAllelicFraction(gds=gdsF,
            gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
            chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
            eProb=0.01, cutOffLOH=c(-5, -3), cutOffHomoScore=-3,
            wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when cutOffLOH is a character string", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'cutOffLOH\' parameter must be a single numeric value."

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
                eProb=0.01, cutOffLOH="-3", cutOffHomoScore=-3,
                wAR=9), error_message, fixed=TRUE)
})




test_that("estimateAllelicFraction() must return error when cutOffLOH is vector of numeric values", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'cutOffHomoScore\' parameter must be a single numeric value."

    expect_error(estimateAllelicFraction(gds=gdsF,
                gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
                eProb=0.01, cutOffLOH=-5, cutOffHomoScore=c(-3, -4),
                wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when cutOffHomoScore is a character string", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'cutOffHomoScore\' parameter must be a single numeric value."

    expect_error(estimateAllelicFraction(gds=gdsF,
                    gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
                    chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
                    eProb=0.01, cutOffLOH=-5, cutOffHomoScore="-3",
                    wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when studyType is not a valid choice", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    expect_error(estimateAllelicFraction(gds=gdsF,
        gdsSample=gdsF, sampleCurrent="test", study.id="DNA_1",
        chrInfo=chrInfo, studyType="TOTO", minCov=10L, minProb=0.1,
        eProb=0.01, cutOffLOH=-5, cutOffHomoScore=-3, wAR=9))
})

#############################################################################
### Tests createStudy2GDS1KG() results
#############################################################################

context("createStudy2GDS1KG() results")


test_that(paste0("createStudy2GDS1KG() must return error when fileNamePED is",
            " a numeric value and pedStudy is NULL"), {

    error_message <- paste0("The \'fileNamePED\' must be a character string ",
            "representing the RDS Sample information file. ",
            "The file must exist.")

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
            fileNamePED=33, pedStudy=NULL, fileNameGDS=NULL,
            batch=1, studyDF=NULL, listSamples=NULL,
            PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when fileNamePED is NULL and pedStudy is NULL", {

    error_message <- paste0("One of the parameter \'fineNamePED\' of ",
                        "\'pedStudy\' must be defined.")

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
                    fileNamePED=NULL, pedStudy=NULL, fileNameGDS=NULL,
                    batch=1, studyDF=NULL, listSamples=NULL,
                    PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when pedDF is missing mandatory column", {

    pedDF <- data.frame(Name.ID = c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID = c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Sample.Type = rep("Primary Tumor", 3),
                    Source = rep("Databank B", 3), stringsAsFactors = FALSE)

    error_message <- paste0("The PED study data frame is incomplete. ",
                            "One or more mandatory columns are missing.")

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
                    fileNamePED=NULL, pedStudy=pedDF, fileNameGDS=NULL,
                    batch=1, studyDF=NULL, listSamples=NULL,
                    PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when fileNameGDS is numerical value", {

    pedDF <- data.frame(Name.ID = c("Sample_01", "Sample_02", "Sample_03"),
                        Case.ID = c("Patient_h11", "Patient_h12", "Patient_h18"),
                        Diagnosis = rep("Cancer", 3),
                        Sample.Type = rep("Primary Tumor", 3),
                        Source = rep("Databank B", 3), stringsAsFactors = FALSE)

    error_message <- paste0("The \'fileNameGDS\' must be a character ",
                "string representing the GDS 1KG file. The file must exist.")

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
                        fileNamePED=NULL, pedStudy=pedDF, fileNameGDS=33,
                        batch=1, studyDF=NULL, listSamples=NULL,
                        PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when batch is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    pedDF <- data.frame(Name.ID = c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID = c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis = rep("Cancer", 3),
                    Sample.Type = rep("Primary Tumor", 3),
                    Source = rep("Databank B", 3), stringsAsFactors = FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
            fileNamePED=NULL, pedStudy=pedDF, fileNameGDS=gdsFile,
            batch="1", studyDF=NULL, listSamples=NULL,
            PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when batch is vector of numerics", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    pedDF <- data.frame(Name.ID = c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID = c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis = rep("Cancer", 3),
                    Sample.Type = rep("Primary Tumor", 3),
                    Source = rep("Databank B", 3), stringsAsFactors = FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
                fileNamePED=NULL, pedStudy=pedDF, fileNameGDS=gdsFile,
                batch=c(1,2), studyDF=NULL, listSamples=NULL,
                PATHSAMPLEGDS=NULL, verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when listSamples is vector of numerics", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    pedDF <- data.frame(Name.ID = c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID = c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis = rep("Cancer", 3),
                    Sample.Type = rep("Primary Tumor", 3),
                    Source = rep("Databank B", 3), stringsAsFactors = FALSE)

    error_message <- paste0("The \'listSamples\' must be a vector ",
                        "of character strings (1 entry or more) or NULL.")

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
            fileNamePED=NULL, pedStudy=pedDF, fileNameGDS=gdsFile,
            batch=1, studyDF=NULL, listSamples=c(1,2),
            PATHSAMPLEGDS=NULL, verbose=TRUE), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when listSamples is numeric", {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")

    pedDF <- data.frame(Name.ID = c("Sample_01", "Sample_02", "Sample_03"),
                Case.ID = c("Patient_h11", "Patient_h12", "Patient_h18"),
                Diagnosis = rep("Cancer", 3),
                Sample.Type = rep("Primary Tumor", 3),
                Source = rep("Databank B", 3), stringsAsFactors = FALSE)

    error_message <- paste0("The \'listSamples\' must be a vector ",
                            "of character strings (1 entry or more) or NULL.")

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
        fileNamePED=NULL, pedStudy=pedDF, fileNameGDS=gdsFile,
        batch=1, studyDF=NULL, listSamples=1,
        PATHSAMPLEGDS=NULL, verbose=TRUE), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when verbose is numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    pedDF <- data.frame(Name.ID = c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID = c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis = rep("Cancer", 3),
                    Sample.Type = rep("Primary Tumor", 3),
                    Source = rep("Databank B", 3), stringsAsFactors = FALSE)

    error_message <- "The \'verbose\' parameter must be a logical (TRUE or FALSE)."

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
            fileNamePED=NULL, pedStudy=pedDF, fileNameGDS=gdsFile,
            batch=1, studyDF=NULL, listSamples=NULL,
            PATHSAMPLEGDS=NULL, verbose=22), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when both  fileNamePED and pedStudy are defined", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    pedDF <- data.frame(Name.ID = c("Sample_01", "Sample_02", "Sample_03"),
                        Case.ID = c("Patient_h11", "Patient_h12", "Patient_h18"),
                        Diagnosis = rep("Cancer", 3),
                        Sample.Type = rep("Primary Tumor", 3),
                        Source = rep("Databank B", 3), stringsAsFactors = FALSE)

    error_message <- paste0("Both \'fileNamePED\' and \'pedStudy\' parameters cannot be ",
                        "defined at the same time.")

    expect_error(createStudy2GDS1KG(PATHGENO=file.path("data", "sampleGeno"),
                fileNamePED=gdsFile, pedStudy=pedDF, fileNameGDS=gdsFile,
                batch=1, studyDF=NULL, listSamples=NULL,
                PATHSAMPLEGDS=NULL, verbose=22), error_message, fixed=TRUE)
})


#############################################################################
### Tests computePoolSyntheticAncestryGr() results
#############################################################################

context("computePoolSyntheticAncestryGr() results")


test_that(paste0("computePoolSyntheticAncestryGr() must return error when gds is numeric value"), {

    error_message <- "The \'gds\' must be an object of class \'gds.class\'"

    expect_error(computePoolSyntheticAncestryGr(gds=33, gdsSample="test.gds",
            sampleRM=c("1", "2"), study.id.syn="Synthetic", np=1L,
            spRef = "test", eigen.cnt = 15L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when gdsSample is character string"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- "The \'gdsSample\' must be an object of class \'gds.class\'"

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample="test.gds",
                sampleRM=c("1", "2"), study.id.syn="Synthetic", np=1L,
                spRef = "test", eigen.cnt = 15L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when sampleRM is vector of numeric values"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- "The \'sampleRM\' parameter must be a vector of character strings."

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
            sampleRM=c(1, 2), study.id.syn="Synthetic", np=1L,
            spRef="test", algorithm="exact", eigen.cnt=32L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when study.id.syn is numeric value"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- "The \'study.id.syn\' parameter must be a character string."

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
        sampleRM=c("Sample01", "Sample02"), study.id.syn=11, np=1L,
        spRef="test", algorithm="exact", eigen.cnt=32L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when listCatPop is numeric value"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- paste0("The \'listCatPop\' parameter must be a vector of ",
                                "character strings.")

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
        sampleRM=c("Sample01", "Sample02"), study.id.syn="Test", np=1L,
        listCatPop=11,
        spRef="test", algorithm="exact", eigen.cnt=32L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when np is character string"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
                sampleRM=c("1", "2"), study.id.syn="Synthetic", np="12",
                spRef = "test", eigen.cnt = 15L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when np is zero"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
                    sampleRM=c("1", "2"), study.id.syn="Synthetic", np=0L,
                    spRef = "test", eigen.cnt = 15L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when kList is a vector with zero"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- paste0("The \'kList\' parameter must be a vector of positive ",
                                "integers.")

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
        sampleRM=c("1", "2"), study.id.syn="Synthetic", np=1L, kList=c(0, 1),
        pcaList=c(1, 2), spRef = "test", eigen.cnt = 15L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when pcaList is a vector with zero"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- paste0("The \'pcaList\' parameter must be a ",
                                "vector of positive integers.")

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
        sampleRM=c("1", "2"), study.id.syn="Synthetic", np=1L,
        pcaList=c(0, 1, 2), spRef = "test", eigen.cnt = 15L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when algorithm is zero"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
                    sampleRM=c("1", "2"), study.id.syn="Synthetic", np=1L,
                    spRef="test", algorithm=22, eigen.cnt=15L), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when algorithm is not in the list of choices"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
                sampleRM=c("1", "2"), study.id.syn="Synthetic", np=1L,
                spRef="test", algorithm="Hello", eigen.cnt=15L))
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when eigen.cnt is character string"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- "The \'eigen.cnt\' parameter must be a single integer."

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
                sampleRM=c("1", "2"), study.id.syn="Synthetic", np=1L,
                spRef="test", algorithm="exact", eigen.cnt="15"), error_message)
})


test_that(paste0("computePoolSyntheticAncestryGr() must return error when missing.rate is character string"), {

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")
    gds1KG <- openfn.gds(gdsFIle)
    withr::defer(closefn.gds(gds1KG), envir=parent.frame())

    error_message <- paste0("The \'missing.rate\' parameter must be a single positive ",
                                    "numeric between zero and one or NaN.")

    expect_error(computePoolSyntheticAncestryGr(gds=gds1KG, gdsSample=gds1KG,
        sampleRM=c("1", "2"), study.id.syn="Synthetic", np=1L, spRef="test",
        algorithm="exact", eigen.cnt=15L, missing.rate="0.02"), error_message)
})


#############################################################################
### Tests computeAncestryFromSyntheticFile() results
#############################################################################


context("computeAncestryFromSyntheticFile() results")


test_that(paste0("computeAncestryFromSyntheticFile() must return error when gds is character string"), {

    gdsFile <- test_path("fixtures", "GDS_Sample_with_study_demo.gds")
    gdsSample <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    error_message <- "The \'gds\' must be an object of class \'gds.class\'"

    expect_error(computeAncestryFromSyntheticFile(gds="test.gds",
        gdsSample=gdsFile, listFiles, sample.ana.id, spRef, study.id.syn, np=1L,
        listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=NaN), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when gdsSample is character string"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'gdsSample\' must be an object of class \'gds.class\'"

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF,
        gdsSample="sample.gds", listFiles, sample.ana.id, spRef, study.id.syn,
        np=1L,  listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=NaN), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when study.id.syn is integer"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'study.id.syn\' parameter must be a character string."

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn=12L, np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=NaN), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when algorithm is numeric"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm=23,
        eigen.cnt=32L,  missing.rate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when np is negative"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn="Synthetic", np=-1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when listCatPop is numeric"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'listCatPop\' parameter must be a vector of ",
                                    "character strings.")

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn="Synthetic", np=1L, listCatPop=c(1, 2, 3),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when missing.rate is negative"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'missing.rate\' must be a single ",
                        "numeric positive value between 0 and 1 or NaN.")

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=-0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when fieldPopIn1KG is numeric"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'fieldPopIn1KG\' parameter must be a ",
                                "character string.")

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG=22, fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when fieldPopInfAnc is vector of strings"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'fieldPopInfAnc\' parameter must be a ",
                                    "character string.")

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="test", fieldPopInfAnc=c("SuperPop", "test"),
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when kList has one negative numeric"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'kList\' parameter must be a vector of ",
                                "positive integers.")

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="test", fieldPopInfAnc="SuperPop",
        kList=c(1, 2, -3, 4), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when pcaList has one negative numeric"), {

    gdsFile <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(gdsFile)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'pcaList\' parameter must be a vector of ",
                                "positive integers.")

    expect_error(computeAncestryFromSyntheticFile(gds=gdsF, gdsSample=gdsF,
        listFiles=gdsFile, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        study.id.syn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="test", fieldPopInfAnc="SuperPop",
        kList=c(1, 2, 3, 4), pcaList=c(2, -15, 1), algorithm="exact",
        eigen.cnt=32L,  missing.rate=0.2), error_message)
})


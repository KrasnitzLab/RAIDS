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

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gdsProfile=fileGDS, listPCA=list(),
                    currentProfile="sample1",
                    np="test"), error_message)
})

test_that("projectSample2PCA() must return error when np is negative integer", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gdsProfile=fileGDS, listPCA=list(),
                                    currentProfile="sample1",
                                    np=-1L), error_message)
})

test_that("projectSample2PCA() must return error when np is zero", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(projectSample2PCA(gdsProfile=fileGDS, listPCA=list(),
                                    currentProfile="sample1",
                                    np=0L), error_message)
})

test_that("projectSample2PCA() must return error when currentProfile is number", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    error_message <- paste0("The \'currentProfile\' ",
                                "parameter must be a character string.")

    expect_error(projectSample2PCA(gdsProfile=fileGDS, listPCA=list(),
                                        currentProfile=101,
                                        np=1L), error_message)
})


#############################################################################
### Tests appendStudy2GDS1KG() results
#############################################################################


context("appendStudy2GDS1KG() results")

test_that("appendStudy2GDS1KG() must return error when pathGeno is a numeric", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                            study.desc="Pancreatic study", study.platform="WES",
                            stringsAsFactors=FALSE)

    error_message <- paste0("The \'pathGeno\' must be a character string ",
                            "representing a path. The path must exist.")

    expect_error(appendStudy2GDS1KG(pathGeno=22, filePedRDS=sampleRDS,
        fileNameGDS=fileGDS, batch=2,
        studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose="TRUE"), error_message, fixed=TRUE)
})


test_that("appendStudy2GDS1KG() must return error when filePedRDS is numeric", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study", study.platform="WES",
                        stringsAsFactors=FALSE)

    error_message <- paste0("The \'filePedRDS\' must be a character string ",
        "representing the RDS Sample information file. The file must exist.")

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
        filePedRDS=21, fileNameGDS=fileGDS, batch=1,
        studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when fileNameGDS is numeric", {

    fileGDS <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- paste0("The \'fileNameGDS\' must be a character string ",
                        "representing the GDS 1KG file. The file must exist.")

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
        filePedRDS=sampleRDS, fileNameGDS=33, batch=1,
        studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when batch is a vector of numerics", {

    fileGDS <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES", study.desc="Pancreatic",
                                study.platform="WES", stringsAsFactors=FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
        filePedRDS=sampleRDS, fileNameGDS=fileGDS, batch=c(1,2),
        studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose=TRUE), error_message)
})

test_that("appendStudy2GDS1KG() must return error when batch is a character string", {

    fileGDS <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures",  "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
        filePedRDS=sampleRDS, fileNameGDS=fileGDS, batch="2",
        studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when studyDF is missing mandatory column", {

    fileGDS <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures",  "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.descption="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- paste0("The \'studyDF\' must be a data.frame and contain ",
        "those 3 columns: \'study.id\', \'study.desc\' and \'study.platform\'.")

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
            filePedRDS=sampleRDS, fileNameGDS=fileGDS, batch=1,
            studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
            genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("appendStudy2GDS1KG() must return error when listSamples is a numeric", {

    fileGDS <- test_path("fixtures",  "1KG_Test.gds")
    sampleRDS <- test_path("fixtures",  "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study", study.platform="WES",
                        stringsAsFactors=FALSE)

    error_message <- paste0("The \'listSamples\' must be a vector ",
        "of character strings (1 entry or more) or NULL.")

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
        filePedRDS=sampleRDS, fileNameGDS=fileGDS, batch=2,
        studyDF=studyInfo, listSamples=33, pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose=FALSE), error_message, fixed=TRUE)
})

test_that("appendStudy2GDS1KG() must return error when genoSource is a numeric", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                            study.desc="Pancreatic study", study.platform="WES",
                            stringsAsFactors=FALSE)

    error_message <- 'The \'genoSource\' parameter must be a character string.'

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
            filePedRDS=sampleRDS, fileNameGDS=fileGDS, batch=2,
            studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
            genoSource=3, verbose=TRUE), error_message, fixed=TRUE)
})


test_that("appendStudy2GDS1KG() must return error when verbose is a character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study", study.platform="WES",
                    stringsAsFactors=FALSE)

    error_message <- 'The \'verbose\' parameter must be a logical (TRUE or FALSE).'

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
        filePedRDS=sampleRDS, fileNameGDS=fileGDS, batch=2,
        studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose="TRUE"), error_message, fixed=TRUE)
})


test_that("appendStudy2GDS1KG() must return error when genoSource not in list of choices", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                            study.desc="Pancreatic study", study.platform="WES",
                            stringsAsFactors=FALSE)

    expect_error(appendStudy2GDS1KG(pathGeno=test_path("fixtures"),
        filePedRDS=sampleRDS, fileNameGDS=fileGDS, batch=2,
        studyDF=studyInfo, listSamples=NULL, pathProfileGDS=NULL,
        genoSource="TOP-pileup", verbose=TRUE))
})

#############################################################################
### Tests pruningSample() results
#############################################################################


context("pruningSample() results")


test_that("pruningSample() must return error when gdsReference is a character string", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    error_message <- "The \'gdsReference\' must be an object of class \'gds.class\'."

    expect_error(pruningSample(gdsReference=fileGDS, method="corr", currentProfile="test",
        studyID="test", listSNP=NULL, slideWindowMaxBP=5e5, thresholdLD=sqrt(0.1),
        np=1, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=".",
        outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when keepPrunedGDS is a character string", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures",  "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- 'The \'keepPrunedGDS\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile="test",
        studyID="test", listSNP=NULL, slideWindowMaxBP=5e5, thresholdLD=sqrt(0.1),
        np=1, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS="YES",
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=".",
        outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when keepFile is a character string", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

        gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- 'The \'keepFile\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile="1KG_Test",
        studyID="test", listSNP=NULL, slideWindowMaxBP=5e5, thresholdLD=sqrt(0.1),
        np=1, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=TRUE,
        pathProfileGDS=dataDir, keepFile="NO", pathPrunedGDS=".",
        outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when np is a character string", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'np\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile="test",
        studyID="test", listSNP=NULL, slideWindowMaxBP=5e5, thresholdLD=sqrt(0.1),
        np="1", verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=".", outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when slideWindowMaxBP is a character string", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'slideWindowMaxBP\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile="test",
        studyID="test", listSNP=NULL, slideWindowMaxBP="4", thresholdLD=sqrt(0.1),
        np=1, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=".", outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when thresholdLD is a character string", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'thresholdLD\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile="test",
        studyID="test", listSNP=NULL, slideWindowMaxBP=4, thresholdLD="3",
        np=1, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=".", outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when thresholdLD is a vector of numerics", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'thresholdLD\' parameter must be a single positive numeric value."

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile="test",
        studyID="test", listSNP=NULL, slideWindowMaxBP=4, thresholdLD=c(3,3),
        np=1, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=".", outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when method is a numeric", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())
    error_message <- "The \'method\' parameter must be a character string."

    expect_error(pruningSample(gdsReference=gdsF, method=3, currentProfile="test",
        studyID="test", listSNP=NULL, slideWindowMaxBP=50000L, thresholdLD=sqrt(0.1),
        np=1, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=".", outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when method is not in the list of choices", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    expect_error(pruningSample(gdsReference=gdsF, method="test", currentProfile="test",
        studyID="test", listSNP=NULL, slideWindowMaxBP=50000L, thresholdLD=sqrt(0.1),
        np=1, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=".", outPrefix="pruned"))
})


test_that("pruningSample() must return error when currentProfile is a numeric", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'currentProfile\' parameter must be a character string."

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile=2,
        studyID="test", listSNP=NULL, slideWindowMaxBP=50000L,
        thresholdLD=sqrt(0.1), np=1, verbose=FALSE, chr=NULL,
        superPopMinAF=NULL, keepPrunedGDS=FALSE, pathProfileGDS=dataDir,
        keepFile=FALSE, pathPrunedGDS=".", outPrefix="pruned"),
        error_message, fixed=TRUE)
})



test_that("pruningSample() must return error when pathPrunedGDS is a numeric", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'pathPrunedGDS\' parameter must be ",
        "a character string representing an existing directory.")

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile="Sample2",
        studyID="test", listSNP=NULL, slideWindowMaxBP=50000L, thresholdLD=sqrt(0.1),
        np=1L, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=dataDir, keepFile=FALSE, pathPrunedGDS=2, outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when pathProfileGDS is a numeric", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'pathProfileGDS\' parameter must be a ",
                    "character string representing an existing directory.")

    expect_error(pruningSample(gdsReference=gdsF, method="corr", currentProfile="Sample2",
            studyID="test", listSNP=NULL, slideWindowMaxBP=50000L, thresholdLD=sqrt(0.1),
            np=1L, verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
            pathProfileGDS=33, keepFile=FALSE, pathPrunedGDS=dataDir, outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when pathProfileGDS is a non existing directory", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'pathProfileGDS\' parameter must be a ",
                    "character string representing an existing directory.")

    expect_error(pruningSample(gdsReference=gdsF, method="corr",
        currentProfile="Sample2", studyID="test", listSNP=NULL,
        slideWindowMaxBP=50000L, thresholdLD=sqrt(0.1), np=1L,
        verbose=FALSE, chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=paste0(dataDir, "_NOT_EXISTING_DIRECTORY"),
        keepFile=FALSE, pathPrunedGDS=dataDir, outPrefix="pruned"),
        error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when verbose is a character string", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- 'The \'verbose\' parameter must be a logical (TRUE or FALSE).'

    expect_error(pruningSample(gdsReference=gdsF, method="corr",
        currentProfile="Sample2", studyID="test", listSNP=NULL,
        slideWindowMaxBP=50000L, thresholdLD=sqrt(0.1), np=1L, verbose="HI",
        chr=NULL, superPopMinAF=NULL, keepPrunedGDS=FALSE,
        pathProfileGDS=pdataDir, keepFile=FALSE, pathPrunedGDS=dataDir,
        outPrefix="pruned"), error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when GDS Sample file does not exist", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    sampleRDS <- test_path("fixtures", "Sample_Info_Test.RDS")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    fileProfileGDS <- file.path(dataDir, "A_File_That_DOES_NOT_EXIST.gds")
    error_message <- paste0("The Profile GDS file \'", fileProfileGDS,
                                    " does not exist.")

    expect_error(pruningSample(gdsReference=gdsF, method="corr",
        currentProfile="A_File_That_DOES_NOT_EXIST",
        studyID="test", listSNP=NULL, slideWindowMaxBP=50000L,
        thresholdLD=sqrt(0.1), np=1L, verbose=FALSE, chr=NULL,
        superPopMinAF=NULL, keepPrunedGDS=FALSE, pathProfileGDS=dataDir,
        keepFile=FALSE, pathPrunedGDS=dataDir, outPrefix="pruned"),
        error_message, fixed=TRUE)
})


test_that("pruningSample() must return error when no SNV left after filtering", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                    study.platform="PLATFORM",
                    stringsAsFactors=FALSE)

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())
    dataDirSample <- test_path("fixtures/sampleGDSforPruning")

    file.copy(file.path(dataDirSample, "ex1_demoForPruning.gds"),
              file.path(dataDirSample, "ex1.gds"))
    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                 envir=parent.frame())

    error_message <- paste0("In pruningSample, the sample ex1 ",
                                    "doesn't have SNPs after filters")

    expect_error(pruningSample(gdsReference=gdsF, method="corr",
        currentProfile="ex1", studyID=studyDF$study.id,
        listSNP=NULL, slideWindowMaxBP=50000L,
        thresholdLD=sqrt(0.1), np=1L, verbose=FALSE, chr=22,
        superPopMinAF=0.41, keepPrunedGDS=TRUE, pathProfileGDS=dataDirSample,
        keepFile=TRUE, pathPrunedGDS=dataDirSample, outPrefix="prunedTest"),
        error_message, fixed=TRUE)
})




test_that("pruningSample() must return error when the study is not found", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                          study.platform="PLATFORM",
                          stringsAsFactors=FALSE)

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    dataDirSample <- test_path("fixtures/sampleGDSforPruning")

    file.copy(file.path(dataDirSample, "ex1_demoForPruning.gds"),
              file.path(dataDirSample, "ex1.gds"))

    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                 envir=parent.frame())

    error_message <- paste0("In pruningSample the profile \'ex1\'",
                        " doesn't exists for the study \'demo\'\n")

    expect_error(pruningSample(gdsReference=gdsF, method="corr",
        currentProfile="ex1", studyID="demo",
        listSNP=NULL, slideWindowMaxBP=50000L,
        thresholdLD=sqrt(0.1), np=1L, verbose=FALSE, chr=22,
        superPopMinAF=0.41, keepPrunedGDS=TRUE, pathProfileGDS=dataDirSample,
        keepFile=TRUE, pathPrunedGDS=dataDirSample, outPrefix="prunedTest"),
        error_message, fixed=TRUE)
})


test_that("pruningSample() must return expect result", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                          study.platform="PLATFORM",
                          stringsAsFactors=FALSE)

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    dataDirSample <- test_path("fixtures/sampleGDSforPruning")

    file.copy(file.path(dataDirSample, "ex1_demoForPruning.gds"),
                      file.path(dataDirSample, "ex1.gds"))
    withr::defer((unlink(file.path(dataDirSample, "prunedTest.Obj.rds"))),
                        envir=parent.frame())
    withr::defer((unlink(file.path(dataDirSample, "prunedTest.rds"))),
                        envir=parent.frame())
    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                 envir=parent.frame())

    result <- pruningSample(gdsReference=gdsF, method="corr",
                currentProfile="ex1", studyID=studyDF$study.id,
                listSNP=NULL, slideWindowMaxBP=50000L,
                thresholdLD=sqrt(0.1), np=1L, verbose=FALSE, chr=NULL,
                superPopMinAF=NULL, keepPrunedGDS=TRUE,
                pathProfileGDS=dataDirSample, keepFile=TRUE,
                pathPrunedGDS=dataDirSample, outPrefix="prunedTest")

    expect_equal(result, 0L)
    expect_true(file.exists(file.path(dataDirSample, "prunedTest.Obj.rds")))
    expect_true(file.exists(file.path(dataDirSample, "prunedTest.rds")))
})


#############################################################################
### Tests add1KG2SampleGDS() results
#############################################################################


context("add1KG2SampleGDS() results")


test_that("add1KG2SampleGDS() must return error when gdsReference is a character string", {

    error_message <- "The \'gdsReference\' must be an object of class \'gds.class\'."

    expect_error(add1KG2SampleGDS(gdsReference="toto.gds", fileProfileGDS="sample.gds",
        currentProfile="sample", studyID="TCGA"), error_message, fixed=TRUE)
})


test_that("add1KG2SampleGDS() must return error when fileProfileGDS is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'fileProfileGDS\' must be a character ",
            "string representing the GDS Sample file. The file must exist.")

    expect_error(add1KG2SampleGDS(gdsReference=gdsF, fileProfileGDS=33,
        currentProfile="sample", studyID="TCGA"), error_message, fixed=TRUE)
})

test_that("add1KG2SampleGDS() must return error when currentProfile is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'currentProfile\' must be a character string."

    expect_error(add1KG2SampleGDS(gdsReference=gdsF, fileProfileGDS=fileGDS,
            currentProfile=33, studyID="TCGA"), error_message, fixed=TRUE)
})


test_that("add1KG2SampleGDS() must return error when studyID is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'studyID\' must be a character string."

    expect_error(add1KG2SampleGDS(gdsReference=gdsF, fileProfileGDS=fileGDS,
        currentProfile="Test", studyID=22), error_message, fixed=TRUE)
})


test_that("add1KG2SampleGDS() must return expect result", {

    dataDir <- test_path("fixtures")

    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                          study.platform="PLATFORM",
                          stringsAsFactors=FALSE)

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    dataDirSample <- test_path("fixtures/sampleGDSforAddingGenotype")


    file.copy(file.path(dataDirSample, "ex1_demoForAddGenotype.gds"),
              file.path(dataDirSample, "ex1.gds"))
    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                 envir=parent.frame())

    result <- add1KG2SampleGDS(gdsReference=gdsF,
                fileProfileGDS=file.path(dataDirSample, "ex1.gds"),
                currentProfile=c("ex1"), studyID=studyDF$study.id)

    expect_equal(result, 0L)

    content <- openfn.gds(file.path(dataDirSample, "ex1.gds"))
    withr::defer((gdsfmt::closefn.gds(content)), envir=parent.frame())

    expect_true(gdsfmt::exist.gdsn(content, "pruned.study"))
    expect_true(gdsfmt::exist.gdsn(content, "geno.ref"))
    expect_true(gdsfmt::exist.gdsn(content, "sample.id"))
    expect_true(gdsfmt::exist.gdsn(content, "snp.id"))
    expect_true(gdsfmt::exist.gdsn(content, "snp.position"))
    expect_true(gdsfmt::exist.gdsn(content, "snp.index"))
    expect_true(gdsfmt::exist.gdsn(content, "genotype"))
    expect_true(gdsfmt::exist.gdsn(content, "lap"))
})


#############################################################################
### Tests computePCARefSample() results
#############################################################################


context("computePCARefSample() results")


test_that("computePCARefSample() must return error when gdsSample isnumeric value", {

    error_message <- "The \'gdsSample\' must be an object of class \'gds.class\'."

    expect_error(computePCARefSample(gdsSample="test.gds", name.id="test",
                            studyIDRef="Ref.1KG", np=1L, algorithm="exact",
                            eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when name.id is numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'name.id\' parameter must be a single character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id=22,
                            studyIDRef="Ref.1KG", np=1L, algorithm="exact",
                            eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when studyIDRef is numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'studyIDRef\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="Synthetic",
                            studyIDRef=33, np=1L, algorithm="exact",
                            eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when np is a character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
                        studyIDRef="Ref.1KG",
                        np="1", algorithm="exact",
                        eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when np is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
                                studyIDRef="Ref.1KG", np=1L, algorithm=33,
                                eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when algorithm is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
            studyIDRef="Ref.1KG", np=1L, algorithm=33,
            eigen.cnt=32L), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when algorithm is not a valid choice", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
            studyIDRef="Ref.1KG", np=1L, algorithm="sun", eigen.cnt=32L))
})


test_that("computePCARefSample() must return error when eigen.cnt is a string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'eigen.cnt\' parameter must be a single integer."

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
        studyIDRef="Ref.1KG", np=1L, algorithm="sun", eigen.cnt="32L"),
        error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when missingRate is negative value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'missingRate\' must be a single numeric ",
                                "positive value between 0 and 1 or NaN.")

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
        studyIDRef="Ref.1KG", np=1L, algorithm="sun", eigen.cnt=32L,
        missingRate=-0.02), error_message, fixed=TRUE)
})


test_that("computePCARefSample() must return error when algorithm is not in the list", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("\'arg\' should be one of \"exact\", ",
                                        "\"randomized\"")

    expect_error(computePCARefSample(gdsSample=gdsF, name.id="TCGA",
        studyIDRef="Ref.1KG", np=1L, algorithm="TITI", eigen.cnt=32L,
        missingRate=0.02), error_message, fixed=TRUE)
})


#############################################################################
### Tests addStudy1Kg() results
#############################################################################


context("addStudy1Kg() results")


test_that("addStudy1Kg() must return error when gdsReference is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'gdsReference\' must be an object of class \'gds.class\'."

    expect_error(addStudy1Kg(gdsReference=33, fileProfileGDS=fileGDS),
                                error_message, fixed=TRUE)
})


test_that("addStudy1Kg() must return error when fileProfileGDS is a numeric value", {

    ## Create a temporary GDS file in an test directory
    fileGDS <- test_path("fixtures", "GDS_TEMP_processStudy_101.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_Sample_file(fileGDS)

    error_message <- paste0("The \'fileProfileGDS\' must be a character ",
        "string representing the GDS Sample file. The file must exist.")

    expect_error(addStudy1Kg(gdsReference=GDS_file_tmp, fileProfileGDS=33),
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
                        study.platform=c("Panel"), stringsAsFactors=FALSE)

    add.gdsn(GDS_file_Sample, "study.list", study.list)

    study.annot <- data.frame(data.id=c("TOTO1"), case.id=c("TOTO1"),
                sample.type=c("Study"), diagnosis=c("Study"),
                source=rep("IGSR"), study.id=c("Study"),
                stringsAsFactors=FALSE)

    add.gdsn(GDS_file_Sample, "study.annot", study.annot)

    sync.gds(GDS_file_Sample)

    closefn.gds(GDS_file_Sample)
    withr::defer((unlink(gdsFileSample, force=TRUE)), envir=parent.frame())

    result0 <- addStudy1Kg(gdsReference=GDS_file_tmp_1KG,
                                fileProfileGDS=gdsFileSample)

    gds_sample_file <- openfn.gds(gdsFileSample, readonly=TRUE)

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
                        study.platform=c("Panel"), stringsAsFactors=FALSE)

    add.gdsn(GDS_file_Sample, "study.list", study.list)

    study.annot <- data.frame(data.id=c("TOTO1"), case.id=c("TOTO1"),
                             sample.type=c("Study"), diagnosis=c("Study"),
                             source=rep("IGSR"), study.id=c("Ref.1KG"),
                             stringsAsFactors=FALSE)

    add.gdsn(GDS_file_Sample, "study.annot", study.annot)

    sync.gds(GDS_file_Sample)

    closefn.gds(GDS_file_Sample)
    withr::defer((unlink(gdsFileSample)), envir=parent.frame())

    result0 <- addStudy1Kg(gdsReference=GDS_file_tmp_1KG, fileProfileGDS=gdsFileSample)

    gds_sample_file <- openfn.gds(gdsFileSample, readonly=TRUE)

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
### Tests createStudy2GDS1KG() results
#############################################################################

context("createStudy2GDS1KG() results")


test_that(paste0("createStudy2GDS1KG() must return error when filePedRDS is",
            " a numeric value and pedStudy is NULL"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")

    error_message <- paste0("The \'filePedRDS\' must be a character string ",
            "representing the RDS Sample information file. ",
            "The file must exist.")

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
            filePedRDS=33, pedStudy=NULL, fileNameGDS=NULL,
            batch=1, studyDF=NULL, listProfiles=NULL,
            pathProfileGDS=NULL, genoSource="snp-pileup", verbose=TRUE),
            error_message)
})


test_that("createStudy2GDS1KG() must return error when filePedRDS is NULL and pedStudy is NULL", {

    error_message <- paste0("One of the parameter \'fineNamePED\' of ",
                        "\'pedStudy\' must be defined.")

    expect_error(createStudy2GDS1KG(pathGeno=file.path("data", "sampleGeno"),
                    filePedRDS=NULL, pedStudy=NULL, fileNameGDS=NULL,
                    batch=1, studyDF=NULL, listProfiles=NULL,
                    pathProfileGDS=NULL,
                    genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when pedDF is missing mandatory column", {

    dataDir <- system.file("extdata/tests", package="RAIDS")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Sample.Type=rep("Primary Tumor", 3),
                    Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    error_message <- paste0("The PED study data frame is incomplete. ",
                            "One or more mandatory columns are missing.")

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
                    filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=NULL,
                    batch=1, studyDF=NULL, listProfiles=NULL,
                    pathProfileGDS=NULL,
                    genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when fileNameGDS is numerical value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis=rep("Cancer", 3),
                    Sample.Type=rep("Primary Tumor", 3),
                    Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    error_message <- paste0("The \'fileNameGDS\' must be a character ",
                "string representing the GDS 1KG file. The file must exist.")

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
                        filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=33,
                        batch=1, studyDF=NULL, listProfiles=NULL,
                        pathProfileGDS=NULL,
                        genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when batch is character string", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis=rep("Cancer", 3),
                    Sample.Type=rep("Primary Tumor", 3),
                    Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
            filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=fileGDS,
            batch="1", studyDF=NULL, listProfiles=NULL,
            pathProfileGDS=NULL,
            genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when batch is vector of numerics", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis=rep("Cancer", 3),
                    Sample.Type=rep("Primary Tumor", 3),
                    Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    error_message <- "The \'batch\' must be a single integer."

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
                filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=fileGDS,
                batch=c(1,2), studyDF=NULL, listProfiles=NULL,
                pathProfileGDS=NULL,
                genoSource="snp-pileup", verbose=TRUE), error_message)
})


test_that("createStudy2GDS1KG() must return error when listSamples is vector of numerics", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis=rep("Cancer", 3),
                    Sample.Type=rep("Primary Tumor", 3),
                    Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                          study.platform="PLATFORM", stringsAsFactors=FALSE)

    error_message <- paste0("The \'listProfiles\' must be a vector ",
                        "of character strings (1 entry or more) or NULL.")

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
            filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=fileGDS,
            batch=1, studyDF=studyDF, listProfiles=c(1,2),
            pathProfileGDS=NULL,
            genoSource="snp-pileup", verbose=TRUE), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when listProfiles is numeric", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                Diagnosis=rep("Cancer", 3),
                Sample.Type=rep("Primary Tumor", 3),
                Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                        study.platform="PLATFORM", stringsAsFactors=FALSE)

    error_message <- paste0("The \'listProfiles\' must be a vector ",
                            "of character strings (1 entry or more) or NULL.")

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
        filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=fileGDS,
        batch=1, studyDF=studyDF, listProfiles=1,
        pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose=TRUE), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when studyDF is missing column", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis=rep("Cancer", 3),
                    Sample.Type=rep("Primary Tumor", 3),
                    Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                         stringsAsFactors=FALSE)

    error_message <- paste0("The study data frame \'studyDF\' is incomplete. ",
        "One or more mandatory columns are missing. The mandatory ",
        "columns are: \'study.id\', \'study.desc\', \'study.platform\'.")

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
        filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=fileGDS,
        batch=1, studyDF=studyDF, listProfiles=1,
        pathProfileGDS=NULL,
        genoSource="snp-pileup", verbose=TRUE), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when verbose is numeric", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                    Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                    Diagnosis=rep("Cancer", 3),
                    Sample.Type=rep("Primary Tumor", 3),
                    Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                    study.platform="PLATFORM", stringsAsFactors=FALSE)

    error_message <- "The \'verbose\' parameter must be a logical (TRUE or FALSE)."

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
            filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=fileGDS,
            batch=1, studyDF=studyDF, listProfiles=NULL,
            pathProfileGDS=dataDir,
            genoSource="snp-pileup", verbose=22), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when pathProfileGDS is numeric", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(dataDir, "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                Diagnosis=rep("Cancer", 3),
                Sample.Type=rep("Primary Tumor", 3),
                Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                study.platform="PLATFORM", stringsAsFactors=FALSE)

    error_message <- paste0("The \'pathProfileGDS\' must be a character ",
        "string representing the path where the Profile GDS files ",
        "will be generated.")

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
        filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=gdsFile,
        batch=1, studyDF=studyDF, listProfiles=NULL,
        pathProfileGDS=33,
        genoSource="snp-pileup", verbose=FALSE), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when both filePedRDS and pedStudy are defined", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                        Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                        Diagnosis=rep("Cancer", 3),
                        Sample.Type=rep("Primary Tumor", 3),
                        Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                    study.platform="PLATFORM", stringsAsFactors=FALSE)

    error_message <- paste0("Both \'filePedRDS\' and \'pedStudy\' parameters ",
                        "cannot be defined at the same time.")

    expect_error(createStudy2GDS1KG(pathGeno=file.path("data", "sampleGeno"),
        filePedRDS=fileGDS, pedStudy=pedDF, fileNameGDS=fileGDS,
        batch=1, studyDF=studyDF, listProfiles=NULL,
        pathProfileGDS=dataDir,
        genoSource="snp-pileup", verbose=TRUE), error_message, fixed=TRUE)
})


test_that("createStudy2GDS1KG() must return error when genoSource not in list of choices", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    pedDF <- data.frame(Name.ID=c("Sample_01", "Sample_02", "Sample_03"),
                        Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                        Diagnosis=rep("Cancer", 3),
                        Sample.Type=rep("Primary Tumor", 3),
                        Source=rep("Databank B", 3), stringsAsFactors=FALSE)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                          study.platform="PLATFORM", stringsAsFactors=FALSE)

    expect_error(createStudy2GDS1KG(pathGeno=dataDir,
        filePedRDS=NULL, pedStudy=pedDF, fileNameGDS=fileGDS,
        batch=1, studyDF=studyDF, listProfiles=NULL,
        pathProfileGDS=dataDir, genoSource="snp-CANADA", verbose=TRUE))
})


test_that("createStudy2GDS1KG() must return expected results when all parameters ok", {

    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")

    withr::defer((unlink(file.path(dataDir, "ex1.gds"))), envir=parent.frame())


    pedDF <- data.frame(Name.ID=c("ex1", "ex2", "ex3"),
                Case.ID=c("Patient_h11", "Patient_h12", "Patient_h18"),
                Diagnosis=rep("Cancer", 3),
                Sample.Type=rep("Primary Tumor", 3),
                Source=rep("Databank B", 3), stringsAsFactors=FALSE)
    rownames(pedDF) <- pedDF$Name.ID

    studyDF <- data.frame(study.id="MYDATA",
                            study.desc="Description",
                            study.platform="PLATFORM",
                            stringsAsFactors=FALSE)

    result <- createStudy2GDS1KG(pathGeno=dataDir,
                pedStudy=pedDF, fileNameGDS=fileGDS,
                batch=1, studyDF=studyDF, listProfiles=c("ex1"),
                pathProfileGDS=dataDir,
                genoSource="snp-pileup", verbose=FALSE)

    expect_true(file.exists(file.path(dataDir, "ex1.gds")))
    expect_equal(result, 0L)
})


#############################################################################
### Tests computePoolSyntheticAncestryGr() results
#############################################################################

context("computePoolSyntheticAncestryGr() results")


test_that("computePoolSyntheticAncestryGr() must return error when gdsProfile is character string", {

    error_message <- "The \'gdsProfile\' must be an object of class \'gds.class\'"

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile="test.gds",
                sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L,
                spRef=spRef, eigenCount=15L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when sampleRM is vector of numeric values", {

    fileGDS <- test_path("fixtures", "GDS_Sample_with_study_demo.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    error_message <- "The \'sampleRM\' parameter must be a vector of character strings."

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
            sampleRM=c(1, 2), studyIDSyn="Synthetic", np=1L,
            spRef=spRef, algorithm="exact", eigenCount=32L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when studyIDSyn is numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- "The \'studyIDSyn\' parameter must be a character string."

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
        sampleRM=c("Sample01", "Sample02"), studyIDSyn=11, np=1L,
        spRef=spRef, algorithm="exact", eigenCount=32L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when listCatPop is numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- paste0("The \'listCatPop\' parameter must be a vector of ",
                                "character strings.")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
        sampleRM=c("Sample01", "Sample02"), studyIDSyn="Test", np=1L,
        listCatPop=11,
        spRef=spRef, algorithm="exact", eigenCount=32L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when np is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
                sampleRM=c("1", "2"), studyIDSyn="Synthetic", np="12",
                spRef=spRef, eigenCount=15L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when np is zero", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
                    sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=0L,
                    spRef=spRef, eigenCount=15L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when kList is a vector with zero", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- paste0("The \'kList\' parameter must be a vector of ",
                                "positive integers.")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
        sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L, kList=c(0, 1),
        pcaList=c(1, 2), spRef=spRef, eigenCount=15L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when pcaList is a vector with zero", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- paste0("The \'pcaList\' parameter must be a ",
                                "vector of positive integers.")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
        sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L,
        pcaList=c(0, 1, 2), spRef=spRef, eigenCount=15L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when algorithm is zero", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
                    sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L,
                    spRef=spRef, algorithm=22, eigenCount=15L), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when algorithm is not in the list of choices", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
                sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L,
                spRef=spRef, algorithm="Hello", eigenCount=15L))
})


test_that("computePoolSyntheticAncestryGr() must return error when eigenCount is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- "The \'eigenCount\' parameter must be a single integer."

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
                sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L,
                spRef=spRef, algorithm="exact", eigenCount="15"), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when missingRate is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- paste0("The \'missingRate\' parameter must be a single ",
                            "positive numeric between zero and one or NaN.")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
        sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L, spRef=spRef,
        algorithm="exact", eigenCount=15L, missingRate="0.02"), error_message)
})


test_that("computePoolSyntheticAncestryGr() must return error when verbose is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- paste0("The \'verbose\' parameter must be a ",
                                "logical (TRUE or FALSE).")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
        sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L, spRef=spRef,
        algorithm="exact", eigenCount=15L, missingRate=0.02, verbose="QC"),
        error_message, fixed=TRUE)
})


test_that("computePoolSyntheticAncestryGr() must return error when spRef is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())


    error_message <- paste0("The \'spRef\' parameter must be a vector of ",
                    "character strings with profile identifiers as names.")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
        sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L, spRef="HI",
        algorithm="exact", eigenCount=15L, missingRate=0.02, verbose=TRUE),
        error_message, fixed=TRUE)
})


test_that("computePoolSyntheticAncestryGr() must return error when fieldPopInfAnc is numeric", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    error_message <- paste0("The \'fieldPopInfAnc\' parameter must be a ",
                                "character string.")

    expect_error(computePoolSyntheticAncestryGr(gdsProfile=gdsSample,
        sampleRM=c("1", "2"), studyIDSyn="Synthetic", np=1L, spRef=spRef,
        fieldPopInfAnc=33, algorithm="exact", eigenCount=15L, missingRate=0.02,
        verbose=TRUE), error_message, fixed=TRUE)
})


#############################################################################
### Tests computeAncestryFromSyntheticFile() results
#############################################################################


context("computeAncestryFromSyntheticFile() results")


test_that("computeAncestryFromSyntheticFile() must return error when gdsReference is character string", {

    fileGDS <- test_path("fixtures", "GDS_Sample_with_study_demo.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsSample), envir=parent.frame())

    error_message <- "The \'gdsReference\' must be an object of class \'gds.class\'"

    expect_error(computeAncestryFromSyntheticFile(gdsReference="test.gds",
        gdsSample=fileGDS, listFiles, sample.ana.id, spRef, studyIDSyn, np=1L,
        listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=NaN), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when gdsSample is character string"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'gdsSample\' must be an object of class \'gds.class\'"

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF,
        gdsSample="sample.gds", listFiles, sample.ana.id, spRef, studyIDSyn,
        np=1L,  listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=NaN), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when studyIDSyn is integer"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'studyIDSyn\' parameter must be a character string."

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn=12L, np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=NaN), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when algorithm is numeric"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'algorithm\' parameter must be a character string."

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm=23,
        eigen.cnt=32L,  missingRate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when np is negative"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=-1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when listCatPop is numeric"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'listCatPop\' parameter must be a vector of ",
                                    "character strings.")

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=1L, listCatPop=c(1, 2, 3),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when missingRate is negative"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'missingRate\' must be a single ",
                        "numeric positive value between 0 and 1 or NaN.")

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=-0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when eigenCount is a vector of integers"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'eigenCount\' parameter must be a single integer."

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=c(2L, 3L),  missingRate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when fieldPopIn1KG is numeric"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'fieldPopIn1KG\' parameter must be a ",
                                "character string.")

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG=22, fieldPopInfAnc="SuperPop",
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when fieldPopInfAnc is vector of strings"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'fieldPopInfAnc\' parameter must be a ",
                                    "character string.")

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="test", fieldPopInfAnc=c("SuperPop", "test"),
        kList=seq(2, 15, 1), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when kList has one negative numeric"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'kList\' parameter must be a vector of ",
                                "positive integers.")

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="test", fieldPopInfAnc="SuperPop",
        kList=c(1, 2, -3, 4), pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=0.2), error_message)
})


test_that(paste0("computeAncestryFromSyntheticFile() must return error when pcaList has one negative numeric"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'pcaList\' parameter must be a vector of ",
                                "positive integers.")

    expect_error(computeAncestryFromSyntheticFile(gdsReference=gdsF, gdsSample=gdsF,
        listFiles=fileGDS, sample.ana.id="sample01", spRef=c("HC01", "HC03"),
        studyIDSyn="Synthetic", np=1L, listCatPop=c("EAS", "EUR", "AFR"),
        fieldPopIn1KG="test", fieldPopInfAnc="SuperPop",
        kList=c(1, 2, 3, 4), pcaList=c(2, -15, 1), algorithm="exact",
        eigen.cnt=32L,  missingRate=0.2), error_message)
})


#############################################################################
### Tests computePoolSyntheticAncestry() results
#############################################################################

context("computePoolSyntheticAncestry() results")


test_that("computePoolSyntheticAncestry() must return error when gdsReference is a character string", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    dataRefDemo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC",
                                                    "SampleD"),
                            pop.group=c("TSI", "TSI", "YRI", "YRI"),
                            superPop=c("EUR", "EUR", "AFR", "AFR"))

    error_message <- paste0("The \'gdsReference\' must be an object of ",
                                        "class \'gds.class\'.")

    expect_error(computePoolSyntheticAncestry(gdsReference="toto.gds", gdsSample=gdsF,
        sample.ana.id="test", dataRef=dataRefDemo, spRef="TODO",
        studyIDSyn="synthetic",
        np=1L, listCatPop="EUR", fieldPopIn1KG="superPop",
        fieldPopInfAnc="SuperPop", kList=seq(2, 15, 1),
        pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L, missingRate=0.025), error_message)
})


test_that(paste0("computePoolSyntheticAncestry() must return error ",
                    "when gds is a numerical value"), {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    dataRefDemo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC",
                                                         "SampleD"),
                            pop.group=c("TSI", "TSI", "YRI", "YRI"),
                            superPop=c("EUR", "EUR", "AFR", "AFR"))

    error_message <- paste0("The \'gdsReference\' must be an object of ",
                                    "class \'gds.class\'.")

    expect_error(computePoolSyntheticAncestry(gdsReference=31, gdsSample=gdsF,
            sample.ana.id="test", dataRef=dataRefDemo, spRef="TODO",
            studyIDSyn="synthetic",
            np=1L, listCatPop="EUR", fieldPopIn1KG="superPop",
            fieldPopInfAnc="SuperPop", kList=seq(2, 15, 1),
            pcaList=seq(2, 15, 1), algorithm="exact",
            eigen.cnt=32L, missingRate=0.025), error_message)
})


test_that("computePoolSyntheticAncestry() must return error when profileGDS is a numerical value", {

    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    dataRefDemo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC",
                                            "SampleD"),
                                pop.group=c("TSI", "TSI", "YRI", "YRI"),
                                superPop=c("EUR", "EUR", "AFR", "AFR"))

    error_message <- paste0("The \'profileGDS\' must be an object of ",
                                    "class \'gds.class\'.")

    expect_error(computePoolSyntheticAncestry(gdsReference=gdsF, gdsSample=33,
            sample.ana.id="test", dataRef=dataRefDemo, spRef="TODO",
            studyIDSyn="synthetic",
            np=1L, listCatPop="EUR", fieldPopIn1KG="superPop",
            fieldPopInfAnc="SuperPop", kList=seq(2, 15, 1),
            pcaList=seq(2, 15, 1), algorithm="exact",
            eigen.cnt=32L, missingRate=0.025), error_message)
})


test_that("computePoolSyntheticAncestry() must return error when dataRef is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'dataRef\' must be a data.frame object."

    expect_error(computePoolSyntheticAncestry(gdsReference=gdsF, gdsSample=gdsF,
        sample.ana.id="test", dataRef="test", spRef="TODO",
        studyIDSyn="synthetic",
        np=1L, listCatPop="EUR", fieldPopIn1KG="superPop",
        fieldPopInfAnc="SuperPop", kList=seq(2, 15, 1),
        pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L, missingRate=0.025), error_message)
})


test_that("computePoolSyntheticAncestry() must return error when studyIDSyn is numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    dataRefDemo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC", "SampleD"),
                              pop.group=c("TSI", "TSI", "YRI", "YRI"),
                              superPop=c("EUR", "EUR", "AFR", "AFR"))

    error_message <- "The \'studyIDSyn\' parameter must be a character string."

    expect_error(computePoolSyntheticAncestry(gdsReference=gdsF, gdsSample=gdsF,
        sample.ana.id="test", dataRef=dataRefDemo, spRef="TODO",
        studyIDSyn=33,
        np=1L, listCatPop="EUR", fieldPopIn1KG="SuperPop",
        fieldPopInfAnc="SuperPop", kList=seq(2, 15, 1),
        pcaList=seq(2, 15, 1), algorithm="exact",
        eigen.cnt=32L, missingRate=0.025), error_message)
})


test_that("computePoolSyntheticAncestry() must return error when fieldPopIn1KG is numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    dataRefDemo <- data.frame(sample.id=c("SampleA", "SampleB", "SampleC", "SampleD"),
                    pop.group=c("TSI", "TSI", "YRI", "YRI"),
                    superPop=c("EUR", "EUR", "AFR", "AFR"))

    error_message <- "The \'fieldPopIn1KG\' parameter must be a character string."

    expect_error(computePoolSyntheticAncestry(gdsReference=gdsF, gdsSample=gdsF,
            sample.ana.id="test", dataRef=dataRefDemo, spRef="TODO",
            studyIDSyn="synthetic",
            np=1L, listCatPop="EUR", fieldPopIn1KG=33,
            fieldPopInfAnc="SuperPop", kList=seq(2, 15, 1),
            pcaList=seq(2, 15, 1), algorithm="exact",
            eigen.cnt=32L, missingRate=0.025), error_message)
})


#############################################################################
### Tests addPhase1KG2SampleGDSFromFile() results
#############################################################################

context("addPhase1KG2SampleGDSFromFile() results")


test_that("addPhase1KG2SampleGDSFromFile() must return error when gdsReference is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'gdsReference\' must be an object of class \'gds.class\'"

    expect_error(addPhase1KG2SampleGDSFromFile(gdsReference=fileGDS,
            pathProfileGDS=test_path("fixtures"), pathGenotest_path("fixtures"),
            fileSNPsRDS="test", verbose="CANADA"), error_message, fixed=TRUE)
})


test_that("addPhase1KG2SampleGDSFromFile() must return error when verbose is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- paste0("The \'verbose\' parameter must be a ",
                                "logical (TRUE or FALSE).")

    expect_error(addPhase1KG2SampleGDSFromFile(gdsReference=gdsF,
            pathProfileGDS=test_path("fixtures"), pathGenotest_path("fixtures"),
            fileSNPsRDS="test", verbose="CANADA"), error_message, fixed=TRUE)
})


#############################################################################
### Tests computePrunedPCARef() results
#############################################################################

context("computePrunedPCARef() results")


test_that("computePrunedPCARef() must return error when gdsProfile is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    error_message <- "The \'gdsProfile\' must be an object of class \'gds.class\'"

    expect_error(computePrunedPCARef(gdsProfile=fileGDS,
        listRef=c("sample1", "sample2"), np=1L, verbose=FALSE),
        error_message, fixed=TRUE)
})


test_that("computePrunedPCARef() must return error when np is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'np\' parameter must be a single positive integer."

    expect_error(computePrunedPCARef(gdsProfile=gdsF,
        listRef=c("sample1", "sample2"), np="1", verbose=FALSE),
        error_message, fixed=TRUE)
})


test_that("computePrunedPCARef() must return error when verbose is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")
    gdsF <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gdsF), envir=parent.frame())

    error_message <- "The \'verbose\' parameter must be logical (TRUE or FALSE)."

    expect_error(computePrunedPCARef(gdsProfile=gdsF,
    listRef=c("sample1", "sample2"), np=1L, verbose="GLUTEN"),
    error_message, fixed=TRUE)
})


#############################################################################
### Tests runExomeAncestry() results
#############################################################################


context("runExomeAncestry() results")


test_that("runExomeAncestry() must return error when pathOut is numeric", {

    pathOut <- test_path("fixtures")
    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")
    gdsFileAnnot <- test_path("fixtures", "ex1_good_small_1KG_Annot_GDS.gds")

    chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                    study.platform="PLATFORM", stringsAsFactors=FALSE)

    ## Pedigree Study data frame
    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                Case.ID=c("TCGA-H01", "TCGA-H02"),
                Sample.Type=c("DNA", "DNA"),
                Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))

    ## Profiles used for synthetic data set
    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
        "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
        superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)

    error_message <- paste0("The \'pathOut\' must be a character string ",
            "representing the path where the output files will be generated.")

    expect_error(runExomeAncestry(pedStudy=ped, studyDF=studyDF,
        pathProfileGDS=pathOut,
        pathGeno=pathOut, pathOut=33, fileReferenceGDS=fileGDS,
        fileReferenceAnnotGDS=gdsFileAnnot, chrInfo=chrInfo,
        syntheticRefDF=syntheticRefDF,
        genoSource="snp-pileup"), error_message)
})


test_that("runExomeAncestry() must return error when fileReferenceGDS is numeric", {

    pathOut <- test_path("fixtures")
    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")
    gdsFileAnnot <- test_path("fixtures", "ex1_good_small_1KG_Annot_GDS.gds")

    chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                        study.platform="PLATFORM", stringsAsFactors=FALSE)

    ## Pedigree Study data frame
    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                Case.ID=c("TCGA-H01", "TCGA-H02"),
                Sample.Type=c("DNA", "DNA"),
                Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))

    ## Profiles used for synthetic data set
    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
            "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
            superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)

    error_message <- paste0("The \'fileReferenceGDS\' must be a character ",
        "string representing the Reference GDS file. The file must exist.")

    expect_error(runExomeAncestry(pedStudy=ped, studyDF=studyDF,
        pathProfileGDS=pathOut, pathGeno=pathOut, pathOut=pathOut,
        fileReferenceGDS=33, fileReferenceAnnotGDS=gdsFileAnnot,
        chrInfo=chrInfo, syntheticRefDF,
        genoSource="snp-pileup"), error_message)
})


test_that("runExomeAncestry() must return error when fileReferenceAnnotGDS is numeric", {

    pathOut <- test_path("fixtures")
    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")
    gdsFileAnnot <- test_path("fixtures", "ex1_good_small_1KG_Annot_GDS.gds")

    chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                          study.platform="PLATFORM", stringsAsFactors=FALSE)

    ## Pedigree Study data frame
    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                      Case.ID=c("TCGA-H01", "TCGA-H02"),
                      Sample.Type=c("DNA", "DNA"),
                      Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))

    ## Profiles used for synthetic data set
    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
            "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
            superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)

    error_message <- paste0("The \'fileReferenceAnnotGDS\' must be a character",
            " string representing the Reference Annotation GDS file. ",
            "The file must exist.")

    expect_error(runExomeAncestry(pedStudy=ped, studyDF=studyDF,
            pathProfileGDS=pathOut, pathGeno=pathOut, pathOut=pathOut,
            fileReferenceGDS=fileGDS, fileReferenceAnnotGDS=32,
            chrInfo=chrInfo, syntheticRefDF=syntheticRefDF,
            genoSource="snp-pileup"), error_message)
})


test_that("runExomeAncestry() must return error when chrInfo is vector of characters", {

    pathOut <- test_path("fixtures")
    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")
    gdsFileAnnot <- test_path("fixtures", "ex1_good_small_1KG_Annot_GDS.gds")

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                          study.platform="PLATFORM", stringsAsFactors=FALSE)

    ## Pedigree Study data frame
    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                      Case.ID=c("TCGA-H01", "TCGA-H02"),
                      Sample.Type=c("DNA", "DNA"),
                      Diagnosis=c("Cancer", "Cancer"),
                      Source=c("TCGA", "TCGA"), stringsAsFactors=FALSE)

    ## Profiles used for synthetic data set
    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
            "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
            superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)

    error_message <- paste0("The 'chrInfo' parameter must be a ",
            "vector of positive integers.")

    expect_error(runExomeAncestry(pedStudy=ped, studyDF=studyDF,
        pathProfileGDS=pathOut, pathGeno=pathOut, pathOut=pathOut,
        fileReferenceGDS=fileGDS, fileReferenceAnnotGDS=gdsFileAnnot,
        chrInfo=c("ALLO", "TEST"), syntheticRefDF=syntheticRefDF,
        genoSource="snp-pileup"), error_message)
})


test_that("runExomeAncestry() must return error when syntheticRefDF missing column", {

    pathOut <- test_path("fixtures")
    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")
    gdsFileAnnot <- test_path("fixtures", "ex1_good_small_1KG_Annot_GDS.gds")

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                          study.platform="PLATFORM", stringsAsFactors=FALSE)

    ## Pedigree Study data frame
    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                      Case.ID=c("TCGA-H01", "TCGA-H02"),
                      Sample.Type=c("DNA", "DNA"),
                      Diagnosis=c("Cancer", "Cancer"),
                      Source=c("TCGA", "TCGA"), stringsAsFactors=FALSE)

    ## Profiles used for synthetic data set
    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
            "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
            stringsAsFactors=FALSE)

    error_message <- paste0("The reference profile data frame ",
        "\'syntheticRefDF\' is incomplete. One or more mandatory columns are ",
        "missing. The mandatory columns are: \'sample.id\', ",
        "\'pop.group\', \'superPop\'.")

    expect_error(runExomeAncestry(pedStudy=ped, studyDF=studyDF,
        pathProfileGDS=pathOut, pathGeno=pathOut, pathOut=pathOut,
        fileReferenceGDS=fileGDS, fileReferenceAnnotGDS=gdsFileAnnot,
        chrInfo=c(100L, 200L), syntheticRefDF=syntheticRefDF,
        genoSource="snp-pileup"), error_message)
})


test_that("runExomeAncestry() must return error when pathGeno does not exist", {

    pathOut <- test_path("fixtures")
    fileGDS <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")
    gdsFileAnnot <- test_path("fixtures", "ex1_good_small_1KG_Annot_GDS.gds")

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                    study.platform="PLATFORM", stringsAsFactors=FALSE)

    ## Pedigree Study data frame
    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                        Case.ID=c("TCGA-H01", "TCGA-H02"),
                        Sample.Type=c("DNA", "DNA"),
                        Diagnosis=c("Cancer", "Cancer"),
                        Source=c("TCGA", "TCGA"), stringsAsFactors=FALSE)

    ## Profiles used for synthetic data set
    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
                "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
                superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)

    error_message <- "The \'pathGeno\' must be an existing directory."

    expect_error(runExomeAncestry(pedStudy=ped, studyDF=studyDF,
        pathProfileGDS=pathOut, pathGeno="DONOTEXISTDIR", pathOut=pathOut,
        fileReferenceGDS=fileGDS, fileReferenceAnnotGDS=gdsFileAnnot,
        chrInfo=c(100L, 200L), syntheticRefDF=syntheticRefDF,
        genoSource="snp-pileup"), error_message)
})



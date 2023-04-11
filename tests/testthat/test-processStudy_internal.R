### Unit tests for processStudy_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests validatePruningSample() results
#############################################################################

context("validatePruningSample() results")


test_that("validatePruningSample() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validatePruningSample(gdsReference=gdsF, method="corr",
        currentProfile="test.gds", studyID="TCGA",
        listSNP=c("sr10103", "sr10202"), slideWindowMaxBP=1000L,
        thresholdLD=0.008, np=1L, verbose=TRUE, chr=1,
        superPopMinAF=0.002, keepPrunedGDS=TRUE, pathProfileGDS=dataDir,
        keepFile=FALSE, pathPrunedGDS=dataDir, outPrefix="test")

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateComputePoolSyntheticAncestryGr() results
#############################################################################

context("validateComputePoolSyntheticAncestryGr() results")


test_that("validateComputePoolSyntheticAncestryGr() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validateComputePoolSyntheticAncestryGr(
        gdsSample=gdsF, sampleRM="TGCA_01", spRef="TCGA",
        studyIDSyn="TCGA", np=1L, listCatPop=c("AFR", "EAS", "SAS"),
        fieldPopIn1KG="SuperPop",  fieldPopInfAnc="Pop", kList=seq_len(3),
        pcaList=seq_len(10), algorithm="exact", eigen.cnt=12L,
        missing.rate=0.02)

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateCreateStudy2GDS1KG() results
#############################################################################

context("validateCreateStudy2GDS1KG() results")


test_that("validateCreateStudy2GDS1KG() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
            Case.ID=c("TCGA-H01", "TCGA-H02"), Sample.Type=c("DNA", "DNA"),
            Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study",
                        study.platform="WES", stringsAsFactors=FALSE)

    result1 <- RAIDS:::validateCreateStudy2GDS1KG(pathGeno=dataDir,
            pedStudy=ped, fileNameGDS=fileGDS, batch=12, studyDF=studyInfo,
            listProfiles=c("TCGA-H01", "TCGA-H02"), pathProfileGDS=dataDir,
            genoSource="snp-pileup", verbose=TRUE)

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateComputeAncestryFromSyntheticFile() results
#############################################################################

context("validateComputeAncestryFromSyntheticFile() results")


test_that("validateComputeAncestryFromSyntheticFile() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validateComputeAncestryFromSyntheticFile(gdsReference=gdsF,
                    gdsSample=gdsF, listFiles, sample.ana.id,
                    spRef, studyIDSyn="Synthetic", np=1L,
                    listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                    fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
                    kList=c(3,4,5), pcaList=c(4,5,6), algorithm="exact",
                    eigen.cnt=32L, missing.rate=0.02)

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateComputePCARefSample() results
#############################################################################

context("validateComputePCARefSample() results")


test_that("validateComputePCARefSample() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validateComputePCARefSample(gdsSample=gdsF,
                name.id="HCC01", studyIDRef="1KG", np=1L, algorithm="exact",
                eigen.cnt=32L, missingRate=0.02)

    expect_identical(result1, 0L)
})



#############################################################################
### Tests validateAppendStudy2GDS1KG() results
#############################################################################

context("validateAppendStudy2GDS1KG() results")


test_that("validateAppendStudy2GDS1KG() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")
    rdsFile <- file.path(dataDir, "mapSNVSelected_Demo.rds")

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study",  study.platform="WES",
                        stringsAsFactors=FALSE)

    result1 <- RAIDS:::validateAppendStudy2GDS1KG(pathGeno=test_path("fixtures"),
            fileNamePED=rdsFile, fileNameGDS=fileGDS,
            batch=1L, studyDF=studyInfo, listSamples=c("HC01", "HC02"),
            pathProfileGDS=test_path("fixtures"), genoSource="snp-pileup",
            verbose=TRUE)

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateAdd1KG2SampleGDS() results
#############################################################################

context("validateAdd1KG2SampleGDS() results")


test_that("validateAdd1KG2SampleGDS() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    fileGDS <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study",  study.platform="WES",
                        stringsAsFactors=FALSE)

    result1 <- RAIDS:::validateAdd1KG2SampleGDS(gdsReference=gdsF,
                gdsProfileFile=fileGDS, currentProfile="Sample01",
                studyID="Synthetic")

    expect_identical(result1, 0L)
})



#############################################################################
### Tests validateRunExomeAncestry() results
#############################################################################

context("validateRunExomeAncestry() results")


test_that("validateRunExomeAncestry() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    gdsRefFile <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
    gdsRefAnnotFile <- file.path(dataDir, "ex1_good_small_1KG_Annot_GDS.gds")

    ## Pedigree Study data frame
    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                 Case.ID=c("TCGA-H01", "TCGA-H02"),
                 Sample.Type=c("DNA", "DNA"),
                 Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                    study.desc="Pancreatic study",  study.platform="WES",
                    stringsAsFactors=FALSE)

    dataRefSyn <- data.frame(sample.id=c("HG00150", "HG00138",
                    "HG00330", "HG00275"), pop.group=c("GBR", "GBR",
                    "FIN", "FIN"), superPop=c("EUR", "EUR", "EUR", "EUR"),
                    stringsAsFactors=FALSE)

    chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L)

    result <- RAIDS:::validateRunExomeAncestry(pedStudy=ped, studyDF=studyInfo,
        pathProfileGDS=dataDir, pathGeno=dataDir, pathOut=dataDir,
        fileReferenceGDS=gdsRefFile, fileReferenceAnnotGDS=gdsRefAnnotFile,
        chrInfo=chrInfo, dataRefSyn=dataRefSyn, genoSource="snp-pileup")

    expect_identical(result, 0L)
})


#############################################################################
### Tests validatePEDStudyParameter() results
#############################################################################

context("validatePEDStudyParameter() results")


test_that("validatePEDStudyParameter() must return expected results when all input are valid", {

    ## Pedigree Study data frame
    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
                Case.ID=c("TCGA-H01", "TCGA-H02"),
                Sample.Type=c("DNA", "DNA"),
                Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))

    result <- RAIDS:::validatePEDStudyParameter(pedStudy=ped)

    expect_identical(result, 0L)
})


#############################################################################
### Tests validateDataRefSynParameter() results
#############################################################################

context("validateDataRefSynParameter() results")


test_that("validateDataRefSynParameter() must return expected results when all input are valid", {

    ## Profiles used for synthetic data set
    dataRefSyn <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
        "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
        superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)

    result <- RAIDS:::validateDataRefSynParameter(dataRefSyn=dataRefSyn)

    expect_identical(result, 0L)
})


#############################################################################
### Tests validateAddStudy1Kg() results
#############################################################################

context("validateAddStudy1Kg() results")


test_that("validateAddStudy1Kg() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    gdsRefFile <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
    gdsRef <- openfn.gds(gdsRefFile)
    withr::defer((gdsfmt::closefn.gds(gdsRef)), envir = parent.frame())

    fileProfileGDS <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")

    result <- RAIDS:::validateAddStudy1Kg(gdsReference=gdsRef,
                        fileProfileGDS=fileProfileGDS, verbose=FALSE)

    expect_identical(result, 0L)
})

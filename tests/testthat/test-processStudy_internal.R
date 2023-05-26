### Unit tests for processStudy_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)

local_GDS_file <- function(path) {
    GDS_file_tmp  <- createfn.gds(filename=path)
    defer_parent(remove_local_GDS_file(path=path))

    add.gdsn(GDS_file_tmp, "Ref.count", rep(10L, 12))
    add.gdsn(GDS_file_tmp, "Alt.count", rep(12L, 12))
    add.gdsn(GDS_file_tmp, "Total.count", rep(22L, 12))
    add.gdsn(GDS_file_tmp, "lap", rep(0.5, 12))
    sync.gds(GDS_file_tmp)

    return(GDS_file_tmp)
}


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

    spRef <- c("EUR", "AFR")
    names(spRef) <- c("HG01", "HG02")

    result1 <- RAIDS:::validateComputePoolSyntheticAncestryGr(
        gdsProfile=gdsF, sampleRM="TGCA_01", spRef=spRef,
        studyIDSyn="TCGA", np=1L, listCatPop=c("AFR", "EAS", "SAS"),
        fieldPopInfAnc="Pop", kList=seq_len(3),
        pcaList=seq_len(10), algorithm="exact", eigenCount=12L,
        missingRate=0.02, verbose=FALSE)

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
                    gdsProfile=gdsF, listFiles, currentProfile="test",
                    spRef=c("EUR", "AFR"), studyIDSyn="Synthetic", np=1L,
                    listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                    fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
                    kList=c(3,4,5), pcaList=c(4,5,6), algorithm="exact",
                    eigenCount=32L, missingRate=0.02, verbose=FALSE)

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

    result1 <- RAIDS:::validateComputePCARefSample(gdsProfile=gdsF,
        currentProfile="HCC01", studyIDRef="1KG", np=1L, algorithm="exact",
        eigenCount=32L, missingRate=0.02, verbose=FALSE)

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
            filePedRDS=rdsFile, fileNameGDS=fileGDS,
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

    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138",
                    "HG00330", "HG00275"), pop.group=c("GBR", "GBR",
                    "FIN", "FIN"), superPop=c("EUR", "EUR", "EUR", "EUR"),
                    stringsAsFactors=FALSE)

    chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L)

    result <- RAIDS:::validateRunExomeAncestry(pedStudy=ped, studyDF=studyInfo,
        pathProfileGDS=dataDir, pathGeno=dataDir, pathOut=dataDir,
        fileReferenceGDS=gdsRefFile, fileReferenceAnnotGDS=gdsRefAnnotFile,
        chrInfo=chrInfo, syntheticRefDF=syntheticRefDF, genoSource="snp-pileup")

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
    syntheticRefDF <- data.frame(sample.id=c("HG00150", "HG00138", "HG00330",
        "HG00275"), pop.group=c("GBR", "GBR","FIN", "FIN"),
        superPop=c("EUR", "EUR", "EUR", "EUR"), stringsAsFactors=FALSE)

    result <- RAIDS:::validateDataRefSynParameter(syntheticRefDF=syntheticRefDF)

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


#############################################################################
### Tests validateEstimateAllelicFraction() results
#############################################################################

context("validateEstimateAllelicFraction() results")


test_that("validateEstimateAllelicFraction() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")
    gdsRefFile <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")
    gdsRef <- openfn.gds(gdsRefFile)
    withr::defer((gdsfmt::closefn.gds(gdsRef)), envir = parent.frame())

    fileProfileGDS <- file.path(dataDir, "GDS_Sample_with_study_demo.gds")
    gdsSample <- openfn.gds(fileProfileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsSample)), envir = parent.frame())

    chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L)

    result <- RAIDS:::validateEstimateAllelicFraction(gdsReference=gdsRef,
        gdsProfile=gdsSample, currentProfile="ex1", studyID="MYDATA",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=9, cutOffAR=3,
        gdsRefAnnot=NULL, blockID=NULL, verbose=FALSE)

    expect_identical(result, 0L)
})


#############################################################################
### Tests validateProfileGDSExist() results
#############################################################################

context("validateProfileGDSExist() results")


test_that("validateProfileGDSExist() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures/sampleGDSforPruning")

    result <- RAIDS:::validateProfileGDSExist(pathProfile=dataDir,
                                        profile="ex1_demoForPruning")

    expect_identical(result, file.path(dataDir, "ex1_demoForPruning.gds"))
})


#############################################################################
### Tests validateComputePCAMultiSynthetic() results
#############################################################################

context("validateComputePCAMultiSynthetic() results")


test_that("validateComputePCAMultiSynthetic() must return expected results when all input are valid", {

    dataDir <- test_path("fixtures")

    pathFile <- test_path("fixtures/sampleGDSforPoolSyntheticAncestry")
    fileGDS <- test_path(pathFile, "ex1.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsSample)), envir = parent.frame())

    pca <- readRDS(test_path(pathFile, "pca1KG.RDS"))

    result <- RAIDS:::validateComputePCAMultiSynthetic(gdsProfile=gdsSample,
        listPCA=pca, sampleRef=c("HG01"), studyIDSyn="MyStudy",
        verbose=FALSE)

    expect_identical(result, 0L)
})


#############################################################################
### Tests validateComputeKNNRefSynthetic() results
#############################################################################

context("validateComputeKNNRefSynthetic() results")


test_that("validateComputeKNNRefSynthetic() must return expected results when all input are valid", {

    pathFile <- test_path("fixtures/sampleGDSforPoolSyntheticAncestry")
    fileGDS <- test_path(pathFile, "ex1.gds")
    gdsSample <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsSample)), envir = parent.frame())

    pca <- readRDS(test_path(pathFile, "pcaSynthetic.RDS"))

    ## The known ancestry for the 1KG reference profiles
    refKnownSuperPop <- readRDS(file.path(pathFile, "knownSuperPop1KG.RDS"))

    result <- RAIDS:::validateComputeKNNRefSynthetic(gdsProfile=gdsSample,
        listEigenvector=pca, listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
        studyIDSyn="MyStudy", spRef=refKnownSuperPop, fieldPopInfAnc="SuperPop",
        kList=c(11, 12), pcaList=c(10, 11, 12, 13))

    expect_identical(result, 0L)
})


#############################################################################
### Tests validateComputeKNNRefSample() results
#############################################################################

context("validateComputeKNNRefSample() results")


test_that("validateComputeKNNRefSample() must return expected results when all input are valid", {

    pathFile <- test_path("fixtures/sampleGDSforPoolSyntheticAncestry")

    pca <- readRDS(test_path(pathFile, "pcaSynthetic.RDS"))
    pca$sample.id <- pca$sample.id[1]

    ## The known ancestry for the 1KG reference profiles
    refKnownSuperPop <- readRDS(file.path(pathFile, "knownSuperPop1KG.RDS"))

    result <- RAIDS:::validateComputeKNNRefSample(listEigenvector=pca,
        listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
        spRef=refKnownSuperPop, fieldPopInfAnc="SuperPop",
        kList=c(11, 12), pcaList=c(10, 11, 12, 13))

    expect_identical(result, 0L)
})


#############################################################################
### Tests selParaPCAUpQuartile() results
#############################################################################

context("selParaPCAUpQuartile() results")


test_that("selParaPCAUpQuartile() must return expected results when all input are valid", {

    pathFile <- test_path("fixtures/sampleGDSforAncestryByFile")

    pca <- readRDS(file.path(pathFile, "listPCASample.RDS"))

    matKNN <- readRDS(file.path(pathFile, "matKNN.RDS"))

    pedSynt <- readRDS(file.path(pathFile, "pedSyn.RDS"))

    result <- RAIDS:::selParaPCAUpQuartile(matKNN=matKNN,
        pedCall=pedSynt, refCall="superPop", predCall="SuperPop",
        listCall=c("EAS", "EUR", "AFR", "AMR", "SAS"), kList=seq(3,5,1),
        pcaList=seq(10,15,1))

    sampleDFPop <- data.frame(D=c(11, 12, 13, 14, 15),
                K=c(4, 5, 5, 3, 5),
                AUROC.min=c(0.539772727272727, 0.528409090909091,
                    0.517045454545455, 0.528409090909091, 0.602272727272727),
                AUROC=c(0.685267857142857, 0.710267857142857,
                    0.722767857142857, 0.726339285714286, 0.738392857142857),
                Accu.CM=c(0.538461538461538, 0.576923076923077,
                    0.596153846153846, 0.596153846153846, 0.615384615384615),
                stringsAsFactors=FALSE)
    rownames(sampleDFPop) <- c(5L, 9L, 12L, 13L, 18L)

    sampleDFPCA <- data.frame(D=c(10, 11, 12, 13, 14, 15),
        median=c(0.539772727272727, 0.539772727272727, 0.528409090909091,
                0.528409090909091, 0.539772727272727, 0.602272727272727),
        mad=c(0.016847727272727, 0.016847727272727, 0.075814772727273,
                0.016847727272727, 0.016847727272727, 0.000000000000000),
        upQuartile=c(0.545454545454545, 0.559659090909091, 0.559659090909091,
                0.559659090909091, 0.559659090909091, 0.602272727272727),
        K=c(4, 5, 4, 4, 4, 4), stringsAsFactors=FALSE)
    rownames(sampleDFPCA) <- c("75%",  "75%1", "75%2", "75%3", "75%4", "75%5")

    sampleDFAUROC <- data.frame(pcaD=c(10, 11, 11, 12, 13),
        K=c(3, 3, 4, 3, 5),
        Call=c("SAS", "AMR", "EUR", "AFR", "AFR"),
        L=c(0.520945878038968, 0.400250400882566, 0.472946268504436,
                    0.718231321124318, 0.784535852000281),
        AUC=c(0.692857142857143, 0.528409090909091, 0.640476190476190,
                    0.840225563909774, 0.889097744360902),
        H=c(0.864768407675317, 0.656567780935616, 0.808006112447945,
                    0.962219806695231, 0.993659636721523),
        stringsAsFactors=FALSE)
    rownames(sampleDFAUROC) <- c(5L, 19L, 22L, 33L, 58L)

    expect_true(is.list(result))

    expect_true(is.data.frame(result$dfPCA))
    expect_identical(colnames(result$dfPCA),
                        c("D", "median", "mad", "upQuartile", "K"))
    expect_identical(nrow(result$dfPCA), 6L)
    expect_equal(result$dfPCA, sampleDFPCA)

    expect_true(is.data.frame(result$dfAUROC))
    expect_identical(colnames(result$dfAUROC),
                     c("pcaD", "K", "Call", "L", "AUC", "H"))
    expect_identical(nrow(result$dfAUROC), 90L)
    expect_equal(result$dfAUROC[c(5, 19, 22, 33, 58),], sampleDFAUROC)

    expect_true(is.data.frame(result$dfPop))
    expect_identical(colnames(result$dfPop),
                        c("D", "K", "AUROC.min", "AUROC", "Accu.CM"))
    expect_identical(nrow(result$dfPop), 18L)
    expect_equal(result$dfPop[c(5, 9, 12, 13, 18),], sampleDFPop)

    expect_identical(result$D, 15)
    expect_identical(result$K, 4)
    expect_identical(result$listD, 15)
})


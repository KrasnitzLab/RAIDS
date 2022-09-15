### Unit tests for processStudy_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests validatePruningSample() results
#############################################################################

context("validatePruningSample() results")


test_that("validatePruningSample() must return epxected results when all input are valid", {

    data.dir <- test_path("fixtures")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validatePruningSample(gds=gdsF, method="corr",
        sampleCurrent="test.gds", study.id="TCGA",
        listSNP=c("sr10103", "sr10202"), slide.max.bp.v=1000L,
        ld.threshold.v=0.008, np=1L, verbose.v=TRUE, chr=1,
        minAF.SuperPop=0.002, keepGDSpruned=TRUE, PATHSAMPLEGDS=data.dir,
        keepFile=FALSE, PATHPRUNED=data.dir, outPref="test")

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateComputePoolSyntheticAncestryGr() results
#############################################################################

context("validateComputePoolSyntheticAncestryGr() results")


test_that("validateComputePoolSyntheticAncestryGr() must return epxected results when all input are valid", {

    data.dir <- test_path("fixtures")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validateComputePoolSyntheticAncestryGr(gds=gdsF,
        gdsSample=gdsF, sampleRM="TGCA_01", spRef="TCGA",
        study.id.syn="TCGA", np=1L, listCatPop=c("AFR", "EAS", "SAS"),
        fieldPopIn1KG="SuperPop",  fieldPopInfAnc="Pop", kList=seq_len(3),
        pcaList=seq_len(10), algorithm="exact", eigen.cnt=12L,
        missing.rate=0.02)

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateCreateStudy2GDS1KG() results
#############################################################################

context("validateCreateStudy2GDS1KG() results")


test_that("validateCreateStudy2GDS1KG() must return epxected results when all input are valid", {

    data.dir <- test_path("fixtures")
    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    ped <- data.frame(Name.ID=c("Sample_01", "Sample_02"),
            Case.ID=c("TCGA-H01", "TCGA-H02"), Sample.Type=c("DNA", "DNA"),
            Diagnosis=c("Cancer", "Cancer"), Source=c("TCGA", "TCGA"))

    studyInfo <- data.frame(study.id="Pancreatic.WES",
                        study.desc="Pancreatic study",
                        study.platform="WES", stringsAsFactors=FALSE)

    result1 <- RAIDS:::validateCreateStudy2GDS1KG(pedStudy=ped,
            fileNameGDS=gdsFile, batch=12, studyDF=studyInfo,
            listSamples=c("TCGA-H01", "TCGA-H02"), PATHSAMPLEGDS=data.dir,
            verbose=TRUE)

    expect_identical(result1, 0L)
})


#############################################################################
### Tests validateComputeAncestryFromSyntheticFile() results
#############################################################################

context("validateComputeAncestryFromSyntheticFile() results")


test_that("validateComputeAncestryFromSyntheticFile() must return epxected results when all input are valid", {

    data.dir <- test_path("fixtures")
    gdsFile <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFile)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    result1 <- RAIDS:::validateComputeAncestryFromSyntheticFile(gds=gdsF,
                    gdsSample=gdsF, listFiles, sample.ana.id,
                    spRef, study.id.syn="Synthetic", np=1L,
                    listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"),
                    fieldPopIn1KG="superPop", fieldPopInfAnc="SuperPop",
                    kList=c(3,4,5), pcaList=c(4,5,6), algorithm="exact",
                    eigen.cnt=32L, missing.rate=0.02)

    expect_identical(result1, 0L)
})

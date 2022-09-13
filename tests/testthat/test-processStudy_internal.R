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


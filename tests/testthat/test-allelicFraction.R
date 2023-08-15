### Unit tests for allelicFraction.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests computeAlleleFraction() results
#############################################################################

context("computeAlleleFraction() results")


test_that("computeAlleleFraction() must return expected results when not imbalanced SNVs", {

    ## Data frame with SNV information for the specified chromosome (chr 1)
    snpInfo <- data.frame(cnt.tot=c(41, 17, 27, 15, 11, 37, 16, 32),
        cnt.ref=c(40, 17, 27, 15, 4, 14, 16, 32),
        cnt.alt=c(0, 0, 0, 0, 7, 23, 0, 0),
        snp.pos=c(3722256, 3722328, 3767522, 3868160, 3869467, 4712655,
            6085318, 6213145),
        snp.chr=c(rep(1, 8)),
        normal.geno=c(rep(1, 8)), pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE,
        TRUE, TRUE),
        pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, rep(TRUE, 3)),
        snp.index=c(160, 162, 204, 256, 259, 288, 366, 465),
        keep=rep(TRUE, 8), hetero=c(rep(FALSE, 4), TRUE, TRUE, rep(FALSE, 2)),
        homo=c(rep(TRUE, 4), FALSE, FALSE, TRUE, TRUE),
        lap=rep(-1, 8), LOH=rep(0, 8), imbAR=rep(-1, 8),
        stringAsFactor=FALSE)

    result <- RAIDS:::computeAlleleFraction(snpPos=snpInfo, w=10, cutOff=-3)
    expect_equal(result, NULL)
})



test_that("computeAlleleFraction() must return expected results when imbalanced SNVs", {

    set.seed(112)

    ## Data frame with SNV information for the specified chromosome (chr 1)
    snpInfo <- data.frame(cnt.tot=c(41, 17, 27, 15, 11, 37, 16, 32),
        cnt.ref=c(40, 37, 37, 25, 4, 14, 26, 32),
        cnt.alt=c(1, 8, 6, 5, 0, 0, 3, 4),
        snp.pos=c(3722256, 3722328, 3767522, 3868160, 3869467, 4712655,
                                    6085318, 6213145),
        snp.chr=c(rep(1, 8)),
        normal.geno=c(rep(1, 8)),
        pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, rep(TRUE, 3)),
        snp.index=c(160, 162, 204, 256, 259, 288, 366, 465),
        keep=rep(TRUE, 8), hetero=c(rep(FALSE, 4), TRUE, TRUE, rep(FALSE, 2)),
        homo=c(rep(TRUE, 4), FALSE, FALSE, TRUE, TRUE),
        lap=rep(-1, 8), LOH=rep(0, 8), imbAR=c(1, 1, 1, 1, 1, 0, 1, 1),
        stringAsFactor=FALSE)

    result <- RAIDS:::computeAlleleFraction(snpPos=snpInfo, w=10, cutOff=-3)
    expect_equal(result, matrix(c(1, 5, 0, 7, 8, NA), byrow=TRUE,
                                        nrow=2))
})



test_that("computeAlleleFraction() must return expected results with small w parameter", {

    set.seed(112)

    ## Data frame with SNV information for the specified chromosome (chr 1)
    snpInfo <- data.frame(cnt.tot=c(41, 17, 27, 15, 11, 37, 29, 36, 34, 38, 29, 26),
                          cnt.ref=c(40, 37, 37, 25,  4, 14, 26, 32, 26, 32, 34, 31),
                          cnt.alt=c(1,   8,  6,  5, 12, 54,  3,  4,  5,  6,  5,  5),
                          snp.pos=c(3722256, 3722328, 3767522, 3868160, 3869467, 4712655,
                                    6085318, 6213145, 6385318, 6313145, 6485318, 6413145),
                          snp.chr=c(rep(1, 12)),
                          normal.geno=c(rep(1, 12)),
                          pruned=rep(TRUE, 12),
                          snp.index=c(160, 162, 204, 256, 259, 288, 366,
                                      465, 467, 477, 488, 489),
                          keep=rep(TRUE, 12), hetero=c(rep(TRUE, 12)),
                          homo=c(rep(FALSE, 12)),
                          lap=rep(-1, 12), LOH=rep(0, 12), imbAR=rep(1, 12),
                          stringAsFactor=FALSE)

    result <- RAIDS:::computeAlleleFraction(snpPos=snpInfo, w=4, cutOff=-3)
    expect_equal(result, matrix(c(1, 12, 0.148714810281518), byrow=TRUE,
                                nrow=1))
})


#############################################################################
### Tests estimateAllelicFraction() results
#############################################################################


context("estimateAllelicFraction() results")


test_that("estimateAllelicFraction() must return error when gdsReference is a character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'gdsReference\' must be an object of class \'gds.class\'."

    expect_error(estimateAllelicFraction(gdsReference="test.gds",
        gdsProfile=gdsF, currentProfile="gds", studyID= "TCGA",
        chrInfo=chrInfo, studyType="DNA", minCov=10, minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when gdsSample is a character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'gdsProfile\' must be an object of class \'gds.class\'."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile="test.gds", currentProfile="gds", studyID= "TCGA",
        chrInfo=chrInfo, studyType="DNA", minCov=10, minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when currentProfile is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'currentProfile\' must be a character string."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile=33, studyID= "TCGA",
        chrInfo=chrInfo, studyType="DNA", minCov=10, minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when studyID is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'studyID\' must be a character string."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID= 33,
        chrInfo=chrInfo, studyType="DNA", minCov=10, minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when studyType is a numeric value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'studyType\' must be a character string."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType=33, minCov=10, minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minCov is vector of numeric values", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'minCov\' must be a single numeric positive value"

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=c(10, 11), minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minCov is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'minCov\' must be a single numeric positive value"

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov="3", minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minProb is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb="0.1",
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minProb is negative value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=-0.1,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minProb is above 1", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=1.01,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when minProb is vector of numeric values", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=c(0.01, 0.02),
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb="0.001", cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is negative value", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb=-0.01, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is above 1", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb=1.01, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when eProb is vector of numeric values", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb=c(0.01, 0.001), cutOffLOH=-5, cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when wAR is vector of numeric values", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'wAR\' parameter must be a single numeric positive value."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb=0.01, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR="9"), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when cutOffLOH is vector of numeric values", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'cutOffLOH\' parameter must be a single numeric value."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb=0.01, cutOffLOH=c(-5, -3), cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when cutOffLOH is a character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'cutOffLOH\' parameter must be a single numeric value."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb=0.01, cutOffLOH="-3", cutOffHomoScore=-3,
        wAR=9), error_message, fixed=TRUE)
})




test_that("estimateAllelicFraction() must return error when cutOffLOH is vector of numeric values", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'cutOffHomoScore\' parameter must be a single numeric value."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb=0.01, cutOffLOH=-5, cutOffHomoScore=c(-3, -4),
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when cutOffHomoScore is a character string", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'cutOffHomoScore\' parameter must be a single numeric value."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.1,
        eProb=0.01, cutOffLOH=-5, cutOffHomoScore="-3",
        wAR=9), error_message, fixed=TRUE)
})


test_that("estimateAllelicFraction() must return error when verbose is not a numeric", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    error_message <- "The \'verbose\' parameter must be a logical (TRUE or FALSE)."

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA",
        chrInfo=chrInfo, studyType="TOTO", minCov=10L, minProb=0.1,
        eProb=0.01, cutOffLOH=-5, cutOffHomoScore=-3, wAR=9, verbose=22),
        error_message, fixed=TRUE)
})

test_that("estimateAllelicFraction() must return error when studyType is not a valid choice", {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    expect_error(estimateAllelicFraction(gdsReference=gdsF,
        gdsProfile=gdsF, currentProfile="test", studyID="DNA_1",
        chrInfo=chrInfo, studyType="TOTO", minCov=10L, minProb=0.1,
        eProb=0.01, cutOffLOH=-5, cutOffHomoScore=-3, wAR=9))
})


test_that("estimateAllelicFraction() must return valid results", {

    file1KG <- test_path("fixtures", "ex1_good_small_1KG_GDS.gds")

    gds1KG <- openfn.gds(file1KG)
    withr::defer((gdsfmt::closefn.gds(gds1KG)), envir=parent.frame())


    dataDirSample <- test_path("fixtures/sampleGDSforEstimAlleFraction")
    file.copy(file.path(dataDirSample, "ex1_demoForEstimAllFrac.gds"),
                file.path(dataDirSample, "ex1.gds"))
    withr::defer((unlink(file.path(dataDirSample, "ex1.gds"))),
                    envir=parent.frame())

    gdsProfile <- openfn.gds(file.path(dataDirSample, "ex1.gds"),
                                readonly=FALSE)
    withr::defer((gdsfmt::closefn.gds(gdsProfile)), envir=parent.frame())

    chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
        181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
        135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
        83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
        156040895L, 57227415L,  16569L)

    result <- estimateAllelicFraction(gdsReference=gds1KG,
        gdsProfile=gdsProfile, currentProfile="ex1", studyID="MYDATA",
        chrInfo=chrInfo, studyType="DNA", minCov=10L, minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=9, cutOffAR=3,
        gdsRefAnnot=NULL, blockID=NULL, verbose=FALSE)

    expect_equal(result, 0L)

    result1 <- read.gdsn(index.gdsn(node=gdsProfile, path="lap"))
    result2 <- read.gdsn(index.gdsn(node=gdsProfile, path="segment"))

    expect_equal(result1, rep(0.5, 61))
    expect_equal(result2, c(rep(1, 6), 2, 3, rep(4, 5), rep(5, 4), rep(6, 3),
        rep(7, 2), 9, rep(10, 5), rep(11, 4), 12, 12, rep(13, 4), rep(14, 4),
        15, rep(16, 3), rep(17, 3), 18, 18, 19, 19, rep(20, 3), 21, rep(22, 4)))
})


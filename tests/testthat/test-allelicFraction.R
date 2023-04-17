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

    result <- RAIDS:::computeAlleleFraction(snp.pos=snpInfo, w=10, cutOff=-3)
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

    result <- RAIDS:::computeAlleleFraction(snp.pos=snpInfo, w=10, cutOff=-3)
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

    result <- RAIDS:::computeAlleleFraction(snp.pos=snpInfo, w=4, cutOff=-3)
    expect_equal(result, matrix(c(1, 12, 0.148714810281518), byrow=TRUE,
                                nrow=1))
})

#############################################################################
### Tests computeAlleleFraction() results
#############################################################################

context("computeLOHBlocksDNAChr() results")


test_that("computeLOHBlocksDNAChr() must return expected results", {

    dataDir <- testthat::test_path("fixtures")
    fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")

    ## Open the reference GDS file
    gds1KG <- snpgdsOpen(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gds1KG)), envir = parent.frame())

    chrInfo <- c(248956422L, 242193529L, 198295559L, 190214555L,
        181538259L, 170805979L, 159345973L, 145138636L, 138394717L, 133797422L,
        135086622L, 133275309L, 114364328L, 107043718L, 101991189L, 90338345L,
        83257441L,  80373285L,  58617616L,  64444167L,  46709983L, 50818468L,
        156040895L, 57227415L,  16569L)


    ## Data frame with SNV information for the specified chromosome (chr 1)
    snpInfo <- data.frame(cnt.tot=c(41, 17, 27, 15, 11, 37, 29, 36, 26, 32, 29, 26),
                    cnt.ref=c(40, 37, 37, 25,  4, 14, 26, 32, 26, 32, 34, 31),
                    cnt.alt=c(1,   8,  6,  5, 12, 54,  3,  1,  0,  0,  2,  5),
                    snp.pos=c(3722256, 3722328, 3767522, 3868160, 3869467, 4712655,
                            6085318, 6213145, 6385318, 6313145, 6485318, 6413145),
                    snp.chr=c(rep(1, 12)),
                    normal.geno=c(rep(1, 12)),
                    pruned=rep(TRUE, 12),
                    snp.index=c(160, 162, 204, 256, 259, 288, 366,
                                    465, 467, 477, 488, 489),
                    keep=rep(TRUE, 12),
                    hetero=c(rep(TRUE, 8), FALSE, TRUE, TRUE, TRUE),
                    homo=c(rep(FALSE, 8), TRUE, FALSE, FALSE, FALSE),
                    stringAsFactor=FALSE)

    result <- RAIDS:::computeLOHBlocksDNAChr(gdsReference=gds1KG,
         chrInfo=chrInfo, snp.pos=snpInfo, chr=1L, genoN=0.0001)

    expected <- data.frame(chr=rep(1, 12),
                           start=c(1, snpInfo$snp.pos[seq_len(8)]+1, 6313146,
                                        6485319, 6413146),
                           end=c(snpInfo$snp.pos[seq_len(8)], 6313145,
                                        6485318, 6413145, 248956422),
                           logLHR=c(rep(0, 9), -100, 0, 0),
                           LH1=c(rep(0, 9), -10.729445965503494, 0, 0),
                           LM1=c(rep(0, 9), -0.912269425836972, 0, 0),
                           homoScore=c(rep(0, 9), -9.817176539666523, 0, 0),
                           nbSNV=c(rep(0, 9), 1, 0, 0),
                           nbPruned=c(rep(0, 9), 1, 0, 0),
                           nbNorm=c(rep(0, 9), 1, 0, 0),
                           LOH=c(rep(0, 12)))

    expect_equal(result, expected)

})


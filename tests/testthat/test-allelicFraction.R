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

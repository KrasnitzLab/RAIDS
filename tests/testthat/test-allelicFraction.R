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
        cnt.alt=c(1, 8, 6, 5, 27, 53, 3, 0),
        snp.pos=c(3722256, 3722328, 3767522, 3868160, 3869467, 4712655,
                                    6085318, 6213145),
        snp.chr=c(rep(1, 8)),
        normal.geno=c(rep(1, 8)), pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE,
                                                             TRUE, TRUE),
        pruned=c(TRUE, TRUE, FALSE, TRUE, FALSE, rep(TRUE, 3)),
        snp.index=c(160, 162, 204, 256, 259, 288, 366, 465),
        keep=rep(TRUE, 8), hetero=c(rep(FALSE, 4), TRUE, TRUE, rep(FALSE, 2)),
        homo=c(rep(TRUE, 4), FALSE, FALSE, TRUE, TRUE),
        lap=rep(-1, 8), LOH=rep(0, 8), imbAR=rep(1, 8),
        stringAsFactor=FALSE)

    result <- RAIDS:::computeAlleleFraction(snp.pos=snpInfo, w=10, cutOff=-3)
    expect_equal(result, matrix(c(1, 8, 0.168993740972557), byrow=FALSE,
                                        nrow=1))
})

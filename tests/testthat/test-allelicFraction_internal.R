### Unit tests for allelicFraction_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests testEmptyBox() results
#############################################################################

context("testEmptyBox() results")


test_that("testEmptyBox() must return expected results", {

    ## Data frame with SNV information
    snpInfo <- snpInfo <- data.frame(
        cnt.ref=c(40, 17, 27, 15, 4, 14, 16, 32),
        cnt.alt=c(2, 4, 5, 10, 7, 23, 0, 0))

    result <- RAIDS:::testEmptyBox(matCov=snpInfo, pCutOff=-3)

    expected <- list()
    expected[["pWin"]] <- c(0.000000000411092, 0.007197380065918,
                                0.000113074202091, 0.424356222152710,
                                0.548828125000000, 0.187741558998824,
                                0.000030517578125, 0.000000000465661)
    expected[["p"]] <- -12.262362633720395
    expected[["pCut"]] <- 0
    expected[["pCut1"]] <- 1

    expect_equal(result, expected)
})


#############################################################################
### Tests computeLOHBlocksDNAChr() results
#############################################################################

context("computeLOHBlocksDNAChr() results")


test_that("computeLOHBlocksDNAChr() must return expected results", {

    dataDir <- testthat::test_path("fixtures")
    fileGDS <- file.path(dataDir, "ex1_good_small_1KG_GDS.gds")

    ## Open the reference GDS file
    gds1KG <- snpgdsOpen(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gds1KG)), envir = parent.frame())

    chromosome <- c(248956422L, 242193529L, 198295559L, 190214555L,
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
                    chrInfo=chromosome, snpPos=snpInfo, chr=1L, genoN=0.0001)

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

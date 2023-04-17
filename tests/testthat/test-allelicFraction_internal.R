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

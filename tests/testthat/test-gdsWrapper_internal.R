### Unit tests for gdsWrapper_internal.R functions

library(RAIDS)
library(withr)
library(testthat)
library(SNPRelate)


#############################################################################
### Tests runIBDKING() results
#############################################################################

context("runIBDKING() results")


test_that("runIBDKING() must return expected results", {

    ## Open genotype GDS file
    genoFile <- snpgdsOpen(snpgdsExampleFileName())
    withr::defer((snpgdsClose(genoFile)), envir = parent.frame())

    samples <- read.gdsn(index.gdsn(genoFile, "sample.id"))
    CEU <- samples[
            read.gdsn(index.gdsn(genoFile, "sample.annot/pop.group"))=="CEU"]
    CEU <- CEU[1:3]

    result <- RAIDS:::runIBDKING(gds=genoFile, profileID=CEU, snpID=NULL,
                                        maf=0.05, verbose=FALSE)

    sampleIDs <- c("NA07034", "NA07055", "NA12814")
    ibs0 <- matrix(c(0.000000000000000, 0.079705823891870, 0.088194306191519,
        0.079705823891870, 0.000000000000000, 0.094035785288270,
        0.088194306191519, 0.094035785288270, 0.000000000000000), nrow=3,
        byrow=TRUE)
    kinship <- matrix(c(0.500000000000000, 0.017237640936687, 0.001294777729823,
        0.017237640936687, 0.500000000000000, -0.014341590612777,
        0.001294777729823, -0.014341590612777, 0.500000000000000), nrow=3,
        byrow=TRUE)

    expect_equal(result$sample.id, sampleIDs)
    expect_equal(result$IBS0, ibs0)
    expect_equal(result$kinship, kinship)
})


#############################################################################
### Tests runLDPruning() results
#############################################################################

context("runLDPruning() results")


test_that("runLDPruning() must return expected results", {

    ## Open genotype GDS file
    genoFile <- snpgdsOpen(snpgdsExampleFileName())
    withr::defer((snpgdsClose(genoFile)), envir = parent.frame())

    set.seed(121)

    result <- RAIDS:::runLDPruning(gds=genoFile, verbose=FALSE)

    expect_equal(class(result), "list")
    expect_equal(length(result), 22)
    expect_equal(names(result) , paste0("chr", 1:22))
})

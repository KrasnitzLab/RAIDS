### Unit tests for gdsWrapper functions

library(RAIDS)
library(withr)
library(testthat)


remove_local_GDS_file <- function(path) {
    unlink(x=path, force=TRUE)
}

local_GDS_file <- function(path) {
    GDS_file_tmp  <- createfn.gds(filename=path)
    defer_parent(remove_local_GDS_file(path=path))

    add.gdsn(GDS_file_tmp, "Ref.count", rep(10L, 12))
    add.gdsn(GDS_file_tmp, "Alt.count", rep(12L, 12))
    add.gdsn(GDS_file_tmp, "Total.count", rep(22L, 12))
    sync.gds(GDS_file_tmp)

    return(GDS_file_tmp)
}


#############################################################################
### Tests addUpdateSegment() results
#############################################################################

context("addUpdateSegment() results")


test_that("addUpdateSegment() must copy the expected entry in \"segment\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "GDS_TEMP_01.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

    ## Vector of segment identifiers
    segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)

    ## Add segments to the GDS file
    RAIDS:::addUpdateSegment(gds=GDS_file_tmp, snp.seg=segments)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="segment"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, segments)
})


test_that("addUpdateSegment() must copy the expected entry in \"segment\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFile <- file.path(data.dir, "GDS_TEMP_02.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(gdsFile)

    ## Vector of segment identifiers
    segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L)
    ## Vector of segment identifiers
    segments2 <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L, 5L, 5L)

    ## Add segments to the GDS file
    RAIDS:::addUpdateSegment(gds=GDS_file_tmp, snp.seg=segments)

    ## Update segments to the GDS file
    RAIDS:::addUpdateSegment(gds=GDS_file_tmp, snp.seg=segments2)

    ## Read segments information from GDS file
    results <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="segment"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results, segments2)
})

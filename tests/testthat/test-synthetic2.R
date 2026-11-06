### Unit tests for synthetic2.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests select1KGPopForSynthetic() results
#############################################################################

context("select1KGPopForSynthetic() results")


test_that("select1KGPopForSynthetic() must return error when fileReferenceGDS is not existing", {

    expected <- paste0("The \'fileReferenceGDS\' must be a character string ",
                "representing the Reference GDS file. The file must exist.")

    expect_error(select1KGPopForSynthetic(fileReferenceGDS="not_a_file.gds",
                    nbProfiles=2L), expected, fixed=TRUE)
})


test_that("select1KGPopForSynthetic() must return error when nbProfiles is not positive", {

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")

    expected <- "The \'nbProfiles\' parameter must be a single positive integer."

    expect_error(select1KGPopForSynthetic(fileReferenceGDS=fileRef,
                    nbProfiles=0L), expected, fixed=TRUE)
    expect_error(select1KGPopForSynthetic(fileReferenceGDS=fileRef,
                    nbProfiles=c(1L, 2L)), expected, fixed=TRUE)
})


test_that("select1KGPopForSynthetic() must return nbProfiles profiles per population", {

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")

    gdsRef <- openfn.gds(fileRef)
    sampleId <- read.gdsn(index.gdsn(gdsRef, "sample.id"))
    sampleAnnot <- read.gdsn(index.gdsn(gdsRef, "sample.annot"))
    closefn.gds(gdsRef)

    set.seed(2)
    result <- select1KGPopForSynthetic(fileReferenceGDS=fileRef,
                                            nbProfiles=2L)

    expect_equal(colnames(result), c("sample.id", "pop.group", "superPop"))
    expect_equal(nrow(result), 2 * length(unique(sampleAnnot$pop.group)))
    expect_true(all(table(result$pop.group) == 2))
    expect_false(any(duplicated(result$sample.id)))

    ## The population information matches the reference
    pos <- match(result$sample.id, sampleId)
    expect_equal(result$pop.group, sampleAnnot$pop.group[pos])
    expect_equal(result$superPop, sampleAnnot$superPop[pos])
    expect_true(noOpenGDS())
})


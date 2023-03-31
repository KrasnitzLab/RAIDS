### Unit tests for allelicFraction.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests computeLOHBlocksDNAChr() results
#############################################################################


context("computeLOHBlocksDNAChr() results")


test_that("computeLOHBlocksDNAChr() must return error when chr is vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'chr\' must be a single integer value representing ",
                                "a chromosome")

    expect_error(RAIDS:::computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                            snp.pos=data.frame(), chr=c(1, 2), genoN=0.001),
                 error_message)
})


test_that("computeLOHBlocksDNAChr() must return error when chr is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'chr\' must be a single integer value representing ",
                                "a chromosome")

    expect_error(RAIDS:::computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                        snp.pos=data.frame(), chr="MONTREAL", genoN=0.001),
                 error_message)
})

test_that("computeLOHBlocksDNAChr() must return error when genoN is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'genoN\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(RAIDS:::computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                                    snp.pos=data.frame(), chr=2,
                                    genoN="TOTO"), error_message)
})



test_that("computeLOHBlocksDNAChr() must return error when genoN is vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'genoN\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(RAIDS:::computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                                        snp.pos=data.frame(), chr=2,
                                        genoN=c(0.2, 0.33)), error_message)
})


test_that("computeLOHBlocksDNAChr() must return error when genoN is negative numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'genoN\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(RAIDS:::computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                                        snp.pos=data.frame(), chr=2,
                                        genoN=-0.33), error_message)
})


test_that("computeLOHBlocksDNAChr() must return error when genoN is numeric above 1", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'genoN\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(RAIDS:::computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                                        snp.pos=data.frame(), chr=2,
                                        genoN=1.33), error_message)
})


test_that("computeLOHBlocksDNAChr() must return error when chr not in chrInfo", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'chr\' must be present in the ",
                                        "\'chrInfo\' parameter.")

    expect_error(RAIDS:::computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                                        snp.pos=data.frame(), chr=4,
                                        genoN=0.03), error_message)
})


test_that("computeLOHBlocksDNAChr() must return error when snp.pos is a numeric", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'snp.pos\' must be a data.frame."

    expect_error(RAIDS:::computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                                        snp.pos=33, chr=2,
                                        genoN=0.03), error_message)
})

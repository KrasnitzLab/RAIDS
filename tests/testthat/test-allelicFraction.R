### Unit tests for allelicFraction.R functions

library(RAIDS)
library(withr)
library(gdsfmt)



#############################################################################
### Tests getTableSNV() results
#############################################################################

context("getTableSNV() results")


test_that("getTableSNV() must return error when minCov is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'minCov\' must be a single numeric positive value."

    expect_error(getTableSNV(gds=gdsF, gdsSample=gdsF,
                                sampleCurrent="TEST", study.id="TEST",
                minCov="CANADA", minProb=0.999, eProb=0.001), error_message)
})


test_that("getTableSNV() must return error when minProb is a character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
            "value between 0 and 1.")

    expect_error(getTableSNV(gds=gdsF, gdsSample=gdsF,
                        sampleCurrent="TEST", study.id="TEST",
                        minCov=12, minProb="QUEBEC", eProb=0.001), error_message)
})


test_that("getTableSNV() must return error when minProb is negative", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                            "value between 0 and 1.")

    expect_error(getTableSNV(gds=gdsF, gdsSample=gdsF,
                        sampleCurrent="TEST", study.id="TEST",
                        minCov=12, minProb=-0.32, eProb=0.001), error_message)
})


test_that("getTableSNV() must return error when minProb is superior to 1", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                            "value between 0 and 1.")

    expect_error(RAIDS:::getTableSNV(gds=gdsF, gdsSample=gdsF,
                        sampleCurrent="TEST", study.id="TEST",
                        minCov=12, minProb=1.32, eProb=0.001), error_message)
})


test_that("getTableSNV() must return error when minProb is vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                            "value between 0 and 1.")

    expect_error(RAIDS:::getTableSNV(gds=gdsF, gdsSample=gdsF,
                                sampleCurrent="TEST", study.id="TEST",
                                minCov=12, minProb=c(0.32, 0.44), eProb=0.001),
                    error_message)
})



test_that("getTableSNV() must return error when eprob is negative", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                        "value between 0 and 1.")

    expect_error(RAIDS:::getTableSNV(gds=gdsF, gdsSample=gdsF,
                    sampleCurrent="TEST", study.id="TEST",
                    minCov=12, minProb=0.02, eProb="0.001"), error_message)
})


test_that("getTableSNV() must return error when eProb is superior to 1", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                        "value between 0 and 1.")

    expect_error(RAIDS:::getTableSNV(gds=gdsF, gdsSample=gdsF,
                        sampleCurrent="TEST", study.id="TEST",
                        minCov=12, minProb=0.32, eProb=1.011), error_message)
})


test_that("getTableSNV() must return error when eProb is vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                            "value between 0 and 1.")

    expect_error(RAIDS:::getTableSNV(gds=gdsF, gdsSample=gdsF,
                sampleCurrent="TEST", study.id="TEST",
                minCov=12, minProb=0.32, eProb=c(0.001, 0.2)), error_message)
})


test_that("getTableSNV() must return error when gds is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'gds\' must be an object of class ",
                                "\'gds.class\'.")

    expect_error(RAIDS:::getTableSNV(gds="HELLO", gdsSample=gdsF,
                        sampleCurrent="TEST", study.id="TEST",
                        minCov=12, minProb=0.32, eProb=0.001), error_message)
})


test_that("getTableSNV() must return error when gdsSample is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'gdsSample\' must be an object of class ",
                            "\'gds.class\'.")

    expect_error(RAIDS:::getTableSNV(gds=gdsF, gdsSample="BABY",
                        sampleCurrent="TEST", study.id="TEST",
                        minCov=12, minProb=0.32, eProb=0.001), error_message)
})


test_that("getTableSNV() must return error when verbose is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'verbose\' parameter must be a logical."

    expect_error(RAIDS:::getTableSNV(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="TEST", study.id="TEST",
        minCov=12, minProb=0.32, eProb=0.001, verbose="blue"), error_message)
})


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

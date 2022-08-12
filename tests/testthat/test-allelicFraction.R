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

    expect_error(getTableSNV(gds=gdsF, gdsSample=gdsF,
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

    expect_error(getTableSNV(gds=gdsF, gdsSample=gdsF,
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

    expect_error(getTableSNV(gds=gdsF, gdsSample=gdsF,
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

    expect_error(getTableSNV(gds=gdsF, gdsSample=gdsF,
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

    expect_error(getTableSNV(gds=gdsF, gdsSample=gdsF,
                sampleCurrent="TEST", study.id="TEST",
                minCov=12, minProb=0.32, eProb=c(0.001, 0.2)), error_message)
})


test_that("getTableSNV() must return error when gds is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'gds\' must be an object of class ",
                                "\'gdsn.class\' or \'gds.class\'")

    expect_error(getTableSNV(gds="HELLO", gdsSample=gdsF,
                        sampleCurrent="TEST", study.id="TEST",
                        minCov=12, minProb=0.32, eProb=0.001), error_message)
})

test_that("getTableSNV() must return error when gdsSample is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'gdsSample\' must be an object of class ",
                            "\'gdsn.class\' or \'gds.class\'")

    expect_error(getTableSNV(gds=gdsF, gdsSample="BABY",
                        sampleCurrent="TEST", study.id="TEST",
                        minCov=12, minProb=0.32, eProb=0.001), error_message)
})


context("computeLOHBlocksDNAChr() results")


test_that("computeLOHBlocksDNAChr() must return error when chr is vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'chr\' must be a single integer value representing ",
                                "a chromosome")

    expect_error(computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
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

    expect_error(computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
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

    expect_error(computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
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

    expect_error(computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
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

    expect_error(computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
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

    expect_error(computeLOHBlocksDNAChr(gds=gdsF, chrInfo=c("chr1", "chr2"),
                                        snp.pos=data.frame(), chr=2,
                                        genoN=1.33), error_message)
})


context("computeAllelicImbDNAChr() results")


test_that("computeAllelicImbDNAChr() must return error when chr is vector of numerics", {

    error_message <- paste0("The \'chr\' must be a single integer value representing ",
                            "a chromosome")

    expect_error(computeAllelicImbDNAChr(snp.pos=data.frame(), chr=c(1, 2),
                                            wAR=10, cutOffEmptyBox=-3),
                 error_message)
})


test_that("computeLOHBlocksDNAChr() must return error when chr is character string", {

    error_message <- paste0("The \'chr\' must be a single integer value representing ",
                            "a chromosome")

    expect_error(computeAllelicImbDNAChr(snp.pos=data.frame(), chr="titi",
                                            wAR=10, cutOffEmptyBox=-3),
                 error_message)
})



test_that("computeLOHBlocksDNAChr() must return error when wAR is character string", {

    error_message <- "The \'wAR\' must be a single numeric positive value."

    expect_error(computeAllelicImbDNAChr(snp.pos=data.frame(), chr=1,
                                            wAR="10", cutOffEmptyBox=-3),
                 error_message)
})



test_that("computeLOHBlocksDNAChr() must return error when wAR is vector of numerics", {

    error_message <- "The \'wAR\' must be a single numeric positive value."

    expect_error(computeAllelicImbDNAChr(snp.pos=data.frame(), chr=1,
                                            wAR=c(11, 21), cutOffEmptyBox=-3),
                 error_message)
})


context("computeAllelicFractionDNA() results")


test_that("computeAllelicFractionDNA() must return error when minCov is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'minCov\' must be a single numeric positive value"

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov="10", minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                 error_message)
})


test_that("computeAllelicFractionDNA() must return error when minCov is vextor of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'minCov\' must be a single numeric positive value"

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov=c(10, 11), minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                error_message)
})


test_that("computeAllelicFractionDNA() must return error when wAR is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'wAR\' must be a single numeric positive value."

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov=10, minProb=0.999,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
        wAR="BETA", verbose=FALSE), error_message)
})


test_that("computeAllelicFractionDNA() must return error when wAR is vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'wAR\' must be a single numeric positive value."

    expect_error(computeAllelicFractionDNA(gds=gdsFIle, gdsSample=gdsF,
                            sampleCurrent="test", study.id="test",
                            chrInfo=c("chr1", "chr2"), minCov=10, minProb=0.999,
                            eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3,
                            wAR=c(11, 21), verbose=FALSE),
                    error_message)
})


test_that("computeAllelicFractionDNA() must return error when minProb is vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                            "value between 0 and 1.")

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample="test",
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov=10, minProb=c(0.22, 0.999),
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                    error_message)
})


test_that("computeAllelicFractionDNA() must return error when minProb is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov=10, minProb="0.2",
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                    error_message)
})


test_that("computeAllelicFractionDNA() must return error when minProb is negative value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov=10, minProb=-0.11,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                    error_message)
})


test_that("computeAllelicFractionDNA() must return error when minProb is negative value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'minProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov=10, minProb=1.001,
        eProb=0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                    error_message)
})


test_that("computeAllelicFractionDNA() must return error when eProb is negative value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(computeAllelicFractionDNA(gds=gdsFIle, gdsSample="test",
            sampleCurrent="test", study.id="test",
            chrInfo=c("chr1", "chr2"), minCov=10, minProb=0.901,
            eProb=-0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                    error_message)
})


test_that("computeAllelicFractionDNA() must return error when eProb is vector of numerics", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov=10, minProb=0.999,
        eProb=c(0.22, 0.999), cutOffLOH=-5, cutOffHomoScore=-3, wAR=10,
        verbose=FALSE), error_message)
})


test_that("computeAllelicFractionDNA() must return error when eProb is character string", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(computeAllelicFractionDNA(gds=gdsF, gdsSample=gdsF,
        sampleCurrent="test", study.id="test", chrInfo=c("chr1", "chr2"),
        minCov=10, minProb=0.2, eProb="0.001", cutOffLOH=-5,
        cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                error_message)
})


test_that("computeAllelicFractionDNA() must return error when eProb is negative value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    error_message <- paste0("The \'eProb\' must be a single numeric positive ",
                                "value between 0 and 1.")

    expect_error(computeAllelicFractionDNA(gds=gdsFIle, gdsSample=gdsF,
        sampleCurrent="test", study.id="test",
        chrInfo=c("chr1", "chr2"), minCov=10, minProb=0.11,
        eProb=-0.001, cutOffLOH=-5, cutOffHomoScore=-3, wAR=10, verbose=FALSE),
                 error_message)
})


test_that("computeAllelicFractionDNA() must return error when verbose is numeric value", {

    data.dir <- system.file("extdata/tests", package="RAIDS")
    gdsFIle <- file.path(data.dir, "1KG_Test.gds")

    error_message <- paste0("The \'verbose\' parameters must be a ",
                                "single logical value (TRUE or FALSE).")

    expect_error(computeAllelicFractionDNA(gds=gdsFIle, gdsSample=gdsF,
        sampleCurrent="test", study.id="test", chrInfo=c("chr1", "chr2"),
        minCov=10, minProb=0.11,  eProb=0.001, cutOffLOH=-5,
        cutOffHomoScore=-3, wAR=10, verbose=22), error_message, fixed=TRUE)
})


context("computeAlleleFraction() results")


test_that("computeAlleleFraction() must return error when w is character string", {

    error_message <- "The \'w\' must be a single numeric positive value."

    expect_error(computeAlleleFraction(snp.pos=data.frame(), chr=1,
                            w="BETA", cutOff=-3), error_message)
})

test_that("computeAlleleFraction() must return error when w is a vector of numerics", {

    error_message <- "The \'w\' must be a single numeric positive value."

    expect_error(computeAlleleFraction(snp.pos=data.frame(), chr=1,
                                    w=c(11, 12), cutOff=-3), error_message)
})

test_that("computeAlleleFraction() must return error when w is negative numeric", {

    error_message <- "The \'w\' must be a single numeric positive value."

    expect_error(computeAlleleFraction(snp.pos=data.frame(), chr=1,
                                    w=-2, cutOff=-3), error_message)
})


test_that("computeAlleleFraction() must return error when chr is character string", {

    error_message <- paste0("The \'chr\' must be a single integer value representing ",
                                "a chromosome")

    expect_error(computeAlleleFraction(snp.pos=data.frame(), chr="1",
                                        w=11, cutOff=-3), error_message)
})



test_that("computeAlleleFraction() must return error when chr is a vector of numerics", {

    error_message <- paste0("The \'chr\' must be a single integer value representing ",
                                "a chromosome")

    expect_error(computeAlleleFraction(snp.pos=data.frame(), chr=c(1,2),
                                        w=11, cutOff=-3), error_message)
})

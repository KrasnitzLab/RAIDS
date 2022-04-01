### Unit tests for process1KG.R functions

library(aicsPaper)




#############################################################################
### Tests prepPed1KG() results
#############################################################################

context("prepPed1KG() results")


test_that("prepPed1KG() must return error when batch.v is a character string", {

    data.dir <- system.file("extdata", package="aicsPaper")

    pedDemoFile <- file.path(data.dir, "PedigreeDemo.ped")

    error_message <- "The batch.v must be an integer."

    expect_error(prepPed1KG(pedFile=pedFile, PATHGENO=data.dir,
                                batch.v="SAVE"), error_message)
})


test_that("prepPed1KG() must return error when batch.v is a float", {

    data.dir <- system.file("extdata", package="aicsPaper")

    pedDemoFile <- file.path(data.dir, "PedigreeDemo.ped")

    error_message <- "The batch.v must be an integer."

    expect_error(prepPed1KG(pedFile=pedFile, PATHGENO=data.dir,
                                batch.v=0.111), error_message)
})



test_that("prepPed1KG() must return error when PATHGENO is not existing", {

    data.dir <- system.file("extdata", package="aicsPaper")

    pedDemoFile <- file.path(data.dir, "PedigreeDemo.ped")

    notExisting <- paste0(data.dir, "/totoTest")

    error_message <- paste0("The path \'", notExisting, "\' does not exist.")

    expect_error(prepPed1KG(pedFile=pedDemoFile, PATHGENO=notExisting,
                                batch.v=0L), error_message)
})


test_that("prepPed1KG() must return error when pedigree file is not existing", {

    data.dir <- system.file("extdata", package="aicsPaper")

    pedDemoFile <- file.path(data.dir, "PedigreeDemoTOTO.ped")

    error_message <- paste0("The file \'", pedDemoFile, "\' does not exist.")

    expect_error(prepPed1KG(pedFile=pedDemoFile, PATHGENO=data.dir,
                            batch.v=0L), error_message)
})

test_that("prepPed1KG() must return the expected output", {

    data.dir <- system.file("extdata", package="aicsPaper")

    pedDemoFile <- file.path(data.dir, "PedigreeDemo.ped")

    expected <- data.frame(
        sample.id=c(paste0("HG0010", 0:9)),
        Name.ID=c(paste0("HG0010", 0:9)),
        sex=as.character(c(1,2,2,1,2,1,2,1,2,2)),
        pop.group=c(rep("ACB", 10)),
        superPop=c(rep("AFR", 10)),
        batch=c(rep(0L, 10)),
        stringsAsFactors=FALSE
    )
    row.names(expected) <- expected$sample.id

    expect_equal(prepPed1KG(pedFile=pedDemoFile, PATHGENO=data.dir,
                                batch.v=0L), expected)
})


#############################################################################
### Tests generateMapSnvSel() results
#############################################################################

context("generateMapSnvSel() results")

test_that("generateMapSnvSel() must return error when SNP file is not existing", {

    data.dir <- system.file("extdata", package="aicsPaper")

    snpFile <- file.path(data.dir, "SNP_TOTO.txt")

    outFile1 <- file.path(data.dir, "SNP_TOTO1.txt")

    outFile2 <- file.path(data.dir, "SNP_TOTO2.txt")

    error_message <- paste0("The file \'", snpFile, "\' does not exist.")

    expect_error(generateMapSnvSel(cutOff=0.01, fileSNV=snpFile,
                    fileLSNP=outFile1 , fileFREQ=outFile2), error_message)
})


test_that("generateMapSnvSel() must return error when cutOff file is a character string", {

    data.dir <- system.file("extdata", package="aicsPaper")

    snpFile <- file.path(data.dir, "PedigreeDemoTOTO.ped")

    outFile1 <- file.path(data.dir, "SNP_TOTO1.txt")

    outFile2 <- file.path(data.dir, "SNP_TOTO2.txt")

    error_message <- "The cutOff must be a single numeric value."

    expect_error(generateMapSnvSel(cutOff="CANADA", fileSNV=snpFile,
                        fileLSNP=outFile1 , fileFREQ=outFile2), error_message)
})


test_that("generateMapSnvSel() must return error when cutOff file is a array of numbers", {

    data.dir <- system.file("extdata", package="aicsPaper")

    snpFile <- file.path(data.dir, "PedigreeDemoTOTO.ped")

    outFile1 <- file.path(data.dir, "SNP_TOTO1.txt")

    outFile2 <- file.path(data.dir, "SNP_TOTO2.txt")

    error_message <- "The cutOff must be a single numeric value."

    expect_error(generateMapSnvSel(cutOff=c(0.01, 0.02), fileSNV=snpFile,
                        fileLSNP=outFile1 , fileFREQ=outFile2), error_message)
})


#############################################################################
### Tests generateGDS1KG() results
#############################################################################

context("generateGDS1KG() results")

test_that("generateGDS1KG() must return error when SNP file is not existing", {

    data.dir <- system.file("extdata", package="aicsPaper")

    fileNot <- file.path(data.dir, "TOTO_Not_Present.rds")

    outFile1 <- file.path(data.dir, "GDS_TEMP.gds")

    error_message <- paste0("The file \'", fileNot, "\' does not exist.")

    expect_error(generateGDS1KG(PATHGENO=data.dir,
                                fileNamePED=fileNot,
                                fileListSNP=fileNot,
                                fileSNPSel=fileNot, fileNameGDS=outFile1,
                                listSamples=NULL), error_message)
})

test_that("generateGDS1KG() must return error when PATHGENO is not existing", {

    data.dir <- system.file("extdata", package="aicsPaper")

    pedDemoFile <- file.path(data.dir, "PedigreeDemo.ped")

    notExisting <- paste0(data.dir, "/totoTest")

    error_message <- paste0("The path \'", notExisting, "\' does not exist.")

    expect_error(generateGDS1KG(PATHGENO=notExisting,
                                fileNamePED=pedDemoFile,
                                fileListSNP=pedDemoFile,
                                fileSNPSel=pedDemoFile, fileNameGDS=outFile1,
                                listSamples=NULL), error_message)
})

### Unit tests for process1KG.R functions

library(RAIDS)
library(withr)



#############################################################################
### Tests prepPed1KG() results
#############################################################################

context("prepPed1KG() results")


test_that("prepPed1KG() must return error when batch.v is a character string", {

    data.dir <- test_path("fixtures")

    pedDemoFile <- test_path("fixtures", "PedigreeDemo.ped")

    error_message <- "The batch.v must be an integer."

    expect_error(prepPed1KG(pedFile=pedDemoFile, PATHGENO=data.dir,
                                batch.v="SAVE"), error_message)
})


test_that("prepPed1KG() must return error when batch.v is a vector of float", {

    ddata.dir <- test_path("fixtures")

    pedDemoFile <- test_path("fixtures", "PedigreeDemo.ped")

    error_message <- "The batch.v must be an integer."

    expect_error(prepPed1KG(pedFile=pedDemoFile, PATHGENO=data.dir,
                                batch.v=c(0.111, 2)), error_message)
})



test_that("prepPed1KG() must return error when PATHGENO is not existing", {

    data.dir <- test_path("fixtures")

    pedDemoFile <- test_path("fixtures", "PedigreeDemo.ped")

    notExisting <- paste0(data.dir, "/totoTest")

    error_message <- paste0("The path \'", notExisting, "\' does not exist.")

    expect_error(prepPed1KG(pedFile=pedDemoFile, PATHGENO=notExisting,
                                batch.v=0L), error_message)
})


test_that("prepPed1KG() must return error when pedigree file is not existing", {

    data.dir <- test_path("fixtures")

    pedDemoFile <- test_path("fixtures", "PedigreeDemoTOTO.ped")

    error_message <- paste0("The file \'", pedDemoFile, "\' does not exist.")

    expect_error(prepPed1KG(pedFile=pedDemoFile, PATHGENO=data.dir,
                                batch.v=0L), error_message)
})

test_that("prepPed1KG() must return the expected output", {

    data.dir <- test_path("fixtures")

    pedDemoFile <- test_path("fixtures", "PedigreeDemo.ped")

    expected <- data.frame(
        sample.id=c(paste0("HG0010", 0:9)),
        Name.ID=c(paste0("HG0010", 0:9)),
        sex=as.character(c(1,2,2,1,2,1,2,1,2,2)),
        pop.group=c(rep("ACB", 10)),
        superPop=c(rep("AFR", 10)),
        batch=c(rep(0L, 10)),
        stringsAsFactors=FALSE)
    row.names(expected) <- expected$sample.id

    expect_equal(prepPed1KG(pedFile=pedDemoFile, PATHGENO=data.dir,
                                batch.v=0L), expected)
})


#############################################################################
### Tests generateMapSnvSel() results
#############################################################################

context("generateMapSnvSel() results")

test_that("generateMapSnvSel() must return error when SNP file is not existing", {

    data.dir <- test_path("fixtures")

    snpFile <- test_path("fixtures", "SNP_TOTO.txt")

    outFile1 <- test_path("fixtures", "SNP_TOTO1.rds")

    outFile2 <- test_path("fixtures", "SNP_TOTO2.rds")

    error_message <- paste0("The file \'", snpFile, "\' does not exist.")

    expect_error(generateMapSnvSel(cutOff=0.01, fileSNV=snpFile,
                    fileLSNP=outFile1 , fileFREQ=outFile2), error_message)
})


test_that("generateMapSnvSel() must return error when cutOff file is a character string", {

    data.dir <- test_path("fixtures")

    snpFile <- test_path("fixtures", "PedigreeDemoTOTO.ped")

    outFile1 <- test_path("fixtures", "SNP_TOTO1.rds")

    outFile2 <- test_path("fixtures", "SNP_TOTO2.rds")

    error_message <- "The cutOff must be a single numeric value."

    expect_error(generateMapSnvSel(cutOff="CANADA", fileSNV=snpFile,
                        fileLSNP=outFile1 , fileFREQ=outFile2), error_message)
})


test_that("generateMapSnvSel() must return error when cutOff file is a array of numbers", {

    data.dir <- test_path("fixtures")

    snpFile <- test_path("fixtures", "PedigreeDemoTOTO.ped")

    outFile1 <- test_path("fixtures", "SNP_TOTO1.rds")

    outFile2 <- test_path("fixtures", "SNP_TOTO2.rds")

    error_message <- "The cutOff must be a single numeric value."

    expect_error(generateMapSnvSel(cutOff=c(0.01, 0.02), fileSNV=snpFile,
                        fileLSNP=outFile1 , fileFREQ=outFile2), error_message)
})



test_that("generateMapSnvSel() must generate expected output", {

    data.dir <- test_path("fixtures")

    snvFile <- test_path("fixtures", "matFreqSNV_Demo.txt.bz2")

    ## Temporary output files
    snpIndexFile <- local_file(test_path("fixtures", "listSNP_TEMP01.rds"))
    filterSNVFile <- local_file(test_path("fixtures", "mapSNVSel_TEMP01.rds"))

    expect_equal(generateMapSnvSel(cutOff=c(0.01), fileSNV=snvFile,
                fileLSNP=snpIndexFile , fileFREQ=filterSNVFile), 0L)

    expect_true(file.exists(snpIndexFile))

    expect_true(file.exists(filterSNVFile))

    snpIndexesExpected <- c(seq_len(4), 6, 8, 9)

    expect_equal(readRDS(snpIndexFile), snpIndexesExpected)

    snpFilteredExpected <- data.frame(CHROM=rep("chr1", 7),
                            POS=c(16102, 51478, 51897, 51927,
                                        54489, 54707, 54715),
                            REF=c("T", "T", "C", "G", "G", "G", "C"),
                            ALT=c("G", "A", "A", "A", "A", "C", "T"),
                            AF=as.character(c(0.02, 0.11, 0.08, 0.07,
                                                0.1, 0.23, 0.21)),
                            EAS_AF=c("0.0", "0.0", "0.05", "0.01",
                                     "0.0", "0.08", "0.07"),
                            EUR_AF=c("0.04", "0.2", "0.12", "0.14",
                                     "0.18", "0.38", "0.34"),
                            AFR_AF=c("0.03", "0.02", "0.05", "0.06",
                                     "0.02", "0.18", "0.16"),
                            AMR_AF=c("0.03", "0.12", "0.06", "0.07",
                                     "0.1", "0.25", "0.24"),
                            SAS_AF=c("0.02", "0.22", "0.1", "0.09",
                                      "0.21", "0.28", "0.27"),
                            stringsAsFactors = FALSE)

    expect_equivalent(readRDS(filterSNVFile), snpFilteredExpected)
})



#############################################################################
### Tests generateGDS1KG() results
#############################################################################

context("generateGDS1KG() results")

test_that("generateGDS1KG() must return error when pedigree file does not exist", {

    data.dir <- test_path("fixtures")

    fileNot <- test_path("fixtures", "TOTO_Not_Present.rds")

    pedDemoFile <- test_path("fixtures", "PedigreeDemo.ped")

    outFile1 <- test_path("fixtures", "GDS_TEMP.gds")

    error_message <- paste0("The file \'", fileNot, "\' does not exist.")

    expect_error(generateGDS1KG(PATHGENO=data.dir,
                                fileNamePED=fileNot,
                                fileListSNP=pedDemoFile,
                                fileSNPSel=pedDemoFile, fileNameGDS=outFile1,
                                listSamples=NULL), error_message)
})

test_that("generateGDS1KG() must return error when PATHGENO is not existing", {

    data.dir <- test_path("fixtures")

    pedDemoFile <- test_path("fixtures", "PedigreeDemo.ped")

    notExisting <- paste0(data.dir, "/totoTest")

    error_message <- paste0("The path \'", notExisting, "\' does not exist.")

    expect_error(generateGDS1KG(PATHGENO=notExisting,
                                fileNamePED=pedDemoFile,
                                fileListSNP=pedDemoFile,
                                fileSNPSel=pedDemoFile, fileNameGDS=outFile1,
                                listSamples=NULL), error_message)
})

test_that("generateGDS1KG() must return error when SNP indexes file does not exist", {

    data.dir <- test_path("fixtures")

    pedDemoFile <- test_path("fixtures", "PedigreeDemo.ped")

    notExisting <- paste0(data.dir, "/totoTest")

    error_message <- paste0("The file \'", notExisting, "\' does not exist.")

    expect_error(generateGDS1KG(PATHGENO=data.dir,
                                fileNamePED=pedDemoFile,
                                fileListSNP=notExisting,
                                fileSNPSel=pedDemoFile, fileNameGDS=outFile1,
                                listSamples=NULL), error_message)
})

test_that("generateGDS1KG() must return error when SNP information file does not exist", {

    data.dir <- test_path("fixtures")

    pedDemoFile <- test_path("fixtures", "PedigreeDemo.ped")

    notExisting <- paste0(data.dir, "/totoTest")

    error_message <- paste0("The file \'", notExisting, "\' does not exist.")

    expect_error(generateGDS1KG(PATHGENO=data.dir,
                            fileNamePED=pedDemoFile,
                            fileListSNP=pedDemoFile,
                            fileSNPSel=notExisting, fileNameGDS=outFile1,
                            listSamples=NULL), error_message)
})


test_that("generateGDS1KG() must create a GDS file", {

    data.dir <- test_path("fixtures")

    pedigreeFile <- test_path("fixtures", "PedigreeDemo.rds")

    snpIndexFile <- test_path("fixtures", "listSNPIndexes_Demo.rds")

    filterSNVFile <- test_path("fixtures",  "mapSNVSelected_Demo.rds")

    ## Temporary GDS file containing 1KG information
    GDS_file <- local_file(test_path("fixtures", "1KG_TOTO.gds"))

    generateGDS1KG(PATHGENO=data.dir, fileNamePED=pedigreeFile,
                            fileListSNP=snpIndexFile,
                            fileSNPSel=filterSNVFile, fileNameGDS=GDS_file,
                            listSamples=NULL)

    expect_true(file.exists(GDS_file))
})



#############################################################################
### Tests identifyRelative() results
#############################################################################

context("identifyRelative() results")


test_that("identifyRelative() must return error when gds is character string", {

    fileIBDFile <- test_path("fixtures", "OUTPUT_01.rds")

    filePartFile <- test_path("fixtures", "OUTPUT_02.rds")

    error_message <- paste0("The \'gds\' parameter must be an object of ",
                                "class \'SNPGDSFileClass\'.")

    expect_error(identifyRelative(gds="test", maf=0.01,
            thresh=2^(-11/2), fileIBD=fileIBDFile, filePart=filePartFile),
            error_message)
})


test_that("identifyRelative() must return error when maf is a vector of numbers", {

    fileInput <- test_path("fixtures", "1KG_Demo.gds")

    fileIBDFile <- test_path("fixtures", "OUTPUT_01.rds")

    filePartFile <- test_path("fixtures", "OUTPUT_02.rds")

    error_message <- "The \'maf\' parameter must be a single numeric value."

    expect_error(identifyRelative(gds=fileInput, maf=c(0.01, 0.02),
                    thresh=2^(-11/2),
                    fileIBD=fileIBDFile, filePart=filePartFile), error_message)
})

test_that("identifyRelative() must return error when maf is a character strings", {

    fileInput <- test_path("fixtures", "1KG_Demo.gds")

    fileIBDFile <- test_path("fixtures", "OUTPUT_01.rds")

    filePartFile <- test_path("fixtures", "OUTPUT_02.rds")

    error_message <- "The \'maf\' parameter must be a single numeric value."

    expect_error(identifyRelative(gds=fileInput, maf="test",
                    thresh=2^(-11/2),
                    fileIBD=fileIBDFile, filePart=filePartFile), error_message)
})

test_that("identifyRelative() must return error when thresh is a character strings", {

    fileInput <- test_path("fixtures", "1KG_Demo.gds")

    fileIBDFile <- test_path("fixtures", "OUTPUT_01.rds")

    filePartFile <- test_path("fixtures", "OUTPUT_02.rds")

    error_message <- "The \'thresh\' parameter must be a single numeric value."

    expect_error(identifyRelative(gds=fileInput, maf=0.05,
                thresh="p-value",
                fileIBD=fileIBDFile, filePart=filePartFile), error_message)
})

test_that("identifyRelative() must return error when thresh is a vector of numerics", {

    fileInput <- test_path("fixtures", "1KG_Test.gds")

    fileIBDFile <- test_path("fixtures", "OUTPUT_01.rds")

    filePartFile <- test_path("fixtures", "OUTPUT_02.rds")

    error_message <- "The \'thresh\' parameter must be a single numeric value."

    expect_error(identifyRelative(gds=fileInput, maf=0.05,
                    thresh=c(0.01, 0.03),
                    fileIBD=fileIBDFile, filePart=filePartFile), error_message)
})


#############################################################################
### Tests addRef2GDS1KG() results
#############################################################################

context("addRef2GDS1KG() results")

test_that("addRef2GDS1KG() must return error when GDS file does not exist", {

    fileNot <- test_path("fixtures", "TOTO_GDS.gds")

    filePartFile <- test_path("fixtures", "Sample_Info_Test.RDS")

    error_message <- paste0("The file \'", fileNot, "\' does not exist.")

    expect_error(addRef2GDS1KG(fileNameGDS=fileNot, filePart=filePartFile),
                    error_message)
})


test_that("addRef2GDS1KG() must return error when RDS file does not exist", {

    fileNot <- test_path("fixtures", "TOTO_RDS.rds")

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    error_message <- paste0("The file \'", fileNot, "\' does not exist.")

    expect_error(addRef2GDS1KG(fileNameGDS=fileGDS, filePart=fileNot),
                 error_message)
})


#############################################################################
### Tests getRef1KGPop() results
#############################################################################

context("getRef1KGPop() results")


test_that("select1KGPop() must return error when gds is a character string", {

    error_message <- "The \'gds\' must be an object of class \'gds.class\'"

    expect_error(getRef1KGPop(gds="test.gds", popName="superPop"),
                    error_message)
})


test_that("select1KGPop() must return error when popName is a numeric", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'popName\' parameter must be a single character string."

    expect_error(getRef1KGPop(gds=gdsF, popName=33), error_message)
})


test_that("select1KGPop() must return error when popName is a numeric", {

    data.dir <- test_path("fixtures")

    gdsFIle <- test_path("fixtures", "1KG_Test.gds")

    gdsF <- openfn.gds(gdsFIle)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir = parent.frame())

    error_message <- "The \'popName\' parameter must be a single character string."

    expect_error(getRef1KGPop(gds=gdsF, popName=33), error_message)
})


test_that("select1KGPop() must return expected result", {

    gdsFIle <- test_path("fixtures", "1KG_TEMP_001.gds")

    gds_1KG <- local_GDS_1KG_file(gdsFIle, env=parent.frame())

    results <- getRef1KGPop(gds_1KG, popName="superPop")

    ## Close the GDS file
    ## The file will automatically be deleted
    closefn.gds(gds_1KG)

    expected <- c("EUR",  "EUR")
    names(expected) <- c("HTT101", "HTT103")

    expect_identical(results, expected)
})




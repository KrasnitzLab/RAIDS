### Unit tests for tools.R functions

library(RAIDS)
library(withr)
library(testthat)

#############################################################################
### Tests snvListVCF() results
#############################################################################

context("snvListVCF() results")


test_that("snvListVCF() must return error when offset is a character string", {

    data.dir <- system.file("extdata", package="RAIDS")

    fileGDS <- local_file(file.path(data.dir, "GDS_TEMP2.gds"))

    tmpFileGDS <- createfn.gds(fileGDS)
    put.attr.gdsn(tmpFileGDS$root, "FileFormat", "SNP_ARRAY")
    add.gdsn(tmpFileGDS, "sample.id", c("test01", "test02"))

    samp.annot <- data.frame(sex=c("1", "2"),
                                pop.group=c("ACB", "ACB"),
                                superPop=c("AFR", "AFR"),
                                batch=c("0", "0"),
                                stringsAsFactors=FALSE)
    add.gdsn(tmpFileGDS, "sample.annot", samp.annot)


    add.gdsn(tmpFileGDS, "snp.id", paste0("s",seq_len(2)))
    add.gdsn(tmpFileGDS, "snp.chromosome", c(1L, 1L), storage = "uint16")
    add.gdsn(tmpFileGDS, "snp.position", c(10L, 12L), storage = "int32")
    add.gdsn(tmpFileGDS, "snp.allele", c("0/1", "1/1"))
    add.gdsn(tmpFileGDS, "genotype", valdim=c(2, 2), storage="bit2")

    closefn.gds(tmpFileGDS)

    gds <- snpgdsOpen(fileGDS)

    error_message <- "The \'offset\' must be a single integer."

    expect_error(snvListVCF(gdsReference=gds, fileOut, offset="HELLO",
                freqCutoff=NULL), error_message)

    closefn.gds(gds)

    deferred_run()
})


test_that("snvListVCF() must return error when gds is a character string", {

    data.dir <- system.file("extdata", package="RAIDS")

    fileOut <- file.path(data.dir, "VCF_TEMP.vcf")

    error_message <- "The \'gdsReference\' must be an object of class \'gds.class\'."

    expect_error(snvListVCF(gdsReference="welcome.txt", fileOut=fileOut, offset=0L,
                            freqCutoff=NULL), error_message)
})


test_that("snvListVCF() must return error when freqCutoff is a character string", {

    data.dir <- test_path("fixtures")

    fileGDS <- local_file(file.path(data.dir, "GDS_TEMP3.gds"))
    withr::defer(unlink(fileGDS, force=TRUE), envir = parent.frame())

    tmpFileGDS <- createfn.gds(fileGDS)
    put.attr.gdsn(tmpFileGDS$root, "FileFormat", "SNP_ARRAY")
    add.gdsn(tmpFileGDS, "sample.id", c("test01", "test02"))

    samp.annot <- data.frame(sex=c("1", "2"),
                                pop.group=c("ACB", "ACB"),
                                superPop=c("AFR", "AFR"),
                                batch=c("0", "0"),
                                stringsAsFactors=FALSE)
    add.gdsn(tmpFileGDS, "sample.annot", samp.annot)

    add.gdsn(tmpFileGDS, "snp.id", paste0("s", seq_len(2)))
    add.gdsn(tmpFileGDS, "snp.chromosome", c(1L, 1L), storage = "uint16")
    add.gdsn(tmpFileGDS, "snp.position", c(10L, 12L), storage = "int32")
    add.gdsn(tmpFileGDS, "snp.allele", c("0/1", "1/1"))
    add.gdsn(tmpFileGDS, "genotype", valdim=c(2,2), storage="bit2")

    closefn.gds(tmpFileGDS)

    gds <- snpgdsOpen(fileGDS)
    withr::defer(closefn.gds(gds), envir = parent.frame())

    fileOut <- file.path(data.dir, "VCF_TEMP.vcf")

    error_message <- "The \'freqCutoff\' must be a single numeric or NULL."

    expect_error(snvListVCF(gdsReference=gds, fileOut=fileOut, offset=0L,
                                freqCutoff="BED"), error_message)
})


test_that("snvListVCF() must return expected results when freqCutoff is NULL", {

    data.dir <- test_path("fixtures")

    fileGDS <- file.path(data.dir, "1KG_Test.gds")

    gds <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gds), envir=parent.frame())

    fileOut <- file.path(data.dir, "VCF_TEMP_01.vcf")
    withr::defer(unlink(fileOut, force=TRUE), envir=parent.frame())

    result1 <- suppressWarnings(snvListVCF(gdsReference=gds, fileOut=fileOut,
                        offset=0L, freqCutoff=NULL))

    ## Read two times the vcf file,
    ## First for the columns names, second for the data
    tmp_vcf <- readLines(fileOut)
    tmp_vcf_data <- read.table(fileOut, stringsAsFactors=FALSE)

    # filter for the columns names
    tmp_vcf <- tmp_vcf[-(grep("#CHROM", tmp_vcf)+1):-(length(tmp_vcf))]
    vcf_names <- unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
    names(tmp_vcf_data) <- vcf_names


    expect_equal(result1, 0L)
    expect_true(file.exists(fileOut))
    expect_equal(nrow(tmp_vcf_data), 7)
    expect_equal(ncol(tmp_vcf_data), 8)
})


test_that("snvListVCF() must return expected results when freqCutoff is 0.3", {

    data.dir <- test_path("fixtures")

    fileGDS <- file.path(data.dir, "1KG_Test.gds")

    gds <- openfn.gds(fileGDS)
    withr::defer(closefn.gds(gds), envir=parent.frame())

    fileOut <- file.path(data.dir, "VCF_TEMP_02.vcf")
    withr::defer(unlink(fileOut, force=TRUE), envir=parent.frame())

    result1 <- suppressWarnings(snvListVCF(gdsReference=gds, fileOut=fileOut, offset=0L,
                                           freqCutoff=0.3))

    ## Read two times the vcf file,
    ## First for the columns names, second for the data
    tmp_vcf <- readLines(fileOut)
    tmp_vcf_data <- read.table(fileOut, stringsAsFactors=FALSE)

    # filter for the columns names
    tmp_vcf <- tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
    vcf_names <- unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
    names(tmp_vcf_data) <- vcf_names


    expect_equal(result1, 0L)
    expect_true(file.exists(fileOut))
    expect_equal(nrow(tmp_vcf_data), 2)
    expect_equal(ncol(tmp_vcf_data), 8)
})


#############################################################################
### Tests groupChr1KGSNV() results
#############################################################################

context("groupChr1KGSNV() results")

test_that("groupChr1KGSNV() must return error when pathGenoChr does not exist", {

    data.dir <- system.file("extdata", package="RAIDS")

    dirNotExisting <- file.path(data.dir, "RED_TOMATO")

    error_message <- paste0("The path \'", dirNotExisting,
                                "\' does not exist.")

    expect_error(groupChr1KGSNV(pathGenoChr=dirNotExisting,
                                    pathOut=data.dir), error_message)
})


test_that("groupChr1KGSNV() must return error when pathOut does not exist", {

    data.dir <- system.file("extdata", package="RAIDS")

    dirNotExisting <- file.path(data.dir, "/RED_TOMATO")

    error_message <- paste0("The path \'", dirNotExisting,
                            "\' does not exist.")

    expect_error(groupChr1KGSNV(pathGenoChr=data.dir,
                                    pathOut=dirNotExisting), error_message)
})


test_that("groupChr1KGSNV() must return expected results", {

    dataDir <- test_path("fixtures/demoGenoChr")

    dirTMP <- withr::local_tempdir(tmpdir="test_groupChr1KGSNV")

    result <- groupChr1KGSNV(pathGenoChr=dataDir, pathOut=dirTMP)

    expect_equal(result, 0L)

    result2 <- list.files(dirTMP)

    expect_true(all(c("NA12003.csv.bz2", "NA12004.csv.bz2",
                        "NA12005.csv.bz2", "NA12006.csv.bz2") %in% result2))

    geno1 <- read.csv2(file.path(dirTMP, "NA12003.csv.bz2"),
                                        sep="\t", row.names=NULL)

    expect_equal(nrow(geno1), 2178L)
    expect_equal(colnames(geno1), "NA12003")
    expect_equal(geno1[2,], "0|0")
    expect_equal(geno1[3,], "0|0")
    expect_equal(geno1[444,], "0|0")
    expect_equal(geno1[2177,], "0|1")
    expect_equal(geno1[2082,], "1|1")

    geno2 <- read.csv2(file.path(dirTMP, "NA12006.csv.bz2"),
                                    sep="\t", row.names=NULL)

    expect_equal(nrow(geno2), 2178L)
    expect_equal(colnames(geno2), "NA12006")
    expect_equal(geno2[2,], "0|0")
    expect_equal(geno2[3,], "0|1")
    expect_equal(geno2[444,], "0|0")
    expect_equal(geno2[2177,], "0|0")
    expect_equal(geno2[2082,], "1|1")

    geno3 <- read.csv2(file.path(dirTMP, "NA12005.csv.bz2"),
                                    sep="\t", row.names=NULL)

    expect_equal(nrow(geno3), 2178L)
    expect_equal(colnames(geno3), "NA12005")
    expect_equal(geno3[2,], "1|0")
    expect_equal(geno3[3,], "0|0")
    expect_equal(geno3[444,], "0|0")
    expect_equal(geno3[2177,], "0|0")
    expect_equal(geno3[2082,], "1|1")
})

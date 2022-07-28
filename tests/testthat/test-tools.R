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

    gdsFile <- local_file(file.path(data.dir, "GDS_TEMP2.gds"))

    gdsFileTMP <- createfn.gds(gdsFile)
    put.attr.gdsn(gdsFileTMP$root, "FileFormat", "SNP_ARRAY")
    add.gdsn(gdsFileTMP, "sample.id", c("test01", "test02"))

    samp.annot <- data.frame(sex=c("1","2"),
                             pop.group=c("ACB", "ACB"),
                             superPop=c("AFR", "AFR"),
                             batch=c("0", "0"),
                             stringsAsFactors=FALSE)
    add.gdsn(gdsFileTMP, "sample.annot", samp.annot)


    add.gdsn(gdsFileTMP, "snp.id", paste0("s",seq_len(2)))
    add.gdsn(gdsFileTMP, "snp.chromosome", c(1L, 1L), storage = "uint16")
    add.gdsn(gdsFileTMP, "snp.position", c(10L, 12L), storage = "int32")
    add.gdsn(gdsFileTMP, "snp.allele", c("0/1", "1/1"))
    add.gdsn(gdsFileTMP, "genotype", valdim=c(2,2), storage="bit2")

    closefn.gds(gdsFileTMP)

    gds <- snpgdsOpen(gdsFile)

    error_message <- "The \'offset\' must be a single integer."

    expect_error(snvListVCF(gds=gds, fileOUT, offset="HELLO",
                freqCutoff=NULL), error_message)

    closefn.gds(gds)

    deferred_run()
})


test_that("snvListVCF() must return error when gds is a character string", {

    data.dir <- system.file("extdata", package="RAIDS")

    fileOUT <- file.path(data.dir, "VCF_TEMP.vcf")

    error_message <- paste0("The \'gds\' parameter must be an object of ",
                                "class \'SNPGDSFileClass\'.")

    expect_error(snvListVCF(gds="welcome.txt", fileOUT=fileOUT, offset=0L,
                            freqCutoff=NULL), error_message)
})


test_that("snvListVCF() must return error when freqCutoff is a character string", {

    data.dir <- system.file("extdata", package="RAIDS")

    gdsFile <- local_file(file.path(data.dir, "GDS_TEMP3.gds"))

    gdsFileTMP <- createfn.gds(gdsFile)
    put.attr.gdsn(gdsFileTMP$root, "FileFormat", "SNP_ARRAY")
    add.gdsn(gdsFileTMP, "sample.id", c("test01", "test02"))

    samp.annot <- data.frame(sex=c("1","2"),
                             pop.group=c("ACB", "ACB"),
                             superPop=c("AFR", "AFR"),
                             batch=c("0", "0"),
                             stringsAsFactors=FALSE)
    add.gdsn(gdsFileTMP, "sample.annot", samp.annot)

    add.gdsn(gdsFileTMP, "snp.id", paste0("s",seq_len(2)))
    add.gdsn(gdsFileTMP, "snp.chromosome", c(1L, 1L), storage = "uint16")
    add.gdsn(gdsFileTMP, "snp.position", c(10L, 12L), storage = "int32")
    add.gdsn(gdsFileTMP, "snp.allele", c("0/1", "1/1"))
    add.gdsn(gdsFileTMP, "genotype", valdim=c(2,2), storage="bit2")

    closefn.gds(gdsFileTMP)

    gds <- snpgdsOpen(gdsFile)

    fileOUT <- file.path(data.dir, "VCF_TEMP.vcf")

    error_message <- "The \'freqCutoff\' must be a single numeric or NULL."

    expect_error(snvListVCF(gds=gds, fileOUT=fileOUT, offset=0L,
                            freqCutoff="BED"), error_message)
})



#############################################################################
### Tests groupChr1KGSNV() results
#############################################################################

context("groupChr1KGSNV() results")

test_that("groupChr1KGSNV() must return error when PATHGENOCHR does not exist", {

    data.dir <- system.file("extdata", package="RAIDS")

    dirNotExisting <- file.path(data.dir, "RED_TOMATO")

    error_message <- paste0("The path \'", dirNotExisting,
                                "\' does not exist.")

    expect_error(groupChr1KGSNV(PATHGENOCHR=dirNotExisting,
                                    PATHOUT=data.dir), error_message)
})

test_that("groupChr1KGSNV() must return error when PATHOUT does not exist", {

    data.dir <- system.file("extdata", package="RAIDS")

    dirNotExisting <- file.path(data.dir, "/RED_TOMATO")

    error_message <- paste0("The path \'", dirNotExisting,
                            "\' does not exist.")

    expect_error(groupChr1KGSNV(PATHGENOCHR=data.dir,
                                PATHOUT=dirNotExisting), error_message)
})


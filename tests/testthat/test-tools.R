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


#############################################################################
### Tests getMatrixPopSynthetic() results
#############################################################################

context("getMatrixPopSynthetic() results")


test_that("getMatrixPopSynthetic() must return one column per subcontinental population", {

    pRAIDS <- local_admixture_pRAIDS()

    result <- getMatrixPopSynthetic(pRAIDS)

    gdsRef <- openfn.gds(fileReferenceGDS(pRAIDS))
    sampleId <- read.gdsn(index.gdsn(gdsRef, "sample.id"))
    sampleAnnot <- read.gdsn(index.gdsn(gdsRef, "sample.annot"))
    closefn.gds(gdsRef)

    gdsProfile <- openfn.gds(file.path(pathProfileGDS(pRAIDS), "ex1.gds"))
    studyAnnot <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    closefn.gds(gdsProfile)
    caseSyn <- studyAnnot$case.id[studyAnnot$study.id == "MYDATA.Synthetic"]

    ## The demo has 52 synthetic profiles from 26 populations
    expect_equal(dim(result), c(2, 26))
    expect_setequal(c(result), caseSyn)

    ## All profiles in a column are from the same population
    popGr <- matrix(sampleAnnot$pop.group[match(result, sampleId)], nrow=2)
    expect_true(all(popGr[1, ] == popGr[2, ]))
    expect_equal(length(unique(popGr[1, ])), 26)
    expect_true(noOpenGDS())
})


#############################################################################
### Tests prepPMatrix() results
#############################################################################

context("prepPMatrix() results")


test_that("prepPMatrix() must return the allele frequencies of each ancestry group", {

    pRAIDS <- local_admixture_pRAIDS()
    matGr <- getMatrixPopSynthetic(pRAIDS)

    gdsRef <- openfn.gds(fileReferenceGDS(pRAIDS))
    anchor <- read.gdsn(index.gdsn(gdsRef, "anchor.ref"),
                            start=c(1, 1), count=c(-1, 1))
    sampleId <- read.gdsn(index.gdsn(gdsRef, "sample.id"))
    closefn.gds(gdsRef)

    gdsProfile <- openfn.gds(file.path(pathProfileGDS(pRAIDS), "ex1.gds"))
    pruned <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    closefn.gds(gdsProfile)

    result <- quietly(prepPMatrix(anchor=anchor, snpId=pruned, matGr=matGr,
                                    pRAIDS=pRAIDS))

    ## One entry per ancestry group, in the order found in anchor
    expect_equal(names(result), setdiff(unique(anchor), "-"))
    expect_true(noOpenGDS())

    gdsRef <- snpgdsOpen(fileReferenceGDS(pRAIDS))
    withr::defer(closefn.gds(gdsRef))

    alleleFreq <- function(samples) {
        geno <- quietly(snpgdsGetGeno(gdsRef, sample.id=samples,
                                        snp.id=pruned, snpfirstdim=FALSE))
        return(colSums(geno) / (2 * nrow(geno)))
    }

    for (ag in names(result)) {
        anchorSamples <- sampleId[anchor == ag]

        ## First frequency uses all anchor profiles of the group
        expect_equal(result[[ag]]$f[[1]], alleleFreq(anchorSamples))

        ## Other frequencies exclude the profiles of each synthetic group
        expect_equal(length(result[[ag]]$index), nrow(matGr))
        for (i in seq_len(nrow(matGr))) {
            kept <- anchorSamples[!anchorSamples %in% matGr[i, ]]
            expect_equal(result[[ag]]$f[[result[[ag]]$index[i]]],
                            alleleFreq(kept))
        }
    }
})


test_that("prepPMatrix() must return only the global frequency when matGr is NULL", {

    pRAIDS <- local_admixture_pRAIDS()

    gdsRef <- openfn.gds(fileReferenceGDS(pRAIDS))
    anchor <- read.gdsn(index.gdsn(gdsRef, "anchor.ref"),
                            start=c(1, 1), count=c(-1, 1))
    pruned <- read.gdsn(index.gdsn(gdsRef, "snp.id"))[seq_len(20)]
    closefn.gds(gdsRef)

    result <- quietly(prepPMatrix(anchor=anchor, snpId=pruned, matGr=NULL,
                                    pRAIDS=pRAIDS))

    expect_equal(length(result), 5)
    for (ag in names(result)) {
        expect_equal(length(result[[ag]]$f), 1)
        expect_equal(length(result[[ag]]$f[[1]]), 20)
        expect_null(result[[ag]]$index)
    }
    expect_true(noOpenGDS())
})


#############################################################################
### Tests writePMatrix() results
#############################################################################

context("writePMatrix() results")


test_that("writePMatrix() must write the P matrix used by each synthetic group", {

    pRAIDS <- local_admixture_pRAIDS()
    pathOut <- withr::local_tempdir()

    ## Two ancestry groups, 3 SNVs and 2 synthetic groups
    ## Each synthetic group removes profiles from agA
    rowF <- list(agA=list(f=list(c(0.1, 0.2, 0.3), c(0.15, 0.25, 0.35),
                                    c(0.12, 0.22, 0.32)),
                            index=c(2L, 3L)),
                    agB=list(f=list(c(0.5, 0.6, 0.7), c(0.55, 0.65, 0.75)),
                            index=c(1L, 2L)))

    result <- writePMatrix(pathOut=pathOut, rowF=rowF, pRAIDS=pRAIDS)

    expect_equal(result, c(1L, 2L))

    readP <- function(name) {
        unname(as.matrix(read.table(file.path(pathOut, name))))
    }

    expect_equal(readP("ex1.P.in"),
                    1 - cbind(rowF$agA$f[[1]], rowF$agB$f[[1]]))
    expect_equal(readP(paste0("ex1.syn.", result[1], ".P.in")),
                    1 - cbind(rowF$agA$f[[2]], rowF$agB$f[[1]]))
    expect_equal(readP(paste0("ex1.syn.", result[2], ".P.in")),
                    1 - cbind(rowF$agA$f[[3]], rowF$agB$f[[2]]))
})


test_that("writePMatrix() must give a valid P matrix to a group without anchor profiles", {

    pRAIDS <- local_admixture_pRAIDS()
    pathOut <- withr::local_tempdir()

    ## The first synthetic group removes profiles from agA
    ## The second synthetic group has no anchor profile so it must use
    ## the P matrix with all anchor profiles
    rowF <- list(agA=list(f=list(c(0.1, 0.2, 0.3), c(0.15, 0.25, 0.35)),
                            index=c(2L, 1L)),
                    agB=list(f=list(c(0.5, 0.6, 0.7)),
                            index=c(1L, 1L)))

    result <- writePMatrix(pathOut=pathOut, rowF=rowF, pRAIDS=pRAIDS)

    readP <- function(name) {
        unname(as.matrix(read.table(file.path(pathOut, name))))
    }

    ## prepSynADMIXTURE() uses the file 'ex1.syn.<index>.P.in' and
    ## the file 'ex1.P.in' when the index is 0
    expect_equal(result, c(1L, 0L))
    expect_equal(readP("ex1.syn.1.P.in"),
                    1 - cbind(rowF$agA$f[[2]], rowF$agB$f[[1]]))
    expect_equal(readP("ex1.P.in"),
                    1 - cbind(rowF$agA$f[[1]], rowF$agB$f[[1]]))
    expect_false(file.exists(file.path(pathOut, "ex1.syn.2.P.in")))
})


#############################################################################
### Tests writeBimPruned(), writeBedProfile() and writeFamProfile() results
#############################################################################

context("writeBimPruned(), writeBedProfile() and writeFamProfile() results")


test_that("the PLINK files of the profile must have the same SNVs and samples", {

    pRAIDS <- local_admixture_pRAIDS()
    pathOut <- withr::local_tempdir()
    prefix <- file.path(pathOut, "ex1")

    quietly(writeBimPruned(fileOut=paste0(prefix, ".bim"), pRAIDS=pRAIDS))
    quietly(writeBedProfile(fileOut=paste0(prefix, ".bed"),
                    listProfile="ex1", profileOnly=TRUE, pRAIDS=pRAIDS))
    quietly(writeFamProfile(fileOut=paste0(prefix, ".fam"),
                    listProfile="ex1", profileOnly=TRUE, pRAIDS=pRAIDS))

    result <- genio::read_plink(prefix, verbose=FALSE)

    gdsProfile <- snpgdsOpen(file.path(pathProfileGDS(pRAIDS), "ex1.gds"))
    pruned <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    geno <- quietly(snpgdsGetGeno(gdsProfile, sample.id="ex1",
                                    snp.id=pruned, snpfirstdim=TRUE))
    closefn.gds(gdsProfile)

    expect_equal(dim(result$X), c(61, 1))
    expect_equal(result$bim$id, pruned)
    expect_equal(result$fam$id, "ex1")
    expect_equal(unname(result$X[, 1]), as.numeric(geno[, 1]))
    expect_true(noOpenGDS())
})


test_that("the PLINK files with the reference must remove the profiles in listRM", {

    pRAIDS <- local_admixture_pRAIDS()
    pathOut <- withr::local_tempdir()
    prefix <- file.path(pathOut, "ex1")
    listRM <- c("HG00246", "HG00138")

    quietly(writeBimPruned(fileOut=paste0(prefix, ".bim"), pRAIDS=pRAIDS))
    quietly(writeBedProfile(fileOut=paste0(prefix, ".bed"),
                    listProfile="ex1", listRM=listRM, pRAIDS=pRAIDS))
    quietly(writeFamProfile(fileOut=paste0(prefix, ".fam"),
                    listProfile="ex1", listRM=listRM, pRAIDS=pRAIDS))

    result <- genio::read_plink(prefix, verbose=FALSE)

    ## 156 reference profiles - 2 removed + 1 profile
    expect_equal(dim(result$X), c(61, 155))
    expect_false(any(listRM %in% result$fam$id))
    expect_equal(tail(result$fam$id, 1), "ex1")
    expect_true(noOpenGDS())
})

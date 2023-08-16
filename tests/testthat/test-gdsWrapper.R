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
    add.gdsn(GDS_file_tmp, "lap", rep(0.5, 12))
    sync.gds(GDS_file_tmp)

    return(GDS_file_tmp)
}


#############################################################################
### Tests generateGDSRefSample() results
#############################################################################

context("generateGDSRefSample() results")


test_that("generateGDSRefSample() must copy the expected entry in \"sample.annot\" node of the GDS file", {

    ## Create a temporary GDS file in an test directory
    dataDir <- system.file("extdata/tests", package="RAIDS")
    fileGDS <- file.path(dataDir, "GDS_TEMP_06.gds")

    ## Create and open a temporary GDS file
    GDS_file_tmp  <- local_GDS_file(fileGDS)

    ## Create PEF information
    pedInformation <- data.frame(sample.id=c("sample_01", "sample_02", "sample_03"),
                            Name.ID=c("sample_01", "sample_02", "sample_03"),
                            sex=c(1, 1, 2),  # 1:Male  2: Female
                            pop.group=c("ACB", "ACB", "ACB"),
                            superPop=c("AFR", "AFR", "AFR"),
                            batch=c(1, 1, 1), stringsAsFactors=FALSE)
    rownames(pedInformation) <- pedInformation$Name.ID

    ## Add samples to the GDS file
    results3 <- RAIDS:::generateGDSRefSample(gdsReference=GDS_file_tmp,  dfPedReference=pedInformation,
                                listSamples=c("sample_01", "sample_02"))

    ## Read sample names from GDS file
    results1 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.id"))

    results2 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.annot"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expected1 <- c("sample_01", "sample_02")

    expected2 <- pedInformation[c("sample_01", "sample_02"),]
    expected2$Name.ID <- NULL
    expected2$sample.id <- NULL
    rownames(expected2) <- NULL

    expect_equal(results1, expected1)
    expect_equal(results2, expected2)
    expect_equal(results3, expected1)
})



#############################################################################
### Tests addBlockInGDSAnnot() results
#############################################################################

context("addBlockInGDSAnnot() results")


test_that("addBlockInGDSAnnot() must return expected result", {

    ## Create and open a temporary GDS Annotation file
    GDS_path  <- test_path("fixtures", "GDS_addBlockInGDSAnnot_Temp_01.gds")
    GDS_file_tmp  <- createfn.gds(filename=GDS_path)
    defer(unlink(x=GDS_path, force=TRUE), envir=parent.frame())

    blockType <- "EAS.0.05.500k"
    blockDescription <- "EAS population blocks based on 500k windows"

    entries <- c(1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5)

    ## Add block to the GDS file
    results1 <- RAIDS:::addBlockInGDSAnnot(gds=GDS_file_tmp, listBlock=entries,
            blockName=blockType, blockDesc=blockDescription)


    blockType2 <- "AFR.0.05.500k"
    blockDescription2 <- "AFR population blocks based on 500k windows"

    entries2 <- c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5)

    ## Add block to the GDS file
    results2 <- RAIDS:::addBlockInGDSAnnot(gds=GDS_file_tmp, listBlock=entries2,
            blockName=blockType2, blockDesc=blockDescription2)


    ## Read block.annot from GDS file
    results3 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="block.annot"))

    ## Read block.annot from GDS file
    results4 <- read.gdsn(index.gdsn(node=GDS_file_tmp, path="block"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(results1, 0L)
    expect_equal(results2, 0L)
    expect_equal(results3, data.frame(block.id=c(blockType, blockType2),
                            block.desc=c(blockDescription, blockDescription2),
                            stringsAsFactors=FALSE))

    expect_equal(results4, matrix(data=c(entries, entries2), ncol=2,
                                    byrow=FALSE))
})

#############################################################################
### Tests generateGDSgenotype() results
#############################################################################

context("generateGDSgenotype() results")

test_that("generateGDSgenotype() must return expected result", {

    ## Create and open a temporary GDS Annotation file
    GDS_path  <- test_path("fixtures", "GDS_generateGDSgenotype_Temp_01.gds")
    GDS_file_tmp  <- createfn.gds(filename=GDS_path)
    defer(unlink(x=GDS_path, force=TRUE), envir=parent.frame())

    put.attr.gdsn(GDS_file_tmp$root, "FileFormat", "SNP_ARRAY")

    pedigree <- data.frame(sample.id=c("HG00100", "HG00101", "HG00102"),
                            Name.ID=c("HG00100", "HG00101", "HG00102"),
                            sex=c(1,2,2), pop.group=c("ACB", "ACB", "ACB"),
                            superPop=c("AFR", "AFR", "AFR"), batch=c(0, 0, 0),
                            stringsAsFactors=FALSE)

    rownames(pedigree) <- pedigree$sample.id

    filterSNVFile <- test_path("fixtures", "mapSNVSelected_Demo.rds")

    ## Add information about samples
    listSampleGDS <- RAIDS:::generateGDSRefSample(gdsReference=GDS_file_tmp,
                dfPedReference=pedigree, listSamples=NULL)

    ## Add SNV information to the Reference GDS
    RAIDS:::generateGDSSNPinfo(gdsReference=GDS_file_tmp,
                    fileFreq=filterSNVFile, verbose=FALSE)

    snpIndexFile <- test_path("fixtures", "listSNPIndexes_Demo.rds")

    ## Add genotype information to the Reference GDS
    result1 <- RAIDS:::generateGDSgenotype(gds=GDS_file_tmp,
        pathGeno=test_path("fixtures"),
        fileSNPsRDS=snpIndexFile, listSamples=listSampleGDS, verbose=FALSE)

    result2 <- read.gdsn(index.gdsn(GDS_file_tmp, "genotype"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(result1, 0L)
    expect_equal(result2, matrix(data=c(0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1,
                        1, 0, 0, 2, 2, 0, 2, 1), nrow=7, byrow=FALSE))

})


#############################################################################
### Tests appendGDSgenotype() results
#############################################################################

context("appendGDSgenotype() results")

test_that("appendGDSgenotype() must return expected result", {

    ## Create and open a temporary GDS Annotation file
    GDS_path  <- test_path("fixtures", "GDS_appendGDSgenotype_Temp_01.gds")
    GDS_file_tmp  <- createfn.gds(filename=GDS_path)
    defer(unlink(x=GDS_path, force=TRUE), envir=parent.frame())

    put.attr.gdsn(GDS_file_tmp$root, "FileFormat", "SNP_ARRAY")

    pedigree <- data.frame(sample.id=c("HG00100", "HG00101", "HG00102"),
                           Name.ID=c("HG00100", "HG00101", "HG00102"),
                           sex=c(1,2,2), pop.group=c("ACB", "ACB", "ACB"),
                           superPop=c("AFR", "AFR", "AFR"), batch=c(0, 0, 0),
                           stringsAsFactors=FALSE)

    rownames(pedigree) <- pedigree$sample.id

    filterSNVFile <- test_path("fixtures", "mapSNVSelected_Demo.rds")

    ## Add information about samples
    listSampleGDS <- RAIDS:::generateGDSRefSample(gdsReference=GDS_file_tmp,
                dfPedReference=pedigree, listSamples=NULL)

    ## Add SNV information to the Reference GDS
    RAIDS:::generateGDSSNPinfo(gdsReference=GDS_file_tmp,
                               fileFreq=filterSNVFile, verbose=FALSE)

    snpIndexFile <- test_path("fixtures", "listSNPIndexes_Demo.rds")

    ## Add genotype information to the Reference GDS
    result1 <- RAIDS:::generateGDSgenotype(gds=GDS_file_tmp,
        pathGeno=test_path("fixtures"),
        fileSNPsRDS=snpIndexFile, listSamples=listSampleGDS[1], verbose=FALSE)


    ## Append genotype information to the Reference GDS
    result1 <- RAIDS:::appendGDSgenotype(gds=GDS_file_tmp,
                    pathGeno=test_path("fixtures"),
                    fileSNPsRDS=snpIndexFile, listSample=listSampleGDS[2],
                    verbose=FALSE)

    result2 <- read.gdsn(index.gdsn(GDS_file_tmp, "genotype"))

    ## Close GDS file
    ## The file will automatically be deleted
    closefn.gds(gdsfile=GDS_file_tmp)

    expect_equal(result1, 0L)
    expect_equal(result2, matrix(data=c(0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1,
                                            1), nrow=7, byrow=FALSE))

})

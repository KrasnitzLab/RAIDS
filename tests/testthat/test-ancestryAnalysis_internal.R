### Unit tests for ancestryAnalysis_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests specificSNVColKeep() results
#############################################################################

context("specificSNVColKeep() results")


test_that("specificSNVColKeep() must return the same object when specificSNV is NULL", {

    pRAIDS <- RAIDSparam(fileReferenceGDS=test_path("fixtures",
                                                    "1KG_Test.gds"))

    result <- RAIDS:::specificSNVColKeep(pRAIDS)

    expect_identical(result, pRAIDS)
})


test_that("specificSNVColKeep() must return the reference index of each SNV", {

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")

    gdsRef <- openfn.gds(fileRef)
    snpChr <- read.gdsn(index.gdsn(gdsRef, "snp.chromosome"))
    snpPos <- read.gdsn(index.gdsn(gdsRef, "snp.position"))
    closefn.gds(gdsRef)

    set.seed(4235)
    selected <- sample(seq_along(snpPos), 100)

    ## The last 10 positions are not in the Reference GDS
    df <- data.frame(snp.chromosome=c(snpChr[selected], rep(1, 10)),
                        snp.position=c(snpPos[selected],
                                        1987486 + 10 + seq_len(10)))

    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef, specificSNV=df)

    result <- specificSNV(RAIDS:::specificSNVColKeep(pRAIDS))

    expect_true("snvKeep" %in% colnames(result))
    expect_equal(nrow(result), 110)
    expect_equal(result$snvKeep, c(selected, rep(-1L, 10)))
    expect_true(noOpenGDS())
})


test_that("specificSNVColKeep() must return the index in the snp.KeepDefault SNVs", {

    fileRef <- local_reference_copy()

    gdsRef <- openfn.gds(fileRef, readonly=FALSE)
    snpChr <- read.gdsn(index.gdsn(gdsRef, "snp.chromosome"))
    snpPos <- read.gdsn(index.gdsn(gdsRef, "snp.position"))
    keep <- seq(1, length(snpPos), by=2)
    add.gdsn(gdsRef, "snp.KeepDefault", keep)
    closefn.gds(gdsRef)

    ## Positions 1, 5, 10 and 20 in the kept SNVs
    ## The last one is not in the kept SNVs
    selected <- c(keep[c(1, 5, 10, 20)], 2)
    df <- data.frame(snp.chromosome=snpChr[selected],
                        snp.position=snpPos[selected])

    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef, specificSNV=df)

    result <- specificSNV(RAIDS:::specificSNVColKeep(pRAIDS))

    expect_equal(result$snvKeep, c(1L, 5L, 10L, 20L, -1L))
    expect_true(noOpenGDS())
})


test_that("specificSNVColKeep() must accept chromosomes with the 'chr' prefix", {

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")

    gdsRef <- openfn.gds(fileRef)
    snpChr <- read.gdsn(index.gdsn(gdsRef, "snp.chromosome"))
    snpPos <- read.gdsn(index.gdsn(gdsRef, "snp.position"))
    closefn.gds(gdsRef)

    selected <- c(12, 3456, 7890)
    df <- data.frame(snp.chromosome=paste0("chr", snpChr[selected]),
                        snp.position=snpPos[selected])

    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef, specificSNV=df)

    result <- specificSNV(RAIDS:::specificSNVColKeep(pRAIDS))

    expect_equal(result$snvKeep, selected)
    expect_equal(result$snp.chromosome, snpChr[selected])
})


test_that("specificSNVColKeep() must remove the duplicated SNVs", {

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")

    gdsRef <- openfn.gds(fileRef)
    snpChr <- read.gdsn(index.gdsn(gdsRef, "snp.chromosome"))
    snpPos <- read.gdsn(index.gdsn(gdsRef, "snp.position"))
    closefn.gds(gdsRef)

    selected <- c(100, 200, 100, 300, 200)
    df <- data.frame(snp.chromosome=snpChr[selected],
                        snp.position=snpPos[selected])

    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef, specificSNV=df)

    result <- specificSNV(RAIDS:::specificSNVColKeep(pRAIDS))

    expect_equal(nrow(result), 3)
    expect_equal(result$snvKeep, c(100L, 200L, 300L))
})


test_that("specificSNVColKeep() must not change an existing snvKeep column", {

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")

    df <- data.frame(snp.chromosome=c(1, 1), snp.position=c(10, 20),
                        snvKeep=c(5L, 8L))

    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef, specificSNV=df)

    result <- RAIDS:::specificSNVColKeep(pRAIDS)

    expect_identical(result, pRAIDS)
})


#############################################################################
### Tests syntheticJSD() results
#############################################################################

context("syntheticJSD() results")


## Create the syntheticQ data.frame used by syntheticJSD()
## Each row is a synthetic profile with the admixture proportions
## of the synthetic profile (ag), of its origin (.p) and of the
## ground truth (.gs)
makeSyntheticQ <- function(nbRow, ag) {
    randomProp <- function() {
        m <- matrix(runif(nbRow * length(ag)), nrow=nbRow)
        return(m / rowSums(m))
    }
    res <- data.frame(randomProp(), randomProp(), randomProp())
    colnames(res) <- c(ag, paste0(ag, ".p"), paste0(ag, ".gs"))
    return(res)
}


test_that("syntheticJSD() must return the JSD between each pair of proportions", {

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")
    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef)

    ag <- c("agEUR", "agEAS", "agAMR", "agAFR", "agSAS")

    set.seed(312)
    syntheticQ <- makeSyntheticQ(nbRow=4, ag=ag)

    result <- quietly(RAIDS:::syntheticJSD(syntheticQ=syntheticQ,
                                            pRAIDS=pRAIDS))

    pairJSD <- function(a, b) {
        vapply(seq_len(nrow(a)), function(i) {
            unname(philentropy::JSD(rbind(a[i, ], b[i, ]), unit="log2"))
        }, numeric(1))
    }
    qSyn <- as.matrix(syntheticQ[, ag])
    qPos <- as.matrix(syntheticQ[, paste0(ag, ".p")])
    qGs <- as.matrix(syntheticQ[, paste0(ag, ".gs")])

    expect_equal(result[, colnames(syntheticQ)], syntheticQ)
    expect_equal(result$JSDSynVsPos, pairJSD(qSyn, qPos))
    expect_equal(result$JSDPosVsGs, pairJSD(qPos, qGs))
    expect_equal(result$JSDSynVsGs, pairJSD(qSyn, qGs))
    expect_true(noOpenGDS())
})


test_that("syntheticJSD() must return an error when there is only one synthetic profile", {

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")
    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef)

    ag <- c("agEUR", "agEAS", "agAMR", "agAFR", "agSAS")

    set.seed(312)
    syntheticQ <- makeSyntheticQ(nbRow=1, ag=ag)

    expect_error(RAIDS:::syntheticJSD(syntheticQ=syntheticQ, pRAIDS=pRAIDS),
        "Can't run syntheticJSD with only one synthetic data.", fixed=TRUE)
    expect_true(noOpenGDS())
})


#############################################################################
### Tests computePmatrix() results
#############################################################################

context("computePmatrix() results")


test_that("computePmatrix() must create the P matrix files", {

    pRAIDS <- local_admixture_pRAIDS()
    pathOut <- withr::local_tempdir()

    result <- quietly(RAIDS:::computePmatrix(pathOut=pathOut, pRAIDS=pRAIDS))

    ## The demo has 2 groups of synthetic profiles
    expect_equal(result, c(1L, 2L))
    expect_setequal(list.files(file.path(pathOut, "matP")),
                    c("ex1.P.in", "ex1.syn.1.P.in", "ex1.syn.2.P.in"))

    ## One row per pruned SNV and one column per ancestry group
    matP <- read.table(file.path(pathOut, "matP", "ex1.P.in"))
    expect_equal(dim(matP), c(61, 5))
    expect_true(all(matP >= 0 & matP <= 1))
    expect_true(noOpenGDS())
})

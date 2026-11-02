### Unit tests for processRef_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests referenceVersion() results
#############################################################################

context("referenceVersion() results")


test_that("referenceVersion() must add and append the reference version", {

    fileRef <- local_reference_copy(removeAdmixture=TRUE)
    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef)

    expect_equal(RAIDS:::referenceVersion(nameVer="v1", pRAIDS=pRAIDS), 0L)
    expect_equal(RAIDS:::referenceVersion(nameVer="v2", pRAIDS=pRAIDS), 0L)

    gdsRef <- openfn.gds(fileRef)
    withr::defer(closefn.gds(gdsRef))

    expect_equal(read.gdsn(index.gdsn(gdsRef, "reference.version")),
                    c("v1", "v2"))
})


#############################################################################
### Tests admixtureMeta(), admixtureMatrix() and admixtureAnchor() results
#############################################################################

context("admixtureMeta(), admixtureMatrix() and admixtureAnchor() results")


test_that("admixture functions must add and append two admixture runs", {

    fileRef <- local_reference_copy(removeAdmixture=TRUE)
    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef)

    ## The demo reference has 156 reference profiles
    set.seed(3)
    prop <- round(runif(156), 3)
    matQ1 <- data.frame(agA=prop, agB=1 - prop)
    matQ2 <- data.frame(agB=1 - prop, agA=prop)
    anchor1 <- data.frame(ag=rep(c("agA", "-", "agB", "agA"), length.out=156))
    anchor2 <- data.frame(ag=rev(anchor1$ag))

    for (run in list(list(id="run1", desc="first", matQ=matQ1,
                            anchor=anchor1, anchorId=1L),
                        list(id="run2", desc="second", matQ=matQ2,
                            anchor=anchor2, anchorId=2L))) {
        desc <- data.frame(id=run$id, description=run$desc,
                            anchorId=run$anchorId)
        expect_equal(RAIDS:::admixtureMeta(agGroups=colnames(run$matQ),
                            admixtureDesc=desc, pRAIDS=pRAIDS), 0L)
        expect_equal(RAIDS:::admixtureMatrix(matQ=run$matQ, admixId=run$id,
                            pRAIDS=pRAIDS), 0L)
        expect_equal(RAIDS:::admixtureAnchor(anchor=run$anchor,
                            admixId=run$id, pRAIDS=pRAIDS), 0L)
    }

    gdsRef <- openfn.gds(fileRef)
    withr::defer(closefn.gds(gdsRef))

    admixDesc <- read.gdsn(index.gdsn(gdsRef, "admixture.desc"))
    expect_equal(admixDesc$id.seq, c(1L, 2L))
    expect_equal(admixDesc$id, c("run1", "run2"))
    expect_equal(admixDesc$description, c("first", "second"))
    expect_equal(admixDesc$K, c(2L, 2L))
    expect_equal(admixDesc$anchor.id, c(1L, 2L))

    admixAG <- read.gdsn(index.gdsn(gdsRef, "admixture.ag"))
    expect_equal(admixAG$id.seq, c(1L, 1L, 2L, 2L))
    expect_equal(admixAG$column, c(1L, 2L, 1L, 2L))
    expect_equal(admixAG$ag, c("agA", "agB", "agB", "agA"))

    admixMat <- read.gdsn(index.gdsn(gdsRef, "admixture.2"))
    expect_equal(dim(admixMat), c(156, 2, 2))
    expect_equal(admixMat[, , 1], unname(as.matrix(matQ1)), tolerance=1e-6)
    expect_equal(admixMat[, , 2], unname(as.matrix(matQ2)), tolerance=1e-6)

    anchorRef <- read.gdsn(index.gdsn(gdsRef, "anchor.ref"))
    expect_equal(dim(anchorRef), c(156, 2))
    expect_equal(anchorRef[, 1], anchor1$ag)
    expect_equal(anchorRef[, 2], anchor2$ag)
})


test_that("admixtureMatrix() must return error when matQ has not one row per reference profile", {

    fileRef <- local_reference_copy(removeAdmixture=TRUE)
    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef)

    matQ <- data.frame(agA=rep(0.5, 10), agB=rep(0.5, 10))

    expect_error(RAIDS:::admixtureMatrix(matQ=matQ, admixId="run1",
                    pRAIDS=pRAIDS),
                    "Error sample.id not same size as matQ", fixed=TRUE)
    expect_true(noOpenGDS())
})


test_that("admixtureMatrix() must return error when there is no admixture.desc", {

    fileRef <- local_reference_copy(removeAdmixture=TRUE)
    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef)

    matQ <- data.frame(agA=rep(0.5, 156), agB=rep(0.5, 156))

    expect_error(RAIDS:::admixtureMatrix(matQ=matQ, admixId="run1",
                    pRAIDS=pRAIDS),
                    "Error there is no field admixture.desc", fixed=TRUE)
    expect_true(noOpenGDS())
})


test_that("admixtureMatrix() must return error when admixId is unknown or already saved", {

    fileRef <- local_reference_copy(removeAdmixture=TRUE)
    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef)

    matQ <- data.frame(agA=rep(0.5, 156), agB=rep(0.5, 156))
    RAIDS:::admixtureMeta(agGroups=colnames(matQ),
                admixtureDesc=data.frame(id="run1", description="first",
                                            anchorId=1L), pRAIDS=pRAIDS)
    RAIDS:::admixtureMatrix(matQ=matQ, admixId="run1", pRAIDS=pRAIDS)

    expect_error(RAIDS:::admixtureMatrix(matQ=matQ, admixId="unknown",
                    pRAIDS=pRAIDS),
                    "Error unknown is not in admixture.desc", fixed=TRUE)
    expect_error(RAIDS:::admixtureMatrix(matQ=matQ, admixId="run1",
                    pRAIDS=pRAIDS),
                    "Error the seq.id run1 doesn't match the posistion",
                    fixed=TRUE)
    expect_true(noOpenGDS())
})


test_that("admixtureAnchor() must return error when anchor has not one row per reference profile", {

    fileRef <- local_reference_copy(removeAdmixture=TRUE)
    pRAIDS <- RAIDSparam(fileReferenceGDS=fileRef)

    expect_error(RAIDS:::admixtureAnchor(anchor=data.frame(ag=rep("agA", 10)),
                    admixId="run1", pRAIDS=pRAIDS),
                    "Error sample.id not same size as anchor", fixed=TRUE)
    expect_true(noOpenGDS())
})

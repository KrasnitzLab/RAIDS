### Unit tests for synthetic_internal.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests validateSyntheticGeno() results
#############################################################################

context("validateSyntheticGeno() results")


test_that("validateSyntheticGeno() must return expected results when all input are valid", {

    data.dir <- test_path("fixtures")
    fileGDS <- file.path(data.dir, "1KG_Test.gds")

    gdsF <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gdsF)), envir=parent.frame())

    result1 <- RAIDS:::validateSyntheticGeno(gds=gdsF, gdsRefAnnot=gdsF,
            fileProfileGDS=fileGDS, data.id.profile="TCGA_001",
            listSampleRef=c("Sample01", "Sample02"), nbSim=1L,
            prefId="TEST", pRecomb=0.2, minProb=0.3, seqError=0.2)

    expect_identical(result1, 0L)
})


#############################################################################
### Tests prepPedSynthetic1KG() results
#############################################################################

context("prepPedSynthetic1KG() results")


test_that(paste0("prepPedSynthetic1KG() must return expected results"), {

    fileGDS <- test_path("fixtures", "1KG_Test.gds")

    gds1KG <- openfn.gds(fileGDS)
    withr::defer((gdsfmt::closefn.gds(gds1KG)), envir=parent.frame())


    fileProfileGDS <- test_path("fixtures", "GDS_Sample_TMP_201.gds")
    gdsSample <- createfn.gds(fileProfileGDS)
    withr::defer(unlink(fileProfileGDS, force=TRUE), envir=parent.frame())
    withr::defer((gdsfmt::closefn.gds(gdsSample)), envir=parent.frame())


    ## Create "study.annot" node (the node must be present)
    dataF <- data.frame(data.id=c("A101TCGA", "HG00101.Synthetic.01",
                                        "HG00101.Synthetic.02",
                                        "HG00102.Synthetic.02",
                                        "HG00109.Synthetic.02"),
                        case.id=c("A101TCGA", "HG00101", "HG00101",
                                        "HG00102", "HG00109"),
                        sample.type=c("TCGA", "Synthetic", "Synthetic",
                                        "Synthetic", "Synthetic"),
                        diagnosis=rep("C", 5),
                        source=c("TCGA", "Synthetic", "Synthetic",
                                        "Synthetic", "Synthetic"),
                        study.id=c("TCGA", "TCGA.Synthetic", "TCGA.Synthetic",
                                        "TCGA.Synthetic", "TCGA.Synthetic"))
    add.gdsn(node=gdsSample, name="study.annot", val=dataF)
    sync.gds(gdsSample)

    result1 <- RAIDS:::prepPedSynthetic1KG(gds=gds1KG,  gdsSample=gdsSample,
                                   studyID="TCGA.Synthetic", popName="superPop")

    expected1 <- data.frame(data.id=c("HG00101.Synthetic.01", "HG00101.Synthetic.02",
                            "HG00102.Synthetic.02", "HG00109.Synthetic.02"),
                            case.id=c("HG00101", "HG00101", "HG00102", "HG00109"),
                            sample.type=rep("Synthetic", 4),
                            diagnosis=rep("C", 4),
                            source=rep("Synthetic", 4),
                            study.id=rep("TCGA.Synthetic", 4),
                            superPop=rep("AFR",4))
    rownames(expected1) <- expected1$data.id

    expect_equal(result1, expected1)
})


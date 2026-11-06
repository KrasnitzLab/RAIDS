### Unit tests for ancestryAnalysis.R functions

library(RAIDS)
library(withr)
library(gdsfmt)


#############################################################################
### Tests pruningProfile2() results
#############################################################################

context("pruningProfile2() results")


## Create a RAIDSparam object with a Profile GDS file that has not been pruned
## The Profile GDS file is copied in a temporary directory
local_pruning_pRAIDS <- function(specificSNV=NULL, env=parent.frame()) {
    pathProfile <- withr::local_tempdir(.local_envir=env)

    file.copy(system.file("extdata/tests/ex1_demo.gds", package="RAIDS"),
                file.path(pathProfile, "ex1.gds"))

    ped <- data.frame(Name.ID="ex1", Case.ID="Patient_h11",
                        Diagnosis="Cancer", Sample.Type="Primary Tumor",
                        Source="Databank B", row.names="ex1",
                        stringsAsFactors=FALSE)

    studyDF <- data.frame(study.id="MYDATA", study.desc="Description",
                            study.platform="PLATFORM", stringsAsFactors=FALSE)

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")

    return(RAIDSparam(studyDF=studyDF, pedStudy=ped,
                        pathProfileGDS=pathProfile, fileReferenceGDS=fileRef,
                        specificSNV=specificSNV))
}


## Return the information about the SNVs genotyped (or not) in the profile
getGenotypedSNV <- function(genotyped=TRUE) {
    gdsProfile <- openfn.gds(system.file("extdata/tests/ex1_demo.gds",
                                            package="RAIDS"))
    geno <- read.gdsn(index.gdsn(gdsProfile, "geno.ref"))
    closefn.gds(gdsProfile)

    gdsRef <- openfn.gds(system.file("extdata/tests/ex1_good_small_1KG.gds",
                                        package="RAIDS"))
    snpInfo <- data.frame(
        snp.id=read.gdsn(index.gdsn(gdsRef, "snp.id")),
        snp.chromosome=read.gdsn(index.gdsn(gdsRef, "snp.chromosome")),
        snp.position=read.gdsn(index.gdsn(gdsRef, "snp.position")),
        stringsAsFactors=FALSE)
    closefn.gds(gdsRef)

    if (genotyped) {
        return(snpInfo[which(geno != 3), ])
    }
    return(snpInfo[which(geno == 3), ])
}

## Return the pruned SNVs saved in the Profile GDS file
getPruned <- function(pRAIDS) {
    gdsProfile <- openfn.gds(file.path(pathProfileGDS(pRAIDS), "ex1.gds"))
    pruned <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    closefn.gds(gdsProfile)
    return(pruned)
}


test_that("pruningProfile2() must keep only specific SNVs when snvKeep is already computed", {

    set.seed(1)
    snpInfo <- getGenotypedSNV()
    specific <- snpInfo[sort(sample(seq_len(nrow(snpInfo)), 50)), ]

    pRAIDS <- local_pruning_pRAIDS(specificSNV=specific[,
                                        c("snp.chromosome", "snp.position")])
    pRAIDS <- RAIDS:::specificSNVColKeep(pRAIDS)

    result <- quietly(pruningProfile2(pRAIDS))
    pruned <- getPruned(pRAIDS)

    expect_s4_class(result, "RAIDSparam")
    expect_true(length(pruned) > 0)
    expect_true(all(pruned %in% specific$snp.id))
    expect_true(noOpenGDS())
})


test_that("pruningProfile2() must compute snvKeep and keep only specific SNVs", {

    local_close_all_gds()

    set.seed(1)
    snpInfo <- getGenotypedSNV()
    specific <- snpInfo[sort(sample(seq_len(nrow(snpInfo)), 50)), ]

    pRAIDS <- local_pruning_pRAIDS(specificSNV=specific[,
                                        c("snp.chromosome", "snp.position")])

    result <- quietly(pruningProfile2(pRAIDS))
    pruned <- getPruned(pRAIDS)

    expect_true("snvKeep" %in% colnames(specificSNV(result)))
    expect_true(length(pruned) > 0)
    expect_true(all(pruned %in% specific$snp.id))
    expect_true(noOpenGDS())
})


test_that("pruningProfile2() must return an error when no specific SNV is genotyped", {

    local_close_all_gds()

    snpInfo <- getGenotypedSNV(genotyped=FALSE)

    pRAIDS <- local_pruning_pRAIDS(specificSNV=snpInfo[seq_len(5),
                                        c("snp.chromosome", "snp.position")])
    pRAIDS <- RAIDS:::specificSNVColKeep(pRAIDS)

    expect_error(quietly(pruningProfile2(pRAIDS)),
        "In pruningSample, the sample ex1 doesn't have SNPs after filters",
        fixed=TRUE)
    expect_true(noOpenGDS())
})

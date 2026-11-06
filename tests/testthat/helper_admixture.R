library(RAIDS)
library(withr)
library(gdsfmt)


## Create a RAIDSparam object for the ADMIXTURE workflow tests
## The Profile GDS file (with pruned SNVs and synthetic profiles) is copied
## in a temporary directory that will be removed automatically
## The Reference GDS file contains the 'admixture.desc', 'admixture.ag' and
## 'anchor.ref' entries ('sp0.95v0.2' anchor)
local_admixture_pRAIDS <- function(env=parent.frame()) {
    pathProfile <- withr::local_tempdir(.local_envir=env)

    file.copy(system.file("extdata/demoAncestryCall/ex1.gds", package="RAIDS"),
                file.path(pathProfile, "ex1.gds"))

    ped <- data.frame(Name.ID="ex1", Case.ID="Patient_h11",
                        Diagnosis="Cancer", Sample.Type="Primary Tumor",
                        Source="Databank B", row.names="ex1",
                        stringsAsFactors=FALSE)

    studyDFSyn <- data.frame(study.id="MYDATA.Synthetic",
                        study.desc="MYDATA synthetic data",
                        study.platform="PLATFORM", stringsAsFactors=FALSE)

    fileRef <- system.file("extdata/tests/ex1_good_small_1KG.gds",
                                package="RAIDS")

    return(RAIDSparam(studyDFSyn=studyDFSyn, pedStudy=ped,
                        pathProfileGDS=pathProfile, fileReferenceGDS=fileRef))
}


## Copy the demo Reference GDS file in a temporary file
## The admixture entries can be removed to start from a clean reference
## The file will be removed automatically
local_reference_copy <- function(removeAdmixture=FALSE, env=parent.frame()) {
    fileRef <- withr::local_tempfile(fileext=".gds", .local_envir=env)

    file.copy(system.file("extdata/tests/ex1_good_small_1KG.gds",
                            package="RAIDS"), fileRef)

    if (removeAdmixture) {
        gdsRef <- openfn.gds(fileRef, readonly=FALSE)
        for (node in c("admixture.desc", "admixture.ag", "admixture.5",
                        "anchor.ref", "reference.version")) {
            delete.gdsn(index.gdsn(gdsRef, node), force=TRUE)
        }
        closefn.gds(gdsRef)
    }

    return(fileRef)
}


## Run an expression without the printed and message output
## (SNPRelate and genio functions are verbose)
quietly <- function(expr) {
    suppressMessages(utils::capture.output(res <- expr))
    return(res)
}


## Return TRUE when no GDS file is left open
noOpenGDS <- function() {
    files <- gdsfmt::showfile.gds(closeall=FALSE, verbose=FALSE)
    return(is.null(files) || all(files$State != "open"))
}


## Close all GDS files at the end of the test so that a GDS file left open
## by a failing function does not affect the next tests
local_close_all_gds <- function(env=parent.frame()) {
    withr::defer(gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE),
                    envir=env)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RAIDS.R
\docType{data}
\name{demoPedigreeEx1}
\alias{demoPedigreeEx1}
\title{The pedigree information about a demo profile called 'ex1'.}
\format{
The \code{data.frame} containing the information about a demo
profile called 'ex1'. the \code{data.frame} has 5 columns:
\describe{
\item{Name.ID}{ a \code{character} string representing the unique
identifier of the profile.}
\item{Case.ID}{ a \code{character} string representing the unique
identifier of the case associated to the profile.}
\item{Sample.Type}{ a \code{character} string describing the type of
profile.}
\item{Diagnosis}{ a \code{character} string describing the diagnosis of the
profile.}
\item{Source}{ a \code{character} string describing the source of the
profile.}
}
}
\usage{
data(demoPedigreeEx1)
}
\value{
The \code{data.frame} containing the information about a demo
profile called 'ex1'. the \code{data.frame} has 5 columns:
\describe{
\item{Name.ID}{ a \code{character} string representing the unique
identifier of the profile.}
\item{Case.ID}{ a \code{character} string representing the unique
identifier of the case associated to the profile.}
\item{Sample.Type}{ a \code{character} string describing the type of
profile.}
\item{Diagnosis}{ a \code{character} string describing the diagnosis of the
profile.}
\item{Source}{ a \code{character} string describing the source of the
profile.}
}
}
\description{
The object is a \code{data.frame}.
}
\details{
This object can be
used to test the \code{\link{runExomeAncestry}} function.
}
\examples{


## Required library for GDS
library(SNPRelate)

## Path to the demo 1KG GDS file is located in this package
dataDir <- system.file("extdata", package="RAIDS")

#################################################################
## Load the information about the profile
#################################################################
data(demoPedigreeEx1)
head(demoPedigreeEx1)

#################################################################
## The 1KG GDS file and the 1KG SNV Annotation GDS file
## need to be located in the same directory
## Note that the 1KG GDS file used for this example is a
## simplified version and CANNOT be used for any real analysis
#################################################################
path1KG <- file.path(dataDir, "tests")

fileReferenceGDS  <- file.path(path1KG, "ex1_good_small_1KG.gds")
fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")

#################################################################
## The Sample SNP pileup files (one per sample) need
## to be located in the same directory.
#################################################################
pathGeno <- file.path(dataDir, "example", "snpPileup")

#################################################################
## The path where the Profile GDS Files (one per sample)
## will be created need to be specified.
#################################################################
pathProfileGDS <- file.path(tempdir(), "out.tmp")

pathOut <- file.path(tempdir(), "res.out")

#################################################################
## A data frame containing general information about the study
## is also required. The data frame must have
## those 3 columns: "studyID", "study.desc", "study.platform"
#################################################################
studyDF <- data.frame(study.id="MYDATA",
                        study.desc="Description",
                        study.platform="PLATFORM",
                        stringsAsFactors=FALSE)

####################################################################
## Fix seed to ensure reproducible results
####################################################################
set.seed(2043)

gds1KG <- snpgdsOpen(fileReferenceGDS)
dataRef <- select1KGPop(gds1KG, nbProfiles=2L)
closefn.gds(gds1KG)

## Required library for this example to run correctly
if (requireNamespace("GenomeInfoDb", quietly=TRUE) &&
     requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {

    ## Chromosome length information
    ## chr23 is chrX, chr24 is chrY and chrM is 25
    chrInfo <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]

    \donttest{
        runExomeAncestry(pedStudy=demoPedigreeEx1, studyDF=studyDF,
            pathProfileGDS=pathProfileGDS,
            pathGeno=pathGeno, pathOut=pathOut,
            fileReferenceGDS=fileReferenceGDS,
            fileReferenceAnnotGDS=fileAnnotGDS,
            chrInfo=chrInfo, syntheticRefDF=dataRef,
            genoSource="snp-pileup")

        unlink(pathProfileGDS, recursive=TRUE, force=TRUE)
        unlink(pathOut, recursive=TRUE, force=TRUE)
    }
}

}
\seealso{
\describe{
\item{\code{\link{runExomeAncestry}}}{ for running runs most
steps leading to the ancestry inference call on a specific exome
profile.}
}
}
\keyword{datasets}

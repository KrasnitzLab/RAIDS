% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process1KG.R
\encoding{UTF-8}
\name{generatePhase1KG2GDS}
\alias{generatePhase1KG2GDS}
\title{Adding the phase information into the Reference GDS file}
\usage{
generatePhase1KG2GDS(
  gdsReference,
  gdsReferencePhase,
  pathGeno,
  fileSNVIndex,
  verbose = FALSE
)
}
\arguments{
\item{gdsReference}{an object of class \link[gdsfmt]{gds.class} (GDS file),
an opened Reference GDS file.}

\item{gdsReferencePhase}{an object of class \link[gdsfmt]{gds.class}
(GDS file), an opened Reference Phase GDS file.}

\item{pathGeno}{a \code{character} string representing the path where
the 1K genotyping files for each sample are located. The name of the
genotyping files must correspond to
the individual identification (Individual.ID) in the pedigree file.
Default: \code{"./data/sampleGeno"}.}

\item{fileSNVIndex}{a \code{character} string representing the path and file
name of the RDS file that contains the indexes of the retained SNPs. The
file must exist. The file must be a RDS file.}

\item{verbose}{a \code{logicial} indicating if the function should
print messages when running. Default: \code{FALSE}.}
}
\value{
The function returns \code{0L} when successful.
}
\description{
The function is adding the phase information into the
Reference Phase GDS file. The phase information is extracted from a Reference
GDS file and is added into a Reference Phase GDS file. An entry called
'phase' is added to the Reference Phase GDS file.
}
\examples{

## Required package
library(gdsfmt)

## Path to the demo pedigree file is located in this package
dataDir <- system.file("extdata", package="RAIDS")

## Path where the demo genotype CSV files are located
pathGeno <- file.path(dataDir, "demoProfileGenotypes")

## The RDS file containing the pedigree information
pedigreeFile <- file.path(dataDir, "PedigreeReferenceDemo.rds")

## The RDS file containing the indexes of the retained SNPs
snpIndexFile <- file.path(dataDir, "listSNPIndexes_Demo.rds")

## The RDS file containing the filtered SNP information
filterSNVFile <- file.path(dataDir, "mapSNVSelected_Demo.rds")

## Temporary Reference GDS file containing reference information
fileReferenceGDS <- file.path(tempdir(), "1KG_TEMP_02.gds")

## Create a temporary Reference GDS file containing information from 1KG
generateGDS1KG(pathGeno=pathGeno, filePedRDS=pedigreeFile,
        fileSNVIndex=snpIndexFile, fileSNVSelected=filterSNVFile,
        fileNameGDS=fileReferenceGDS, listSamples=NULL)

## Temporary Phase GDS file that will contain the 1KG Phase information
fileRefPhaseGDS <- file.path(tempdir(), "1KG_TEMP_Phase_02.gds")

## Create Reference Phase GDS file
gdsPhase <- createfn.gds(fileRefPhaseGDS)

## Open Reference GDS file
gdsRef <- openfn.gds(fileReferenceGDS)

## Fill temporary Reference Phase GDS file
if (FALSE) {
    generatePhase1KG2GDS(gdsReference=gdsRef,
            gdsReferencePhase=gdsPhase,
            pathGeno=pathGeno, fileSNVIndex=snpIndexFile,
            verbose=FALSE)
}

## Close Reference Phase information file
closefn.gds(gdsPhase)

## Close Reference information file
closefn.gds(gdsRef)

## Remove temporary files
unlink(fileReferenceGDS, force=TRUE)
unlink(fileRefPhaseGDS, force=TRUE)

}
\author{
Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
}

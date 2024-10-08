% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process1KG.R
\encoding{UTF-8}
\name{generateMapSnvSel}
\alias{generateMapSnvSel}
\title{Generate the filter SNP information file in RDS format}
\usage{
generateMapSnvSel(cutOff = 0.01, fileSNV, fileSNPsRDS, fileFREQ)
}
\arguments{
\item{cutOff}{a single \code{numeric} value, the cut-off
for the frequency in at least one super population. Default: \code{0.01}.}

\item{fileSNV}{a \code{character} string representing the path and
file name of the bulk SNP information file from Reference. The file must
be in text format. The file must exist.}

\item{fileSNPsRDS}{a \code{character} string representing the path and
file name of the RDS file that will contain the indexes of the retained
SNPs. The file extension must be '.rds'.}

\item{fileFREQ}{a \code{character} string representing the path and
file name of the RDS file that will contain the filtered SNP information.
The file extension must be '.rds'.}
}
\value{
The integer \code{0} when successful.
}
\description{
The function applies a cut-off filter to the SNP information
file to retain only the SNP that have a frequency superior or equal to the
specified cut-off in at least one super population. The information about
the retained SNPs is saved in a RDS format file. A RDS file containing the
indexes of the retained SNP is also created.
}
\details{
The filtered SNP information RDS file (parameter \code{fileFREQ}), contains
a \code{data.frame} with those columns:
\describe{
\item{CHROM}{ a \code{character} string representing the chromosome where
the SNV is located.}
\item{POS}{ a \code{character} string representing the SNV position on the
chromosome.}
\item{REF}{ a \code{character} string representing the reference DNA base
for the SNV.}
\item{ALT}{ a \code{character} string representing the alternative DNA base
for the SNV.}\
\item{EAS_AF}{ a \code{character} string representing the allele frequency
of the EAS super population.}
\item{AFR_AF}{ a \code{character} string representing the allele frequency
of the AFR super population.}
\item{AMR_AF}{ a \code{character} string representing the allele frequency
of the AMR super population.}
\item{SAS_AF}{ a \code{character} string representing the allele frequency
of the SAS super population.}
}
}
\examples{

## Path to the demo pedigree file is located in this package
dataDir <- system.file("extdata", package="RAIDS")

## Demo SNV information file used as input
snvFile <- file.path(dataDir, "matFreqSNV_Demo.txt.bz2")

## Temporary output files
## The first file contains the indexes of the retained SNPs
## The second file contains the filtered SNP information
snpIndexFile <- file.path(tempdir(), "listSNP_TEMP.rds")
filterSNVFile <- file.path(tempdir(), "mapSNVSel_TEMP.rds")

## Create a data.frame containing the information of the retained
## samples (samples with existing genotyping files)
generateMapSnvSel(cutOff=0.01, fileSNV=snvFile,
        fileSNPsRDS=snpIndexFile, fileFREQ=filterSNVFile)

## Remove temporary files
unlink(snpIndexFile, force=TRUE)
unlink(filterSNVFile, force=TRUE)

}
\author{
Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
}

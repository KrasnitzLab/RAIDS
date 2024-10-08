% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allelicFraction_internal.R
\encoding{UTF-8}
\name{tableBlockAF}
\alias{tableBlockAF}
\title{Compile the information about the SNVs
for each block}
\usage{
tableBlockAF(snpPos)
}
\arguments{
\item{snpPos}{a \code{data.frame} with lower allelic fraction (lap) for
the SNVs with coverage > \code{minCov}, for a specific chromosome.}
}
\value{
a \code{data.frame} containing only heterozygote
SNV information. The
\code{data.frame} contain those columns:
\describe{
\item{block}{ a single \code{integer} representing the unique identifier
of the block.}
\item{aRF}{ a single \code{numeric} representing the final allelic
fraction; not computed yet, \code{-1} value assigned to all entries.}
\item{aFraction}{a single \code{integer} representing the possible allelic
fraction in absence of loss of heterozygosity (LOH).}
\item{lR}{ a single \code{integer} representing the coverage for
the alternative allele.}
\item{nPhase}{ a single \code{integer} representing the number of SNV
phases.}
\item{sumAlleleLow}{ a single \code{integer} representing the sum of the
alleles with the less coverage.}
\item{sumAlleleHigh}{ a single \code{integer} representing the sum of
the alleles with more coverage.}
\item{lH}{ a single \code{numeric} for the homozygotes log10 of the product
frequencies of the allele not found in the profile (not a probability).}
\item{lM}{ a single \code{numeric} log10 product frequency allele
in population.}
\item{lRhomo}{a single \code{numeric} representing the score
\code{lH} - \code{lM}.}
\item{nbHomo}{ a single \code{integer} representing the number of
homozygote SNVs per block.}
\item{nbKeep}{ a single \code{integer} representing the number of
SNVs retained per block.}
\item{nbHetero}{ a single \code{integer} representing the number of
heterozygote SNVs per block.}
}
}
\description{
The function evaluates a score
about loss of heterozygosity and allelic fraction for each block. It
generates specific information about the variants in the block, like the
number of homozygotes or heterozygotes.
In the case of RNA-seq, the blocks are genes.
}
\examples{

## Loading demo dataset containing SNV information
data(snpPositionDemo)

## Retain SNVs on chromosome 1
subset <- snpPositionDemo[which(snpPositionDemo$snp.chr == 1),]

##Compile the information about the SNVs for each block
result <- RAIDS:::tableBlockAF(subset)
head(result)

}
\author{
Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
}
\keyword{internal}

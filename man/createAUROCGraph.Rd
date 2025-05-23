% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\encoding{UTF-8}
\name{createAUROCGraph}
\alias{createAUROCGraph}
\title{Generate accuracy graph}
\usage{
createAUROCGraph(
  dfAUROC,
  title = "",
  selectD = c(3, 7, 11),
  selectColor = c("#5e688a", "#cd5700", "#CC79A7")
)
}
\arguments{
\item{dfAUROC}{a \code{data.frame} corresponding to res$paraSample$dfAUROC
where res is the result of inferAncestry() or inferAncestryGeneAware()
functions.}

\item{title}{a \code{character} string representing the title of the graph.
Default: \code{""}.}

\item{selectD}{a \code{array} of \code{integer} representing the selected
PCA dimensions to plot. The length of the \code{array} cannot be more than
5 entries. The dimensions must tested by RAIDS (i.e. be present in the
RDS file). Default: \code{c(3,7,11)}.}

\item{selectColor}{a \code{array} of \code{character} strings representing
the selected colors for the associated PCA dimensions to plot. The length
of the \code{array} must correspond to the length of the \code{selectD}
parameter. In addition, the length of the \code{array} cannot be more than
5 entries.
Default: \code{c("#5e688a", "#cd5700", "#CC79A7")}.}
}
\value{
a \code{ggplot} object containing the graphic representation of the
accuracy for different values of PCA dimensions and K-neighbors through
all tested ancestries.
}
\description{
This function extracts the required information from an
output generated by RAIDS to create a graphic representation of the
accuracy for different values of PCA dimensions and K-neighbors through
all tested ancestries.
}
\examples{

## Required library
library(ggplot2)

## Path to RDS file with ancestry information generated by RAIDS (demo file)
dataDir <- system.file("extdata", package="RAIDS")
fileRDS <- file.path(dataDir, "TEST_01.infoCall.RDS")
info <- readRDS(fileRDS)
dfAUROC <- info$paraSample$dfAUROC

## Some of the column names must be updated to fit new standards
colnames(dfAUROC) <- c("D", "K", "Call", "L", "AUROC", "H")

## Create accuracy graph
accuracyGraph <- createAUROCGraph(dfAUROC=dfAUROC, title="Test 01",
    selectD=c(3, 6, 9, 12, 15),
    selectColor=c("steelblue", "darkorange", "violet", "pink", "gray40"))

accuracyGraph

}
\author{
Astrid Deschênes and Pascal Belleau
}

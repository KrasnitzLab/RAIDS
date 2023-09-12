
######################################################################
## How to create the PopulationReferenceSNVAnnotationDemo.gds file
## This file is a small population reference SNV Annotation GDS file
## with limited number of samples and SNVs
######################################################################

## Required librairies
library(SNPRelate)
library(gdsfmt)

## Create a GDS file
gdsRefNew <- createfn.gds("PopulationReferenceSNVAnnotationDemo.gds")

## The entry 'phase' contain the phase of the SNVs in the
## Population Annotation GDS file
## 0 means the first allele is a reference; 1 means the first allele is
## the alternative and 3 means unknown
## The SNVs (rows) and samples (columns) in phase are in the same order as
## in the Population Annotation GDS file.
phase <- matrix(data=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                        0, 0, 0, 1, 1, 0, 0, 0, 1, 1,
                        0, 0, 0, 1, 1, 0, 0, 0, 0, 1,
                        1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                        0, 1, 0, 1, 1, 0, 1, 1, 1, 1,
                        0, 1, 0, 1, 1, 0, 1, 1, 1, 1,
                        0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                       0, 1, 0, 1, 1, 0, 1, 1, 1, 1), ncol=10, byrow=TRUE)
add.gdsn(node=gdsRefNew, name="phase", val=phase, storage="bit2", check=TRUE)

## The entry 'blockAnnot' contain the information for each group of blocks
## that are present in the 'block' entry.
blockAnnot <- data.frame(block.id=c("EAS.0.05.500k", "EUR.0.05.500k",
                    "AFR.0.05.500k", "AMR.0.05.500k", "SAS.0.05.500k"),
                block.desc=c(
                    "EAS populationblock base on SNP 0.05 and windows 500k",
                    "EUR populationblock base on SNP 0.05 and windows 500k",
                    "AFR populationblock base on SNP 0.05 and windows 500k",
                    "AMR populationblock base on SNP 0.05 and windows 500k",
                    "SAS populationblock base on SNP 0.05 and windows 500k"),
                stringsAsFactors=FALSE)
add.gdsn(node=gdsRefNew, name="block.annot", val=blockAnnot, check=TRUE)

## The entry 'block' contain the block information for the SNVs in the
##  Population Annotation GDS file.
## The SNVs (rows) are in the same order as in
##  the Population Annotation GDS file.
## The block groups (columns) are in the same order as in
##  the 'block.annot' entry.
block <- matrix(data=c(-1, -1, -1, -1, -1,
                       -2, -2,  1, -2, -2,
                       -2,  1,  1,  1, -2,
                       -2,  1,  1,  1, -2,
                       -2, -3, -2, -3, -2,
                        1,  2,  2,  2,  1,
                        1,  2,  2,  2,  1,
                       -3, -4, -3, -4, -3,
                        2, -4,  3, -4, -3,
                        2, -4,  3, -4, -3), ncol=5, byrow=TRUE)
add.gdsn(node=gdsRefNew, name="block", val=block, storage="int32", check=TRUE)

closefn.gds(gdsRefNew)

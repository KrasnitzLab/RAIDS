
############################################################
## How to create the PopulationReferenceDemo.gds file
## This file is a small population reference GDS file
## with limited number of samples and SNVs
############################################################

## Required librairies
library(SNPRelate)
library(gdsfmt)

## Create a GDS file
gdsRefNew <- createfn.gds("PopulationReferenceDemo.gds")

## The entry 'sample.id' contain the unique identifiers of 10 samples
## that constitute the reference dataset
sample.id <- c("HG00100", "HG00101", "HG00102", "HG00103", "HG00104",
               "HG00105", "HG00106", "HG00107", "HG00108", "HG00109")
add.gdsn(node=gdsRefNew, name="sample.id", val=sample.id,
         storage="string", check=TRUE)

## A data frame containing the information about the 10 samples
## (in the same order than in the 'sample.id') is created and added to
## the 'sample.annot' entry
## The data frame must contain those columns:
##     'sex': '1'=male, '2'=female
##     'pop.group': acronym for the population (ex: GBR, CDX, MSL, ASW, etc..)
##     'superPop': acronym for the super-population (ex: AFR, EUR, etc...)
##     'batch': number identifying the batch of provenance
sampleInformation <- data.frame(sex=c("1", "1", "1", "1", "1",
        "2", "2", "2", "2", "2"),
    pop.group=c("GBR", "GIH", "CDX", "GBR", "LWK",
        "LWK", "LWK", "GBR", "GIH", "PEL"),
    superPop=c("EUR", "SAS", "EAS", "EUR", "AFR", "AFR", "AFR",
        "EUR", "SAS", "AMR"), batch=rep(0, 10), stringsAsFactors=FALSE)
add.gdsn(node=gdsRefNew, name="sample.annot", val=sampleInformation,
         check=TRUE)

## The identifier of each SNV is added in the 'snp.id' entry
snvID <- c("s1", "s2", "s3", "s4", "s5", "s6", "s7")
add.gdsn(node=gdsRefNew, name="snp.id", val=snvID,
         check=TRUE)

## The chromosome of each SNV is added to the 'snp.chromosome' entry
## The order of the SNVs is the same than in the 'snp.id' entry
snvChrom <- c(rep(1, 7))
add.gdsn(node=gdsRefNew, name="snp.chromosome", val=snvChrom, storage="uint16",
         check=TRUE)

## The position on the chromosome of each SNV is added to
## the 'snp.position' entry
## The order of the SNVs is the same than in the 'snp.id' entry
snvPos <- c(16102, 51478, 51897, 51927, 54489, 54707, 54715)
add.gdsn(node=gdsRefNew, name="snp.position", val=snvPos, storage="int32",
         check=TRUE)

## The allele information of each SNV is added to the 'snp.allele' entry
## The order of the SNVs is the same than in the 'snp.allele' entry
snvAllele <- c("T/G", "T/A", "C/A", "G/A", "G/A", "G/C", "C/T")
add.gdsn(node=gdsRefNew, name="snp.allele", val=snvAllele, storage="string",
         check=TRUE)

## The allele frequency in the general population (between 0 and 1) of each
## SNV is added to the 'snp.AF' entry
## The order of the SNVs is the same than in the 'snp.id' entry
snvAF <- c(0.02, 0.11, 0.08, 0.07, 0.10, 0.23, 0.21)
add.gdsn(node=gdsRefNew, name="snp.AF", val=snvAF, storage="packedreal24",
         check=TRUE)

## The allele frequency in the East Asian population (between 0 and 1) of each
## SNV is added to the 'snp.EAS_AF' entry
## The order of the SNVs is the same than in the 'snp.id' entry
snvAF <- c(0.00, 0.00, 0.05, 0.01, 0.00, 0.08, 0.07)
add.gdsn(node=gdsRefNew, name="snp.EAS_AF", val=snvAF, storage="packedreal24",
         check=TRUE)

## The allele frequency in the European population (between 0 and 1) of each
## SNV is added to the 'snp.EUR_AF' entry
## The order of the SNVs is the same than in the 'snp.id' entry
snvAF <- c(0.04, 0.20, 0.12, 0.14, 0.18, 0.38, 0.34)
add.gdsn(node=gdsRefNew, name="snp.EUR_AF", val=snvAF, storage="packedreal24",
         check=TRUE)

## The allele frequency in the African population (between 0 and 1) of each
## SNV is added to the 'snp.AFR_AF' entry
## The order of the SNVs is the same than in the 'snp.id' entry
snvAF <- c(0.03, 0.02, 0.05, 0.06, 0.02, 0.18, 0.16)
add.gdsn(node=gdsRefNew, name="snp.AFR_AF", val=snvAF, storage="packedreal24",
         check=TRUE)

## The allele frequency in the American population (between 0 and 1) of each
## SNV is added to the 'snp.AMR_AF' entry
## The order of the SNVs is the same than in the 'snp.id' entry
snvAF <- c(0.03, 0.12, 0.06, 0.07, 0.10, 0.25, 0.24)
add.gdsn(node=gdsRefNew, name="snp.AMR_AF", val=snvAF, storage="packedreal24",
         check=TRUE)

## The allele frequency in the South Asian population (between 0 and 1) of each
## SNV is added to the 'snp.SAS_AF' entry
## The order of the SNVs is the same than in the 'snp.id' entry
snvAF <- c(0.02, 0.22, 0.10, 0.09, 0.21, 0.28, 0.27)
add.gdsn(node=gdsRefNew, name="snp.SAS_AF", val=snvAF, storage="packedreal24",
         check=TRUE)

## The genotype of each SNV for each sample is added to the 'genotype' entry
## The genotype correspond to the number of A alleles
## The rows represent the SNVs is the same order than in 'snp.id' entry
## The columns represent the samples is the same order than in 'sample.id' entry
genotypeInfo <- matrix(data=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 1, 1, 1, 0, 1, 1, 0,
                              1, 1, 2, 0, 0, 0, 1, 1, 0, 0,
                              1, 1, 2, 0, 0, 0, 1, 1, 0, 0,
                              0, 0, 0, 1, 1, 1, 0, 1, 1, 0,
                              1, 1, 2, 0, 1, 1, 1, 1, 0, 0,
                              1, 1, 1, 0, 1, 1, 1, 1, 0, 0), ncol=10, byrow=TRUE)
add.gdsn(node=gdsRefNew, name="genotype", val=genotypeInfo,
         storage="bit2", check=TRUE)

## The entry 'sample.ref' is filled with 0 or 1 indicating that 8 samples out
## of 10 samples are retained to be used as reference
## The order of the samples is the same than in the 'sample.id' entry
add.gdsn(node=gdsRefNew, name="sample.ref", val=c(rep(1L, 7), 0, 0, 1),
         storage="bit1", check=TRUE)

closefn.gds(gdsRefNew)

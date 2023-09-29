###############################################################################
## README file
##
## This file describes the objects present in the inst/extdata directory.
###############################################################################


###############################################################################
## extdata/demoGenoChr.tar
###############################################################################

The files present in the extdata/demoGenoChr.tar files are the first
100 genotypes of 3 reference profiles (NA12003, NA12004 and NA12005) for
each chromosome. The genotypes are directly extracted from the 1KG VCF files
that is available at:
 http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV


###############################################################################
## tests/ex1_good_small_1KG.gds
## tests/ex1_good_small_1KG_Annot.gds
###############################################################################

Both files are based on subset of the genotypes from 1000 Genomes described in:
 https://wellcomeopenresearch.org/articles/4-50/

The 500 SNVs from each chromosome habe been retained (11000 for the
22 autosomal chromosomes). For each of the 26 sub-continental populations
in this version of 1000 Genomes, 6 profiles have been randomly selected
(total of 156 profiles) to generate those reference GDS files.


###############################################################################
## example/snpPileup/ex1.vcf.gz
## example/snpPileup/ex1.generic.txt.gz
## example/snpPileup/ex1.txt.gz
## example/snpPileupRNA/ex1.generic.txt.gz
## example/snpPileupRNA/ex1.txt.gz
## demoAncestryCall/ex1.gds
## demoKNNSynthetic/ex1.gds
###############################################################################

This demonstration profile contains genotype from a synthetic profile
generated using a reference profile from 1000 Genome randomly selected
with an AFR ancestry (this profile is not present in the reference GDS file
tests/ex1_good_small_1KG.gds) and one PDAC cancer exome profile
from Tiriac et al 2018 (this exome is restricted access in dbGAP; it was
only use to extract the total read deep at the SNV positions and to estimate
allelic fraction of the regions). Only 500 SNVs per chromosome have been kept.
The pipeline to generate the synthetic data are describe in Belleau et al 2023.


###############################################################################
## demoPorfileGenotypes/HG00*.csv.bz2
###############################################################################

Each file contains the first 10 SNVs from chromosome 1 of profiles
present in 1000 Genomes as described in:
 https://wellcomeopenresearch.org/articles/4-50/


###############################################################################
## matFreqSNV_Demo.txt.bz2
###############################################################################

This is a text file with the frequency for the 5 super populations of the
10 first SNVs from chromosome 1 in 1000 Genomes as described in:
 https://wellcomeopenresearch.org/articles/4-50/


###############################################################################
## listSNPIndexes_Demo.rds
## mapSNVSelected_Demo.rds
###############################################################################

Those 2 files are created using the information from the
matFreqSNV_Demo.txt.bz2 file. The listSNPIndexes_Demo.rds is an index of the
SNVs that have a frequency for at least one of the 1000 Genomes super-population
bigger than 1%. The mapSNVSelected_Demo.rds is a subset of
matFreqSNV_Demo.txt.bz2 containing the selected SNVs present in
listSNPIndexes_Demo.rds.


###############################################################################
## unrelatedPatientsInfo_Demo.rds
###############################################################################

The file contains the list of unrelated profiles from 1000 Genomes.
The values does not represent the real relation between profiles; the
values have been modified for demonstration purpose.
The profiles from 1000 Genomes are described in:
 https://wellcomeopenresearch.org/articles/4-50/



# Generate the 1000 genomes GDS file #

Preparation of the 1000 genomes (1KG) GDS file for RAIDS

## Authors ##

[Pascal Belleau](http://ca.linkedin.com/in/pascalbelleau "Pascal Belleau"),
[Astrid Desch&ecirc;nes](http://ca.linkedin.com/in/astriddeschenes "Astrid Desch&ecirc;nes") and
[Alexander Krasnitz](https://www.cshl.edu/research/faculty-staff/alexander-krasnitz/ "Alexander Krasnitz")

## Download needs ##

you need to download the vcf files from GRCh38
The files used are (add ref):

    ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr*.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz*


You need a pedigree file for the description of the samples

    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped

## Prepare a rds file in R from the pedigree file ##

In R:

    library(RAIDS)
    ## pedigree file from 1KG with the path
    ## If the file is in the current directory
    filePED1KGOri <- "20130606_g1k.ped"
    ped <- prepPed1KG(filePED1KGOri)
    ## save the RDS file here we 
    ##save the file in data/metadata/ped1KG.rds
    filePED1KG <- file.path("data", "metadata","ped1KG.rds")
    saveRDS(ped,  filePED1KG)
    


## Prepare intermediate files based on 1KG vcf for GDS

We generate some intermediate files to prepare the data to be import in R GDS file (Bioconductor gdsfmt)

Generate a file with the SNP position and the frequence in each super population.

PATHVCF is the path to the vcf files from 1KG
The script generates the matFreq.${chr}.txt
We will change extract1000gFreq.py soon.

    
    for i in `ls PATHVCF/*shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz`
    do
    chr=$(echo $FILECUR|perl -n -e '/ALL\.(chr[^\.]+)/;print $1')
    python PATH2SCRIPT/extract1000gFreq.py ${FILECUR} matFreq.${chr}
    done
    for i in `seq 1 22`
    do
    cat matFreq.chr${i}.txt >matFreqSNV.txt
    bzip2 matFreqSNV.txt
    done


Generate a list of SNP with frequency higher then a cutoff (here 0.01) for at least 1 super population.

The file is mapSNVSel.rds

In R:

    library(RAIDS)

    fileSNV.v <- file.path(PATHGENO, "matFreqSNV.txt.bz2")
    fileLSNP.v <- file.path(PATHSEL, "listSNP.rds")
    fileFREQ.v <- file.path(PATHSEL, "mapSNVSel.rds")
    generateMapSnvSel(cutOff = 0.01,
                      fileSNV = fileSNV.v,
                      fileLSNP = fileLSNP.v,
                      fileFREQ = fileFREQ.v)



## Split the VCF by sample

NOTE: Should be implemented

## Generate the base GDS file with 1KG

    library(RAIDS)
    
    filePED1KG <- file.path(PATHMETA,"ped1KG.rds")
    fileSNV.v <- file.path(PATHGENO, "matFreqSNV.txt.bz2")
    fileLSNP.v <- file.path(PATHSEL, "listSNP.rds")
    fileFREQ.v <- file.path(PATHSEL, "mapSNVSel.rds")
    
    fileNameGDS <- "matGeno1000g.gds"
    fileGDS <- file.path(PATHGDS, fileNameGDS)
    fileIBD <- file.path(PATHMETA,"ibd.All.0.05.rds")
    filePart <- file.path(PATHMETA,"part.All.0.05.rds")
    
    # Generate the base structure of the gds file
    generateGDS1KG(PATHGENO = file.path("data", "sampleGeno"),
                   fileNamePED = filePED1KG,
                   fileSNPSel = fileFREQ.v,
                   fileListSNP = fileLSNP.v,
                   fileNameGDS = fileGDS)
    
    # Identify the individual related
    gds <- snpgdsOpen(fileGDS)
    identifyRelative(gds,
                     maf = 0.05,
                     thresh = 2^(-11/2),
                     fileIBD = fileIBD,
                     filePart = filePart)
    closefn.gds(gds)
    
    # create the ref sample list which correspond
    # to the list of sample unrelated
    
    addRef2GDS1KG(fileGDS,
                  filePart)

    
## The base GDS for 1KG
    
        File: fileGDS (15.4G)
        +    [  ] *
        # The id for the genotype
        |--+ sample.id   { Str8 2548, 19.9K }
        # data.frame for the annotation of sample
        # for all genotype good for the reference 
        # samples but sn other data.frame will be 
        # define to describe the samples 
        # in the studies
        |--+ sample.annot   [ data.frame ] *
        |  |--+ sex   { Str8 2548, 5.0K }
        # population 
        |  |--+ pop.group   { Str8 2548, 10.0K }
        # The 5 super populations in 
        |  |--+ superPop   { Str8 2548, 10.0K }
        |  \--+ batch   { Float64 2548, 19.9K }
        |--+ snp.id   { Str8 24542710, 223.5M }
        |--+ snp.chromosome   { UInt16 24542710, 46.8M }
        |--+ snp.position   { Int32 24542710, 93.6M }
        |--+ snp.allele   { Str8 24542710, 93.6M }
        |--+ snp.AF   { PackedReal24 24542710, 70.2M }
        |--+ snp.EAS_AF   { PackedReal24 24542710, 70.2M }
        |--+ snp.EUR_AF   { PackedReal24 24542710, 70.2M }
        |--+ snp.AFR_AF   { PackedReal24 24542710, 70.2M }
        |--+ snp.AMR_AF   { PackedReal24 24542710, 70.2M }
        |--+ snp.SAS_AF   { PackedReal24 24542710, 70.2M }
        |--+ genotype   { Bit2 24542710x2548, 14.6G }
        \--+ sample.ref   { Bit1 2548, 319B }
    


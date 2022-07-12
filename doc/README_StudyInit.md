# Initialize a study

You need a mapped BAM with the same genome than the 
reference here (hg38 for 1KG)
for each sample.

## Extract the read at each in the GDS

We used here snp-pileup from Facet

First you a vcf file with the SNV you want to keep

## Create a VCF file with the SNV you want tokeep
In R
    library(aicsPaper)
    
    gds <- snpgdsOpen(fileGDS1kg)
    snvListVCF(gds, fileOUT, offset = 1)
    closefn.gds(gds)
    
You should compress and indexing the vcf file
You need to install [HTSlib](http://www.htslib.org/download/)

in a terminal

     bgzip fileOUT
     tabix -p vcf fileOUT.gz
     


You need a ped file with the column:

"sample.id"* Id of the sample but the sample can be genotyped more than one
"sex"*
"pop.group"* can be self declare or something else but must be there
"superPop"* can be self declare or something else but must be there
"Sample.Type"* Ex Primary Tumor or Blood Derived Normal 
"Diagnosis"* C or N (cancer normal) this is for the sample not the patient
"Source"* The tissue ex Ovary 
"Recurrent"  
"Case.ID" Patient identifiant     
"Name.ID"* This the unique id of the table
"batch" Number

You need the extra code in the package Facet to extract the SNV
[snp-pileup](https://github.com/mskcc/facets/tree/master/inst/extcode)

The output from snp-pileup should be Name.ID.txt 
(snp-pileup add .gz to the filename) 


    snp-pileup -g -d5000 -q15 -Q20  -r0 VCFGenerateAbove snvSel0.01.vcf.gz OUTPUT.txt FILEBAM.bam
    
    

TODO change the path to something generic

You add the genotype call from the SNP-pileup to
the gds 1KG and create a GDS file for each sample. 
with the covereage information

in R:

    library(aicsPaper)
    
    PATHGENO <- file.path("data", "snpCancer")
    PATHMETA <- file.path("data", "metadata")
    PATHGDS <- file.path("data", "genoGDS1000gAF.All.TCGA_OV.WXS")
    
    fileNamePED <- "pedTCGA_OV_WXS_C.rds"
    fileNameGDS <- "matGeno1000g.gds"
    
    filePED <- file.path(PATHMETA, fileNamePED)
    fileNameGDS <- file.path(PATHGDS, fileNameGDS)
    
    ped <- readRDS(filePED)
    
    studyDF <- data.frame(study.id = "TCGA-OV.WXS.C",
                              study.desc = "Ovarian example",
                              study.platform = "WXS",
                              stringsAsFactors = FALSE)
    
    listSamples <- ped[, "Name.ID"]
    appendStudy2GDS1KG(PATHGENO = PATHGENO,
                       fileNamePED = filePED,
                       fileNameGDS = fileNameGDS,
                       listSamples = listSamples,
                       studyDesc = studyDF,
                       batch = 1)





The gds per sample at this step
    +    [  ]
    |--+ Ref.count   { SparseInt16 24516859x1, 7.2M }
    |--+ Alt.count   { SparseInt16 24516859x1, 1.0M }
    |--+ Total.count   { SparseInt16 24516859x1, 7.6M }
    |--+ sampleStudy   { Str8 1, 29B }

It the coverage it define at each position in the GDS 1KG

We select a subset of SNV pruned SNV

I include a R script with the possibility to run it in parallele
on mutliple instance.

    Rscript ${PATHRSCRIPT}/runPruningStudy1KG.R ${PATHPKG} ${PATHGDS} matGeno1000g.gds $SGE_TASK_ID $PATHOUT $FILEPREF $STUDYNAME

or in R

    library(aicsPaper)
    
    # fileNameGDS with the sample and 1KG genotype
    # listSamples list of sample.id from the study
    # PATHSAMPLEGDS is the path where the gds specific 
    # to the sample is created
    
    gds <- snpgdsOpen(fileNameGDS)
    
    for(i in seq_len(length(listSamples))){
        print(system.time(pruned <- pruningSample(gds=gds,
                                            method="corr",
                                            sampleCurrent = listSamples[i],
                                            listSNP = NULL,
                                            slide.max.bp.v = 5e5,
                                            ld.threshold.v=sqrt(0.1),
                                            np = 1,
                                            verbose.v=FALSE,
                                            chr = NULL,
                                            minAF.SuperPop = NULL,
                                            keepGDSpruned = TRUE,
                                            PATHSAMPLEGDS = PATHSAMPLEGDS,
                                            keepFile = FALSE)))
    }
    
    Sys.time()
    closefn.gds(gds)

The function pruningSample add the snp.id of the snp selected 
in the pruning process to the gds Sample

    +    [  ]
    |--+ Ref.count   { SparseInt16 24516859x1, 7.2M }
    |--+ Alt.count   { SparseInt16 24516859x1, 1.0M }
    |--+ Total.count   { SparseInt16 24516859x1, 7.6M }
    |--+ sampleStudy   { Str8 1, 29B }
    |--+ pruned.study   { Str8 237908, 2.2M }


TODO
Add genotype only for the snp selected in pruning to gdsSample

    library(aicsPaper)
    
    # fileNameGDS1KG the file to the gds1KG
    # PATHSAMPLEGDS the path where the files gdsSample 
    
    gds <- snpgdsOpen(fileNameGDS1KG)
    listGDSSample <- dir(PATHSAMPLEGDS, pattern = ".+.gds")
    
    for(gdsSampleFile in listGDSSample){
        print(system.time(add1KG2SampleGDS(gds, file.path(PATHSAMPLEGDS, gdsSampleFile) )))
    }


see GDSSample

    +    [  ]
    # Ref.count, Alt.count, Total.count are 
    # all the snv from 1KG
    |--+ Ref.count   { SparseInt16 24516859x1, 7.2M }
    |--+ Alt.count   { SparseInt16 24516859x1, 1.0M }
    |--+ Total.count   { SparseInt16 24516859x1, 7.6M }
    # Genotype from pruning
    |--+ sampleStudy   { Str8 1, 29B }    
    |--+ pruned.study   { Str8 237908, 2.2M }
    |--+ sample.id   { Str8 2470, 19.3K }
    |--+ snp.id   { Str8 237908, 2.2M }
    |--+ snp.chromosome   { Int32 237908, 929.3K }
    |--+ snp.position   { Int32 237908, 929.3K }
    |--+ snp.index   { Int32 237908, 929.3K }
    |--+ genotype   { Bit2 237908x2470, 140.1M }
    |--+ SamplePos   { Float64 1, 8B }
    |--+ lap    { PackedReal8 237908, 232.3K }
    \--+ segment   { UInt32 237908, 929.3K } # not sure I keep it


Section annotation study
    |--+ study.list   [ data.frame ] *
    |  |--+ study.id   { Str8 1, 9B }
    |  |--+ study.desc   { Str8 1, 13B }
    |  \--+ study.platform   { Str8 1, 4B }
    \--+ study.annot   [ data.frame ] *
       |--+ data.id   { Str8 188, 1.1K }
       |--+ case.id   { Str8 188, 1.1K }
       |--+ sample.type   { Str8 188, 376B }
       |--+ diagnosis   { Str8 188, 376B }
       |--+ source   { Str8 188, 1.7K }
       \--+ study.id   { Str8 188, 1.7K }
Not yet the good struct for option 1 and 2

option 1 all the hetero from ex varscan

    |--+ Ref.count.o   { SparseInt16 24516859x1, 7.2M }
    |--+ Alt.count.o   { SparseInt16 24516859x1, 1.0M }
    |--+ Total.count.o   { SparseInt16 24516859x1, 7.6M }
    |--+ snp.id.o   { Str8 237908, 2.2M }
    |--+ snp.allele.o   { Str8 237908, 2.2M }
    |--+ snp.chromosome.o   { Int32 237908, 929.3K }
    |--+ snp.position.o   { Int32 237908, 929.3K }
    |--+ normal.geno   { Bit2 237908x2470, 140.1M }
    |--+ cancer.geno   { Bit2 237908x2470, 140.1M }
    
option 2 all the hetero in normal from ex varscan we extract 
the coverage in cancer at the position where the hetero
are (no overlap with anp.id)

    |--+ Ref.count.N   { SparseInt16 24516859x1, 7.2M }
    |--+ Alt.count.N   { SparseInt16 24516859x1, 1.0M }
    |--+ Total.count.N   { SparseInt16 24516859x1, 7.6M }
    |--+ snp.id.N   { Str8 237908, 2.2M }
    |--+ snp.chromosome.N   { Int32 237908, 929.3K }
    |--+ snp.position.N   { Int32 237908, 929.3K }

Block

    File: /mnt/wigclust5/data/unsafe/belleau/process1000G/samples1000gUnrelated/data/genoGDS1KG.2022.04.18/testBlock.gds (69.2M)
    +    [  ]
    |--+ block.annot   [ data.frame ] *
    |  |--+ block.id   { Str8 2, 28B }
    |  \--+ block.desc   { Str8 2, 108B }
    \--+ block   { Int32 9076010x2 LZ4_ra(3.75%), 2.6M }


We create another GDS with the snp.id, 
snp.index the index in GDS of the snp in at least one sample
phase information for the snp.id of the 1KG 
in R

   
    fileGDS1KG <- file.path(PATH1KG, PATHGDS, "matGeno1000g.gds")
    gds <- openfn.gds(fileGDS1KG)
    fileGDSPhase <- file.path(PATH1KG, PATHGDS, "matPhase1000g.gds")
    gdsPhase <- openfn.gds(fileGDSPhase)
    
    PATHSAMPLEGDS <- file.path("data", "genoGDS1KG.TCGA-OV.samples")
    
    
    addPhase1KG2SampleGDSFromGDS(gds,
                                 gdsPhase,
                                 PATHSAMPLEGDS)
    closefn.gds(gdsPhase)
    closefn.gds(gds)



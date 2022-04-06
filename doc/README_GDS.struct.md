# GDS file structure


## The base for 1KG

    File: /mnt/wigclust5/data/unsafe/belleau/process1000G/samples1000gUnrelated/data/genoGDS1KG.2022.03.22/matGeno1000g.gds (15.4G)
    +    [  ] *
    |--+ sample.id   { Str8 2548, 19.9K }
    |--+ sample.annot   [ data.frame ] *
    |  |--+ sex   { Str8 2548, 5.0K }
    |  |--+ pop.group   { Str8 2548, 10.0K }
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

If we add the information relative to specific studies we want to
infer the ancestry.

    File: /mnt/wigclust15/data/unsafe/belleau/TCGA-OV-WXS/data/genoGDS1KG.TCGA-OV.2022.03.30/matGeno1000g.gds (15.5G)
    +    [  ] *
    |--+ sample.id   { Str8 2558, 20.2K }
    |--+ sample.annot   [ data.frame ] *
    |  |--+ sex   { Str8 2558, 5.0K }
    |  |--+ pop.group   { Str8 2558, 10.2K }
    |  |--+ superPop   { Str8 2558, 10.0K }
    |  \--+ batch   { Float64 2558, 20.0K }
    |--+ snp.id   { Str8 24516859, 223.2M }
    |--+ snp.chromosome   { UInt16 24516859, 46.8M }
    |--+ snp.position   { Int32 24516859, 93.5M }
    |--+ snp.allele   { Str8 24516859, 93.5M }
    |--+ snp.AF   { PackedReal24 24516859, 70.1M }
    |--+ snp.EAS_AF   { PackedReal24 24516859, 70.1M }
    |--+ snp.EUR_AF   { PackedReal24 24516859, 70.1M }
    |--+ snp.AFR_AF   { PackedReal24 24516859, 70.1M }
    |--+ snp.AMR_AF   { PackedReal24 24516859, 70.1M }
    |--+ snp.SAS_AF   { PackedReal24 24516859, 70.1M }
    |--+ genotype   { Bit2 24516859x2558, 14.6G }
    |--+ sample.ref   { Bit1 2548, 319B }
    |--+ study.offset   { Int32 1, 4B }
    |--+ study.list   [ data.frame ] *
    |  |--+ study.id   { Str8 1, 14B }
    |  |--+ study.desc   { Str8 1, 16B }
    |  \--+ study.platform   { Str8 1, 4B }
    \--+ study.annot   [ data.frame ] *
       |--+ data.id   { Str8 10, 170B }
       |--+ case.id   { Str8 10, 130B }
       |--+ sample.type   { Str8 10, 140B }
       |--+ diagnosis   { Str8 10, 20B }
       |--+ source   { Str8 10, 60B }
       \--+ study.id   { Str8 10, 140B }


## The phase GDS for 1KG

It includes the linkage disequilibrium and the gene info
for the allelic of the RNA seq (Not show yet)

    File: /mnt/wigclust6/data/unsafe/belleau/process1000G/samples1000gUnrelated/data/genoGDS1KG.2022.04.03/matPhase1000g.gds (5.1G)
    +    [  ]
    \--+ phase   { Bit2 24516859x2548 LZ4_ra(35.0%), 5.1G }


## A group of GDS specific to each sample we want to infor

Principaly use for the simulation

    File: /mnt/wigclust15/data/unsafe/belleau/TCGA-OV-WXS/data/genoGDS1KG.TCGA-OV.samples/TCGA-25-2404-01A-01W-0799-08.gds (164.1M)
    +    [  ]
    |--+ Ref.count   { SparseInt16 24516859x1, 7.2M }
    |--+ Alt.count   { SparseInt16 24516859x1, 1.0M }
    |--+ Total.count   { SparseInt16 24516859x1, 7.6M }
    |--+ pruned.study   { Str8 237908, 2.2M }
    |--+ sampleStudy   { Str8 1, 29B }
    |--+ sample.id   { Str8 2470, 19.3K }
    |--+ snp.id   { Str8 237908, 2.2M }
    |--+ snp.chromosome   { Int32 237908, 929.3K }
    |--+ snp.position   { Int32 237908, 929.3K }
    |--+ snp.index   { Int32 237908, 929.3K }
    |--+ genotype   { Bit2 237908x2470, 140.1M }
    |--+ SamplePos   { Float64 1, 8B }
    |--+ lap   { PackedReal8 237908, 232.3K }
    \--+ segment   { UInt32 237908, 929.3K }


## A GDS phase specific for a study with with a smaller number of SNV

The snp.id are the union the snp.id of GDSSample

File: /mnt/wigclust15/data/unsafe/belleau/TCGA-OV-WXS/data/genoGDS1KG.TCGA-OV.samples/phase1KG.gds (434.6M)
+    [  ]
|--+ snp.id   { Str8 722473, 6.6M }
|--+ snp.index   { Int32 722473, 2.8M }
\--+ phase   { Bit2 722473x2469, 425.3M }



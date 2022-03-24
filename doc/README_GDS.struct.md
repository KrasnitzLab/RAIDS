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




############################################################
## How to create the mapSNVSelected_Demo.rds file
## This file contains a data.frame with the SNV allelic
## frequency in the continental populations
############################################################

selectedSNVs <- data.frame(CHROM=rep("chr1", 7),
                    POS=c(16102, 51478, 51897, 51927, 54489, 54707, 54715),
                    REF=c("T", "T", "C", "G", "G", "G", "C"),
                    ALT=c("G", "A", "A", "A", "A", "C", "T"),
                    AF=c(0.02, 0.11, 0.08, 0.07, 0.1, 0.23, 0.21),
                    EAS_AF=c(0.0, 0.0, 0.05, 0.01, 0.0, 0.08, 0.07),
                    EUR_AF=c(0.04, 0.2, 0.12, 0.14, 0.18, 0.38, 0.34),
                    AFR_AF=c(0.03, 0.02, 0.05, 0.06, 0.02, 0.18, 0.16),
                    AMR_AF=c(0.03, 0.12, 0.06, 0.07, 0.1, 0.25, 0.24),
                    SAS_AF=c(0.02, 0.22, 0.1, 0.09, 0.21, 0.28, 0.27),
                    stringsAsFactors=FALSE)

saveRDS(selectedSNVs, file="mapSNVSelected_Demo.rds")

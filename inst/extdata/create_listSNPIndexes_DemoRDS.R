
############################################################
## How to create the listSNPIndexes_Demo.rds file
## This file contains a vector containing indexes
## of the retained SNVs
############################################################

retainedSNVs <- c(1, 2, 3, 4, 6, 8, 9)

saveRDS(retainedSNVs, file="listSNPIndexes_Demo.rds")

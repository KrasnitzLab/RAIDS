
############################################################
## How to create the unrelatedPatientsInfo_Demo.rds file
## This file contains a list with the unrelated/related
## status for Reference samples
############################################################

unrelatedPatients <- list()
unrelatedPatients$rels <- c("HG00109", NA)

unrelatedPatients$unrels <- rep(NA, 9)

saveRDS(unrelatedPatients, file="unrelatedPatientsInfo_Demo.rds")

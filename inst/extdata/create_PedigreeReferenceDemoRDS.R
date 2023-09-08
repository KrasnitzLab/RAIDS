
############################################################
## How to create the PedigreeReferenceDemo.rds file
## This file contains a small pedigree reference data.frame
############################################################

pedigree <- data.frame(sample.id=c("HG00100", "HG00101", "HG00102",
                    "HG00103", "HG00104", "HG00105", "HG00106", "HG00107",
                    "HG00108", "HG00109"),
                Name.ID=c("HG00100", "HG00101", "HG00102", "HG00103",
                    "HG00104", "HG00105", "HG00106", "HG00107",
                    "HG00108", "HG00109"),
                sex=c("1", "2", "2", "1", "2", "1", "2", "1", "2", "2"),
                pop.group=rep("ABC", 10), superPop=rep("AFR", 10),
                batch=rep(0, 10), stringsAsFactors=FALSE)

saveRDS(pedigree, file="PedigreeReferenceDemo.rds")


############################################################
## How to create the PedigreeDemo.ped file
## This file contains a small pedigree reference data.frame
############################################################

pedigreeDemo <- data.frame(Family.ID=c("BB01", "BB01", "BB01", "BB02", "BB02",
                        "BB02", "BB03", "BB03", "BB03", "BB04"),
                Individual.ID=c("HG00100", "HG00101", "HG00102", "HG00103",
                        "HG00104", "HG00105", "HG00106", "HG00107", "HG00108",
                        "HG00109"),
                Paternal.ID=c("0", "0", "HG00100", "0", "0", "HG00103",
                        "HG00107", "0", "0", "0"),
                Maternal.ID=c("0", "0", "HG00101", "0", "0", "HG00104",
                              "HG00108", "0", "0", "0"),
                Gender=c(1, 2, 2, 1, 2, 1, 2, 1, 2, 2),
                Phenotype=rep(0, 10),
                Population=rep("ACB", 10),
                Relationship=c("father", "mother", "child", "father",
                        "mother", "child", "child", "father", "mother",
                        "mother"),
                Siblings=rep(0, 10),
                Second.Order=rep(0, 10),
                Third.Order=rep(0, 10),
                Other.Comments=rep(0, 10),
                stringsAsFactors=FALSE)

write.table(pedigreeDemo, file="PedigreeDemo.ped", quote=FALSE,
                row.names=FALSE, col.names=TRUE, sep="\t")

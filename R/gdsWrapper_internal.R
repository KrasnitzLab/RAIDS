#' @title Initialization of the section related to the sample
#' information in the GDS file
#'
#' @description This function initializes the section related to the sample
#' information in the GDS file. The information is extracted from
#' the \code{data.frame} passed to the function. The nodes "sample.id" and
#' "sample.annot" are created in the GDS file.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param pedDF a \code{data.frame} containing the information related to the
#' samples. It must have those columns: "sample.id", "Name.ID", "sex",
#' "pop.group", "superPop" and "batch". All columns, except "sex" and batch",
#' are \code{character} strings. The "batch" and "sex" columns are
#' \code{integer}. The unique identifier
#' of this \code{data.frame} is the "Name.ID" column. The row names of the
#' \code{data.frame} must correspond to the identifiers present in the
#' "Name.ID" column.
#'
#' @param listSamples a \code{vector} of \code{character} string representing
#' the identifiers of the selected samples. If \code{NULL}, all samples are
#' selected. Default: \code{NULL}.
#'
#' @return a \code{vector} of \code{character} string with the identifiers of
#' the samples saved in the GDS file.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP_10.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create "sample.annot" node (the node must be present)
#' pedInformation <- data.frame(sample.id=c("sample_01", "sample_02"),
#'     Name.ID=c("sample_01", "sample_02"),
#'     sex=c(1,1),  # 1:Male  2: Female
#'     pop.group=c("ACB", "ACB"),
#'     superPop=c("AFR", "AFR"),
#'     batch=c(1, 1),
#'     stringsAsFactors=FALSE)
#'
#' ## The row names must be the sample identifiers
#' rownames(pedInformation) <- pedInformation$Name.ID
#'
#' ## Add information about 2 samples to the GDS file
#' RAIDS:::generateGDSSample(gds=GDS_file_tmp, pedDF=pedInformation,
#'     listSamples=NULL)
#'
#' ## Read sample identifier list
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.id"))
#'
#' ## Read sample information from GDS file
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.annot"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
generateGDSSample <- function(gds, pedDF, listSamples=NULL) {

    if(!(is.null(listSamples))){
        pedDF <- pedDF[listSamples,]
    }

    add.gdsn(node=gds, name="sample.id", val=pedDF[, "Name.ID"])

    ## Create a data.frame containing the information form the samples
    samp.annot <- data.frame(sex=pedDF[, "sex"],
        pop.group=pedDF[, "pop.group"], superPop=pedDF[, "superPop"],
        batch=pedDF[, "batch"],  stringsAsFactors=FALSE)

    ## Add the data.frame to the GDS object
    add.gdsn(node=gds, name="sample.annot", val=samp.annot)

    ## Return the vector of saved samples
    return(pedDF[, "Name.ID"])
}


#' @title Create a "sample.ref" node i a GDS file with the information about
#' the related/unrelated state of the reference samples
#'
#' @description This function creates a "sample.ref" node in the GDS file.
#' The node contains a vector of integers with value of 1 when
#' the samples are used as references and 0 otherwise.
#' The information used to fill the "sample.ref" node comes from the RDS file
#' that contains the information about the unrelated reference samples.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param filePart a \code{character} string representing the path and file
#' name of a RDS file containing the information about the related and
#' unrelated samples in the reference dataset. The RDS file must exist. The
#' RDS file must contains a \code{vector} of \code{character} strings called
#' "unrels" with the name of the unrelated samples.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Locate RDS with unrelated/related status for 1KG samples
#' data.dir <- system.file("extdata", package="RAIDS")
#' rdsFilePath <- file.path(data.dir, "unrelatedPatientsInfo_Demo.rds")
#'
#' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP_11.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#
#' ## Create "sample.id" node (the node must be present)
#' sampleIDs <- c("HG00104", "HG00109", "HG00110")
#' add.gdsn(node=GDS_file_tmp, name="sample.id", val=sampleIDs)
#'
#' ## Create  "sample.ref" node in GDS file using RDS information
#' RAIDS:::addGDSRef(gds=GDS_file_tmp, filePart=rdsFilePath)
#'
#' ## Read sample reference data.frame
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.ref"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
addGDSRef <- function(gds, filePart) {

    part <- readRDS(filePart)

    sampleGDS <- index.gdsn(gds, "sample.id")
    df <- data.frame(sample.id=read.gdsn(sampleGDS),
                        sample.ref=0, stringsAsFactors=FALSE)

    # The order of part$unrels is not the same than df$sample.id
    df[df$sample.id %in% part$unrels, "sample.ref"] <- 1
    add.gdsn(gds, "sample.ref", df$sample.ref, storage="bit1")

    ## Success
    return(0L)
}


#' @title Appends the genotype information for sepcific samples
#' (1 column == 1 sample) into a GDS file
#'
#' @description This function appends the genotype information into a
#' GDS file. More specifically, the genotype information is added to the
#' "genotype" node. The "genotype" node must already be present in the
#' GDS file. The genotype information is a matrix with the rows corresponding
#' to SNVs and columns corresponding to samples.
#' The number of rows of the new genotype information must
#' correspond to the number of rows of the matrix present in the
#' "genotype" node.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS Sample file.
#'
#' @param matG a \code{matrix} of \code{integer} representing the genotypes
#' of the SNVs for one or multiple samples. The rows correspond to SNVs and
#' the columns correspond to samples. The number of rows must
#' correspond to the number of rows of the matrix present in the
#' "genotype" node.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP_06.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create a "genotype" node with initial matrix
#' geno_initial <- matrix(rep(0L, 10), nrow=2)
#'
#' add.gdsn(node=GDS_file_tmp, name="genotype", val=geno_initial)
#' sync.gds(GDS_file_tmp)
#'
#' ## New genotype information to be added
#' new_genotype <- matrix(rep(1L, 6), nrow=2)
#'
#' ## Add segments to the GDS file
#' RAIDS:::appendGDSgenotypeMat(gds=GDS_file_tmp, matG=new_genotype)
#'
#' ## Read genotype information from GDS file
#' ## The return matrix should be a combinaison of both initial matrix
#' ## and new matrix (column binded)
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="genotype"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils read.csv2
#' @encoding UTF-8
#' @keywords internal
appendGDSgenotypeMat <- function(gds, matG) {

    geno.var <- index.gdsn(gds, "genotype")
    append.gdsn(node=geno.var, val=matG, check=TRUE)

    return(0L)
}


#' @title TODO This function append the genotype and the file related to the
#' pileup
#'
#' @description TODO
#'
#' @param pathGeno a \code{character} string representing the path to a
#' directory with the genotype files for the profiles, as generated by
#' snp-pileup. The path must exist.
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to the sample identifiers that will have a GDS Sample file created. The
#' sample identifiers must be present in the "Name.ID" column of the
#' \code{data.frame} passed to the \code{pedDF} parameter.
#'
#' @param listPos  a \code{array} TODO
#'
#' @param offset  a \code{integer} to adjust if the genome start at 0 or 1
#'
#' @param minCov  a single positive \code{integer} representing the minimum
#' coverage needed to keep the SNPs in the analysis. Default: \code{10}.
#'
#' @param minProb  a \code{array} with the sample to keep
#'
#' @param seqError a single positive \code{numeric} between 0 and 1
#' representing the sequencing error rate. Default: \code{0.001}.
#'
#' @param pedDF a \code{data.frame} with the information about the sample(s).
#' Those are mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", "Source". All columns must be in
#' \code{character} strings. The \code{data.frame}
#' must contain the information for all the samples passed in the
#' \code{listSamples} parameter.
#'
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore.
#'
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings.
#'
#' @param PATHGDSSAMPLE TODO a PATH to a directory where a gds specific
#' to the samples with coverage info is keep
#'
#' @param genoSource a \code{stirng} with two possible values:
#' snp-pileup and generic. It specify if the genotype files
#' are generate by snp-pileup(Facets) or generic format csv
#' with the column at least the columns:
#' Chromosome,Position,Ref,Alt,Count,File1R,File1A
#' where Count is the deep at the position,
#' FileR is the deep of the reference allele, and
#' File1A is the deep of the specific alternative allele
#'
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running.
#'
#' @return The  function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#' gds <- "Demo GDS TODO"
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn write.gdsn openfn.gds
#' @importFrom stats qbinom
#' @importFrom utils read.csv
#' @encoding UTF-8
#' @keywords internal
generateGDS1KGgenotypeFromSNPPileup <- function(pathGeno,
    listSamples, listPos, offset, minCov=10, minProb=0.999,
    seqError=0.001, pedStudy, batch, studyDF, PATHGDSSAMPLE=NULL,
    genoSource=c("snp-pileup", "generic"), verbose) {

    # File with the description of the SNP keep
    listMat <- dir(pathGeno, pattern = ".+.txt.gz")
    listSampleFile <- gsub(".txt.gz", "", listMat)

    g <- as.matrix(rep(-1, nrow(listPos)))

    for(i in seq_len(length(listSamples))) {
        pos <- which(listSampleFile == listSamples[i])

        if(verbose) { message(listSamples[i], " ", Sys.time(), " ", i) }

        if(length(pos) == 1) {

            #
            if(genoSource == "snp-pileup") {
                matSample <- readSNVPileupFile(file.path(pathGeno,
                                        listMat[pos]), offset)
            } else if(genoSource == "generic") {
                matSample <- readSNVFileGeneric(file.path(pathGeno,
                                        listMat[pos]), offset)
            }
            # matAll <- merge(matSample[,c( "Chromosome", "Position",
            #                               "File1R",  "File1A",
            #                               "count" )],
            #                 listPos,
            #                 by.x = c("Chromosome", "Position"),
            #                 by.y = c("snp.chromosome", "snp.position"),
            #                 all.y = TRUE,
            #                 all.x = FALSE)
            #
            # below same as the merge above but faster

            z <- cbind(c(listPos$snp.chromosome, matSample$Chromosome,
                            matSample$Chromosome),
                        c(listPos$snp.position, matSample$Position,
                            matSample$Position),
                        c(rep(1,nrow(listPos)), rep(0,nrow(matSample)),
                            rep(2,nrow(matSample))),
                        c(rep(0,nrow(listPos)), matSample[, "File1R"],
                            -1 * matSample[, "File1R"]),
                        c(rep(0,nrow(listPos)), matSample[, "File1A"],
                            -1 * matSample[, "File1A"]),
                        c(rep(0,nrow(listPos)), matSample[, "count"],
                            -1 * matSample[, "count"]))
            rm(matSample)
            z <- z[order(z[,1], z[,2], z[,3]),]

            matAll <- data.frame(Chromosome=z[z[, 3] == 1, 1],
                Position=z[z[, 3] == 1, 2], File1R=cumsum(z[, 4])[z[, 3] == 1],
                File1A=cumsum(z[,5])[z[, 3] == 1],
                count=cumsum(z[, 6])[z[, 3] == 1])
            rm(z)

            if(is.null(PATHGDSSAMPLE)){
                stop("PATHGDSSAMPLE is NULL in ",
                        "generateGDS1KGgenotypeFromSNPPileup\n")
            } else{
                if(! dir.exists(PATHGDSSAMPLE)) {
                    dir.create(PATHGDSSAMPLE)
                }
            }
            fileGDSSample <- file.path(PATHGDSSAMPLE,
                                        paste0(listSamples[i], ".gds"))
            if(file.exists(fileGDSSample)) {
                gdsSample <- openfn.gds(fileGDSSample, readonly=FALSE)
            } else{
                gdsSample <- createfn.gds(fileGDSSample)
            }

            if (! "Ref.count" %in% ls.gdsn(gdsSample)) {
                var.Ref <- add.gdsn(gdsSample, "Ref.count", matAll$File1R,
                                        valdim=c( nrow(listPos), 1),
                                        storage="sp.int16")
                var.Alt <- add.gdsn(gdsSample, "Alt.count", matAll$File1A,
                                        valdim=c( nrow(listPos), 1),
                                        storage="sp.int16")
                var.Count <- add.gdsn(gdsSample, "Total.count",
                                matAll$count, valdim=c( nrow(listPos), 1),
                                storage="sp.int16")
            } else {
                # you must append
                var.Ref <- append.gdsn(index.gdsn(gdsSample, "Ref.count"),
                                            matAll$File1R)
                var.Alt <- append.gdsn(index.gdsn(gdsSample, "Alt.count"),
                                            matAll$File1A)
                var.Count <- append.gdsn(index.gdsn(gdsSample, "Total.count"),
                                            matAll$count)
            }

            listSampleGDS <- addStudyGDSSample(gdsSample, pedDF=pedStudy,
                                batch=batch, listSamples=c(listSamples[i]),
                                studyDF=studyDF)

            listCount <- table(matAll$count[matAll$count >= minCov])
            cutOffA <-
                data.frame(count=unlist(vapply(as.integer(names(listCount)),
                    FUN=function(x, minProb, eProb){
                        return(max(2,qbinom(minProb, x,eProb))) },
                    FUN.VALUE=numeric(1), minProb=minProb, eProb=2 * seqError)),
                        allele=unlist(vapply(as.integer(names(listCount)),
                    FUN=function(x, minProb, eProb){
                        return(max(2,qbinom(minProb, x,eProb))) },
                    FUN.VALUE=numeric(1), minProb=minProb, eProb=seqError)))

            row.names(cutOffA) <- names(listCount)
            # Initialize the genotype array at -1

            # Select the position where the coverage of the 2 alleles is enough
            listCov <- which(rowSums(matAll[, c("File1R", "File1A")]) >= minCov)

            matAllC <- matAll[listCov,]

            # The difference  depth - (nb Ref + nb Alt) can be realistically
            # explain by sequencing error
            listCov <- listCov[(matAllC$count -
                                    (matAllC$File1R + matAllC$File1A)) <
                            cutOffA[as.character(matAllC$count), "count"]]

            matAllC <- matAll[listCov,]
            rm(matAll)

            g <- as.matrix(rep(-1, nrow(listPos)))
            # The sample is homozygote if the other known allele have a
            # coverage of 0
            g[listCov][which(matAllC$File1A == 0)] <- 0
            g[listCov][which(matAllC$File1R == 0)] <- 2

            # The sample is heterozygote if explain the coverage of
            # the minor allele by sequencing error is not realistic.
            g[listCov][which(matAllC$File1A >=
                        cutOffA[as.character(matAllC$count), "allele"] &
                        matAllC$File1R >= cutOffA[as.character(matAllC$count),
                                                        "allele"])] <- 1

            #g <- as.matrix(g)
            if("geno.ref" %in% ls.gdsn(gdsSample)){
                var.geno <- index.gdsn(gdsSample, "geno.ref")

                compression.gdsn(var.geno, compress="")
                #write.gdsn(cur, x, start=c(1, 2), count=c(-1,1))
                append.gdsn(var.geno,g, check=TRUE)
                compression.gdsn(var.geno, compress="LZMA_RA.fast")
                readmode.gdsn(var.geno)

            }else{
                var.geno <- add.gdsn(gdsSample, "geno.ref",
                                valdim=c(length(g), 1),
                                g, storage="bit2", compress = "LZMA_RA.fast")
                readmode.gdsn(var.geno)
            }

            rm(g)
            closefn.gds(gdsfile=gdsSample)

            if (verbose) {
                message(listMat[pos], " ", i, " ", Sys.time())
            }

        } else {
            stop("Missing samples ", listSamples[i])
        }
    }
    ## Success
    return(0L)
}




#' @title Identity-by-descent (IBD) analysis
#'
#' @description This function calculates the IDB coefficients by KING method
#' of moment using the
#' \code{\link[SNPRelate:snpgdsIBDKING]{SNPRelate::snpgdsIBDKING}}
#' function.
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, an
#' opened SNP GDS file.
#'
#' @param profileID  a \code{vector} of \code{character} strings representing
#' the samples to keep for the analysis. If \code{NULL}, all samples are used.
#' Default: \code{NULL}.
#'
#' @param snpID  a \code{vector} of \code{character} strings representing
#' the SNPs to keep for the analysis. If \code{NULL}, all SNPs are used.
#' Default: \code{NULL}.
#'
#' @param maf  a single \code{numeric} representing the threshold for the minor
#' allele frequency. Only the SNPs with ">= maf" are retained.
#' Default: \code{0.05}.
#'
#' @param verbose a \code{logical} indicating if information is shown during
#' the process in the \code{\link[SNPRelate]{snpgdsIBDKING}}() function.
#'
#' @return a \code{list} containing:
#' \itemize{
#'     \item{sample.id}{a \code{character} string representing the sample
#'     ids used in the analysis}
#'     \item{snp.id}{a \code{character} string representing the SNP ids
#'     used in the analysis}
#'     \item{k0}{a \code{numeric}, the IBD coefficient, the probability of
#'     sharing zero IBD}
#'     \item{k1}{a \code{numeric}, the IBD coefficient, the probability of
#'     sharing one IBD}
#'     \item{IBS0}{a \code{numeric}, the proportion of SNPs with zero IBS}
#'     \item{kinship}{a \code{numeric}, the proportion of SNPs with zero IBS,
#'     if the parameter kinship=TRUE}
#' }
#'
#' @examples
#'
#' ## Required
#' library(SNPRelate)
#'
#' ## Open an example dataset (HapMap)
#' genoFile <- snpgdsOpen(snpgdsExampleFileName())
#'
#' ## Extract CEU population
#' samples <- read.gdsn(index.gdsn(genoFile, "sample.id"))
#' CEU <- samples[
#'     read.gdsn(index.gdsn(genoFile, "sample.annot/pop.group"))=="CEU"]
#'
#' ## Infer the presence of population stratification
#' ibd.robust <- RAIDS:::runIBDKING(gds=genoFile, profileID=CEU, snpID=NULL,
#'      maf=0.05, verbose=FALSE)
#'
#' ## close the genotype file
#' snpgdsClose(genoFile)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom SNPRelate snpgdsIBDKING
#' @encoding UTF-8
#' @keywords internal
runIBDKING <- function(gds, profileID=NULL, snpID=NULL, maf=0.05, verbose) {

    ## Calculate IBD coefficients by KING method of moment
    ibd.robust <- snpgdsIBDKING(gdsobj=gds, sample.id=profileID,
                        snp.id=snpID, maf=maf,
                        type="KING-robust",
                        verbose=verbose)

    ## Return result
    return(ibd.robust)
}


#' @title SNP pruning based on linkage disequilibrium (LD)
#'
#' @description This function is a wrapper for the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function that generates a pruned
#' subset of SNPs that are in approximate linkage equilibrium.
#'
#' @param gds an \code{object} of class
#' \code{\link[SNPRelate]{SNPGDSFileClass}}, a SNP GDS file.
#'
#' @param method a \code{character} string that represents the method that will
#' be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime" and "composite". Default: \code{"corr"}.
#'
#' @param listSamples a \code{vector} of \code{character} strings
#' corresponding to the sample identifiers used in LD pruning done by the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. If \code{NULL}, all
#' samples are used. Default: \code{NULL}.
#'
#' @param listKeep a \code{vector} of SNVs identifiers specifying selected;
#' if \code{NULL}, all SNVs are used in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}} function. Default: \code{NULL}.
#'
#' @param slideWindowMaxBP a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done in the \code{\link[SNPRelate]{snpgdsLDpruning}}()
#' function.
#' Default: \code{500000L}.
#'
#' @param thresholdLD a single \code{numeric} value that represents the LD
#' threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#' Default: \code{sqrt(0.1)}.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param verbose a \code{logical} indicating if information is shown during
#' the process in the \code{\link[SNPRelate]{snpgdsLDpruning}}() function.
#'
#' @return a \code{list} of SNP identifiers stratified by chromosomes as
#' generated by \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#'
#' @details
#'
#' The SNP pruning is based on linkage disequilibrium (LD) and is done by the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function in the
#' SNPRelate package (https://bioconductor.org/packages/SNPRelate/).
#'
#' @examples
#'
#' ## Open an example dataset (HapMap)
#' genoFile <- snpgdsOpen(snpgdsExampleFileName())
#'
#' ## Fix seed to get reproducible results
#' set.seed(1000)
#'
#' ## Get linkage Disequilibrium (LD) based SNP pruning
#' snpSet <- RAIDS:::runLDPruning(gds=genoFile, verbose=FALSE)
#' names(snpSet)
#'
#' ## Get SNP ids
#' snp.id <- unlist(unname(snpSet))
#'
#' ## Close the genotype file
#' snpgdsClose(genoFile)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#'
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpruning
#' @importFrom gdsfmt closefn.gds
#' @encoding UTF-8
#' @keywords internal
runLDPruning <- function(gds, method=c("corr", "r", "dprime", "composite"),
        listSamples=NULL, listKeep=NULL, slideWindowMaxBP=500000L,
        thresholdLD=sqrt(0.1), np=1L, verbose) {

    ## Call SNP LD pruning
    result <- snpgdsLDpruning(gdsobj=gds, method="corr",
                                sample.id=listSamples,
                                snp.id=listKeep, slide.max.bp=slideWindowMaxBP,
                                ld.threshold=thresholdLD, num.thread=np,
                                verbose=verbose)

    return(result)
}


#' @title Append fields related to samples into a GDS file
#'
#' @description This function appends the fields related to samples into
#' a GDS file. The information is extracted from the \code{data.frame} passed
#' to the function and is added to the "sample.annot" and "sample.id" nodes.
#' The "sample.id" and "sample.annot" nodes must already exist.
#' If the samples are part of a study, the function
#' addStudyGDSSample() must be used.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param pedDF a \code{data.frame} with the information about the sample(s).
#' The \code{data.frame} must have the columns: "sample.id", "Name.ID", "sex",
#' "pop.group" and "superPop". The unique identifier for the sample(s) is
#' the "Name.ID" column and the row names of the \code{data.frame} must
#' correspond to the "Name.ID" column.
#'
#' @param batch a \code{integer} representing the batch identifier.
#'
#' @param listSamples a \code{vector} of \code{character} string with the
#' selected sample(s). If \code{NULL}, all samples are used.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{TRUE}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP_03.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create "sample.id" node (the node must be present)
#' add.gdsn(node=GDS_file_tmp, name="sample.id", val=c("sample_01",
#'     "sample_02"))
#'
#' ## Create "sample.annot" node (the node must be present)
#' add.gdsn(node=GDS_file_tmp, name="sample.annot", val=data.frame(
#'   Name.ID=c("sample_01", "sample_02"),
#'     sex=c(1,1),  # 1:Male  2: Female
#'     pop.group=c("ACB", "ACB"),
#'     superPop=c("AFR", "AFR"),
#'     batch=c(1, 1),
#'     stringsAsFactors=FALSE))
#'
#' sync.gds(gdsfile=GDS_file_tmp)
#'
#' ## Create a data.frame with information about samples
#' sample_info <- data.frame(Name.ID=c("sample_04", "sample_05", "sample_06"),
#'     sex=c(1,2,1),  # 1:Male  2: Female
#'     pop.group=c("ACB", "ACB", "ACB"),
#'     superPop=c("AFR", "AFR", "AFR"),
#'     stringsAsFactors=FALSE)
#'
#' ## The row names must be the sample identifiers
#' rownames(sample_info) <- sample_info$Name.ID
#'
#' ## Add information about 2 samples to the GDS file
#' RAIDS:::appendGDSSample(gds=GDS_file_tmp, pedDF=sample_info, batch=2,
#'     listSamples=c("sample_04", "sample_06"), verbose=FALSE)
#'
#' ## Read sample identifier list
#' ## Only "sample_04" and "sample_06" should have been added
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.id"))
#'
#' ## Read sample information from GDS file
#' ## Only "sample_04" and "sample_06" should have been added
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.annot"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn append.gdsn
#' @encoding UTF-8
#' @keywords internal
appendGDSSample <- function(gds, pedDF, batch=1, listSamples=NULL,
                                verbose=TRUE) {

    ## Only keep selected samples
    if(!(is.null(listSamples))){
        pedDF <- pedDF[listSamples,]
    }

    ## Append sample identifiers to the "sample.id" node
    sampleGDS <- index.gdsn(gds, "sample.id")
    append.gdsn(sampleGDS, val=pedDF$Name.ID, check=TRUE)

    ## Create the data.frame with the sample information
    samp.annot <- data.frame(sex = pedDF[, "sex"],
                                pop.group=pedDF[, "pop.group"],
                                superPop=pedDF[, "superPop"],
                                batch=rep(batch, nrow(pedDF)),
                                stringsAsFactors=FALSE)

    if(verbose) { message("Annot") }

    ## Append data.frame to "sample.annot" node
    curAnnot <- index.gdsn(gds, "sample.annot/sex")
    append.gdsn(curAnnot, samp.annot$sex, check=TRUE)
    curAnnot <- index.gdsn(gds, "sample.annot/pop.group")
    append.gdsn(curAnnot, samp.annot$pop.group, check=TRUE)
    curAnnot <- index.gdsn(gds, "sample.annot/superPop")
    append.gdsn(curAnnot, samp.annot$superPop, check=TRUE)
    curAnnot <- index.gdsn(gds, "sample.annot/batch")
    append.gdsn(curAnnot, samp.annot$batch, check=TRUE)

    if(verbose) { message("Annot done") }

    return(0L)
}

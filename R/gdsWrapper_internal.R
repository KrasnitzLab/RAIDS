#' @title Initialization of the section related to the profile
#' information in the GDS file
#'
#' @description This function initializes the section related to the profile
#' information in the GDS file. The information is extracted from
#' the \code{data.frame} passed to the function. The nodes "sample.id" and
#' "sample.annot" are created in the GDS file.
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param dfPedReference a \code{data.frame} containing the information
#' related to the
#' samples. It must have those columns: "sample.id", "Name.ID", "sex",
#' "pop.group", "superPop" and "batch". All columns, except "sex" and batch",
#' are \code{character} strings. The "batch" and "sex" columns are
#' \code{integer}. The unique identifier
#' of this \code{data.frame} is the "Name.ID" column. The row names of the
#' \code{data.frame} must correspond to the identifiers present in the
#' "Name.ID" column.
#'
#' @param listSamples a \code{vector} of \code{character} string representing
#' the identifiers of the selected profiles. If \code{NULL}, all profiles are
#' selected. Default: \code{NULL}.
#'
#' @return a \code{vector} of \code{character} string with the identifiers of
#' the profiles saved in the GDS file.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Temporary GDS file in current directory
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_10.gds")
#'
#' ## Create and open the GDS file
#' tmpGDS  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create "sample.annot" node (the node must be present)
#' pedInformation <- data.frame(sample.id=c("sample_01", "sample_02"),
#'         Name.ID=c("sample_01", "sample_02"),
#'         sex=c(1,1),  # 1:Male  2: Female
#'         pop.group=c("ACB", "ACB"),
#'         superPop=c("AFR", "AFR"),
#'         batch=c(1, 1),
#'         stringsAsFactors=FALSE)
#'
#' ## The row names must be the sample identifiers
#' rownames(pedInformation) <- pedInformation$Name.ID
#'
#' ## Add information about 2 samples to the GDS file
#' RAIDS:::generateGDSRefSample(gdsReference=tmpGDS,
#'         dfPedReference=pedInformation, listSamples=NULL)
#'
#' ## Read sample identifier list
#' read.gdsn(index.gdsn(node=tmpGDS, path="sample.id"))
#'
#' ## Read sample information from GDS file
#' read.gdsn(index.gdsn(node=tmpGDS, path="sample.annot"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=tmpGDS)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
generateGDSRefSample <- function(gdsReference, dfPedReference,
                                        listSamples=NULL) {

    if(!(is.null(listSamples))){
        dfPedReference <- dfPedReference[listSamples,]
    }

    add.gdsn(node=gdsReference, name="sample.id",
                    val=dfPedReference[, "Name.ID"])

    ## Create a data.frame containing the information form the samples
    samp.annot <- data.frame(sex=dfPedReference[, "sex"],
        pop.group=dfPedReference[, "pop.group"],
        superPop=dfPedReference[, "superPop"],
        batch=dfPedReference[, "batch"],  stringsAsFactors=FALSE)

    ## Add the data.frame to the GDS object
    add.gdsn(node=gdsReference, name="sample.annot", val=samp.annot)

    ## Return the vector of saved samples
    return(dfPedReference[, "Name.ID"])
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
#' @param gdsReference an object of class
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
#' ## Required library
#' library(gdsfmt)
#'
#' ## Locate RDS with unrelated/related status for Reference samples
#' dataDir <- system.file("extdata", package="RAIDS")
#' rdsFilePath <- file.path(dataDir, "unrelatedPatientsInfo_Demo.rds")
#'
#' ## Temporary GDS file
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_11.gds")
#'
#' ## Create and open the GDS file
#' tmpGDS  <- createfn.gds(filename=gdsFilePath)
#
#' ## Create "sample.id" node (the node must be present)
#' sampleIDs <- c("HG00104", "HG00109", "HG00110")
#' add.gdsn(node=tmpGDS, name="sample.id", val=sampleIDs)
#'
#' ## Create "sample.ref" node in GDS file using RDS information
#' RAIDS:::addGDSRef(gdsReference=tmpGDS, filePart=rdsFilePath)
#'
#' ## Read sample reference data.frame
#' read.gdsn(index.gdsn(node=tmpGDS, path="sample.ref"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=tmpGDS)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
addGDSRef <- function(gdsReference, filePart) {

    part <- readRDS(filePart)

    sampleGDS <- index.gdsn(gdsReference, "sample.id")
    df <- data.frame(sample.id=read.gdsn(sampleGDS),
                        sample.ref=0, stringsAsFactors=FALSE)

    # The order of part$unrels is not the same than df$sample.id
    df[df$sample.id %in% part$unrels, "sample.ref"] <- 1
    add.gdsn(gdsReference, "sample.ref", df$sample.ref, storage="bit1")

    ## Success
    return(0L)
}


#' @title Appends the genotype information for specific samples
#' (1 column == 1 profile) into a GDS file
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
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Profile GDS file.
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
#' ## Required library
#' library(gdsfmt)
#'
#' ## Create a temporary GDS file
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_06.gds")
#'
#' ## Create and open the GDS file
#' tmpGDS  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create a "genotype" node with initial matrix
#' genoInitial <- matrix(rep(0L, 10), nrow=2)
#'
#' add.gdsn(node=tmpGDS, name="genotype", val=genoInitial)
#' sync.gds(tmpGDS)
#'
#' ## New genotype information to be added
#' newGenotype <- matrix(rep(1L, 6), nrow=2)
#'
#' ## Add segments to the GDS file
#' RAIDS:::appendGDSgenotypeMat(gds=tmpGDS, matG=newGenotype)
#'
#' ## Read genotype information from GDS file
#' ## The return matrix should be a combination of both initial matrix
#' ## and new matrix (column binded)
#' read.gdsn(index.gdsn(node=tmpGDS, path="genotype"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=tmpGDS)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
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


#' @title Append the genotype information from a profile into the associated
#' Profile GDS File
#'
#' @description This function append the genotype information from a specific
#' profile into the Profile GDS file. The genotype information is extracted
#' from a SNV file as generated by SNP-pileup or other tools.
#'
#' @param pathGeno a \code{character} string representing the path to the
#' directory containing the VCF output of SNP-pileup for each sample. The
#' SNP-pileup files must be compressed (gz files) and have the name identifiers
#' of the samples. A sample with "Name.ID" identifier would have an
#' associated file called
#' if genoSource is "VCF", then "Name.ID.vcf.gz",
#' if genoSource is "generic", then "Name.ID.generic.txt.gz"
#' if genoSource is "snp-pileup", then "Name.ID.txt.gz".
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to the sample identifiers that will have a Profile GDS file created. The
#' sample identifiers must be present in the "Name.ID" column of the
#' \code{data.frame} passed to the \code{dfPedProfile} parameter.
#'
#' @param listPos a \code{data.frame} containing 2 columns. The first column,
#' called "snp.chromosome" contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position" contains the
#' position of the SNV on the chromosome.
#'
#' @param offset a \code{integer} to adjust if the genome start at 0 or 1.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' coverage needed to keep the SNVs in the analysis. Default: \code{10}.
#'
#' @param minProb a single positive \code{numeric} between 0 and 1
#' representing the probability that the base change at the SNV position
#' is not an error.
#' Default: \code{0.999}.
#'
#' @param seqError a single positive \code{numeric} between 0 and 1
#' representing the sequencing error rate. Default: \code{0.001}.
#'
#' @param dfPedProfile a \code{data.frame} with the information about
#' the sample(s).
#' Those are mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis" and "Source". All columns must be in
#' \code{character} strings format. The \code{data.frame}
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
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Sample files will be created.
#'
#' @param genoSource a \code{character} string with two possible values:
#' 'snp-pileup', 'generic' or 'VCF'. It specifies if the genotype files
#' are generated by snp-pileup (Facets) or are a generic format CSV file
#' with at least those columns:
#' 'Chromosome', 'Position', 'Ref', 'Alt', 'Count', 'File1R' and 'File1A'.
#' The 'Count' is the depth at the specified position;
#' 'FileR' is the depth of the reference allele and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally the file can be a VCF file with at least those genotype
#' fields: GT, AD, DP.
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running.
#'
#' @return The  function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Current directory
#' dataDir <- file.path(tempdir())
#'
#' ## Copy required file into current directory
#' file.copy(from=file.path(system.file("extdata/tests", package="RAIDS"),
#'                     "ex1.txt.gz"), to=dataDir)
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyDF <- data.frame(study.id = "MYDATA",
#'                         study.desc = "Description",
#'                         study.platform = "PLATFORM",
#'                         stringsAsFactors = FALSE)
#'
#' ## The data.frame containing the information about the samples
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' samplePED <- data.frame(Name.ID=c("ex1", "ex2"),
#'                     Case.ID=c("Patient_h11", "Patient_h12"),
#'                     Diagnosis=rep("Cancer", 2),
#'                     Sample.Type=rep("Primary Tumor", 2),
#'                     Source=rep("Databank B", 2), stringsAsFactors=FALSE)
#' rownames(samplePED) <- samplePED$Name.ID
#'
#' ## List of SNV positions
#' listPositions <- data.frame(snp.chromosome=c(rep(1, 10)),
#'         snp.position=c(3467333, 3467428, 3469375, 3469387, 3469502, 3469527,
#'         3469737, 3471497, 3471565, 3471618))
#'
#' ## Append genotype information to the Profile GDS file
#' result <- RAIDS:::generateGDS1KGgenotypeFromSNPPileup(pathGeno=dataDir,
#'             listSamples=c("ex1"), listPos=listPositions,
#'             offset=-1, minCov=10, minProb=0.999, seqError=0.001,
#'             dfPedProfile=samplePED, batch=1, studyDF=studyDF,
#'             pathProfileGDS=dataDir, genoSource="snp-pileup",
#'             verbose=FALSE)
#'
#' ## The function returns OL when successful
#' result
#'
#' ## The Profile GDS file 'ex1.gds' has been created in the
#' ## specified directory
#' list.files(dataDir)
#'
#' ## Unlink Profile GDS file (created for demo purpose)
#' unlink(file.path(dataDir, "ex1.gds"))
#' unlink(file.path(dataDir, "ex1.txt.gz"))
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
    seqError=0.001, dfPedProfile, batch, studyDF, pathProfileGDS,
    genoSource, verbose) {

    # File with the description of the SNP keep
    if(genoSource == "VCF"){
        listMat <- dir(pathGeno, pattern = ".+.vcf.gz")
        listSampleFile <- gsub(".vcf.gz", "", listMat)
    }else{
        listMat <- dir(pathGeno, pattern = ".+.txt.gz")
        listSampleFile <- gsub(".txt.gz", "", listMat)
    }


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
            } else if(genoSource == "VCF") {
                tmpProfile <- gsub(".vcf.gz", "",listMat[pos])
                matSample <- readSNVVCF(file.path(pathGeno,
                                        listMat[pos]),
                                        profileName=tmpProfile, offset)
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

            if(is.null(pathProfileGDS)){
                stop("pathProfileGDS is NULL in ",
                        "generateGDS1KGgenotypeFromSNPPileup\n")
            } else{
                if(! dir.exists(pathProfileGDS)) {
                    dir.create(pathProfileGDS)
                }
            }
            fileGDSSample <- file.path(pathProfileGDS,
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

            listSampleGDS <- addStudyGDSSample(gdsSample,
                pedProfile=dfPedProfile, batch=batch,
                listSamples=c(listSamples[i]), studyDF=studyDF, verbose=verbose)

            listCount <- table(matAll$count[matAll$count >= minCov])
            cutOffA <-
                data.frame(count=unlist(vapply(as.integer(names(listCount)),
                    FUN=function(x, minProb, eProb){
                        return(max(2,qbinom(minProb, x, eProb))) },
                    FUN.VALUE=numeric(1), minProb=minProb, eProb=2 * seqError)),
                        allele=unlist(vapply(as.integer(names(listCount)),
                    FUN=function(x, minProb, eProb){
                        return(max(2,qbinom(minProb, x, eProb))) },
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
                append.gdsn(var.geno, g, check=TRUE)
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

#' @title Append the genotype information from a profile into the associated
#' Profile GDS File
#'
#' @description This function append the genotype information from a specific
#' profile into the Profile GDS file. The genotype information is extracted
#' from a SNV file as generated by SNP-pileup or other tools.
#'
#' @param profileFile a \code{character} string representing the path and the
#' file name of the genotype file or the bam if genoSource is snp-pileup the
#' fine extension must be .txt.gz, if VCF the extension must be .vcf.gz
#'
#' @param profileName a  \code{character} string representing the profileName
#'
#' @param listPos a \code{data.frame} containing 2 columns. The first column,
#' called "snp.chromosome" contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position" contains the
#' position of the SNV on the chromosome.
#'
#' @param offset a \code{integer} to adjust if the genome start at 0 or 1.
#'
#' @param minCov a single positive \code{integer} representing the minimum
#' coverage needed to keep the SNVs in the analysis. Default: \code{10}.
#'
#' @param minProb a single positive \code{numeric} between 0 and 1
#' representing the probability that the base change at the SNV position
#' is not an error.
#' Default: \code{0.999}.
#'
#' @param seqError a single positive \code{numeric} between 0 and 1
#' representing the sequencing error rate. Default: \code{0.001}.
#'
#' @param dfPedProfile a \code{data.frame} with the information about
#' the sample(s).
#' Those are mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis" and "Source". All columns must be in
#' \code{character} strings format. The \code{data.frame}
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
#' @param pathProfileGDS a \code{character} string representing the path to
#' the directory where the GDS Sample files will be created.
#'
#' @param genoSource a \code{character} string with two possible values:
#' 'snp-pileup', 'generic' or 'VCF'. It specifies if the genotype files
#' are generated by snp-pileup (Facets) or are a generic format CSV file
#' with at least those columns:
#' 'Chromosome', 'Position', 'Ref', 'Alt', 'Count', 'File1R' and 'File1A'.
#' The 'Count' is the depth at the specified position;
#' 'FileR' is the depth of the reference allele and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally the file can be a VCF file with at least those genotype
#' fields: GT, AD, DP.
#'
#' @param paramProfileGDS a \code{list} parameters ...
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running.
#'
#' @return The  function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Current directory
#' dataDir <- file.path(tempdir())
#'
#' ## Copy required file into current directory
#' file.copy(from=file.path(system.file("extdata/tests", package="RAIDS"),
#'                     "ex1.txt.gz"), to=dataDir)
#'
#' ## The data.frame containing the information about the study
#' ## The 3 mandatory columns: "study.id", "study.desc", "study.platform"
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' studyDF <- data.frame(study.id = "MYDATA",
#'                         study.desc = "Description",
#'                         study.platform = "PLATFORM",
#'                         stringsAsFactors = FALSE)
#'
#' ## The data.frame containing the information about the samples
#' ## The entries should be strings, not factors (stringsAsFactors=FALSE)
#' samplePED <- data.frame(Name.ID=c("ex1", "ex2"),
#'                     Case.ID=c("Patient_h11", "Patient_h12"),
#'                     Diagnosis=rep("Cancer", 2),
#'                     Sample.Type=rep("Primary Tumor", 2),
#'                     Source=rep("Databank B", 2), stringsAsFactors=FALSE)
#' rownames(samplePED) <- samplePED$Name.ID
#'
#' ## List of SNV positions
#' listPositions <- data.frame(snp.chromosome=c(rep(1, 10)),
#'         snp.position=c(3467333, 3467428, 3469375, 3469387, 3469502, 3469527,
#'         3469737, 3471497, 3471565, 3471618))
#'
#' ## Append genotype information to the Profile GDS file
#' result <- RAIDS:::generateProfileGDS(profileFile=file.path(dataDir, "ex1.txt.gz"),
#'             profileName="ex1", listPos=listPositions,
#'             offset=-1, minCov=10, minProb=0.999, seqError=0.001,
#'             dfPedProfile=samplePED, batch=1, studyDF=studyDF,
#'             pathProfileGDS=dataDir, genoSource="snp-pileup",
#'             verbose=FALSE)
#'
#' ## The function returns OL when successful
#' result
#'
#' ## The Profile GDS file 'ex1.gds' has been created in the
#' ## specified directory
#' list.files(dataDir)
#'
#' ## Unlink Profile GDS file (created for demo purpose)
#' unlink(file.path(dataDir, "ex1.gds"))
#' unlink(file.path(dataDir, "ex1.txt.gz"))
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn write.gdsn openfn.gds
#' @importFrom stats qbinom
#' @importFrom utils read.csv
#' @encoding UTF-8
#' @keywords internal
generateProfileGDS <- function(profileFile, profileName,
                            listPos, offset, minCov=10, minProb=0.999,
                            seqError=0.001, dfPedProfile, batch, studyDF, pathProfileGDS,
                            genoSource, paramProfileGDS, verbose) {

    # File with the description of the SNP keep
    # if(genoSource == "VCF"){
    #     listMat <- dir(pathGeno, pattern = ".+.vcf.gz")
    #     listSampleFile <- gsub(".vcf.gz", "", listMat)
    # }else{
    #     listMat <- dir(pathGeno, pattern = ".+.txt.gz")
    #     listSampleFile <- gsub(".txt.gz", "", listMat)
    # }




    # for(i in seq_len(length(listSamples))) {
    #pos <- which(listSampleFile == listSamples[i])

    if(verbose) { message("generateProfileGDS start ", " ", Sys.time()) }

#    if(length(pos) == 1) {

        #
    if(genoSource == "snp-pileup") {
        matSample <- readSNVPileupFile(profileFile, offset)
    } else if(genoSource == "generic") {
        matSample <- readSNVFileGeneric(profileFile, offset)
    } else if(genoSource == "VCF") {
        # tmpProfile <- gsub(".vcf.gz", "",listMat[pos])
        matSample <- readSNVVCF(profileFile,
                                profileName=profileName, offset)
    } else if(genoSource == "bam"){

        matSample <- readSNVBAM(fileName=profileFile,
                                varSelected=listPos,
                                paramSNVBAM=paramProfileGDS,
                                offset,
                                verbose=verbose)
        # listPos <- do.call(rbind, listPos)
        colnames(listPos)[1:2] <- c("snp.chromosome", "snp.position")

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

    if(verbose) {message("End read ", Sys.time())}
    g <- as.matrix(rep(-1, nrow(listPos)))
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

    if(is.null(pathProfileGDS)){
        stop("pathProfileGDS is NULL in ",
             "generateGDS1KGgenotypeFromSNPPileup\n")
    } else{
        if(! dir.exists(pathProfileGDS)) {
            dir.create(pathProfileGDS)
        }
    }
    fileGDSSample <- file.path(pathProfileGDS,
                               paste0(profileName, ".gds"))
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

    listSampleGDS <- addStudyGDSSample(gdsSample,
                                       pedProfile=dfPedProfile, batch=batch,
                                       listSamples=c(profileName), studyDF=studyDF, verbose=verbose)

    listCount <- table(matAll$count[matAll$count >= minCov])
    cutOffA <-
        data.frame(count=unlist(vapply(as.integer(names(listCount)),
                                       FUN=function(x, minProb, eProb){
                                           return(max(2,qbinom(minProb, x, eProb))) },
                                       FUN.VALUE=numeric(1), minProb=minProb, eProb=2 * seqError)),
                   allele=unlist(vapply(as.integer(names(listCount)),
                                        FUN=function(x, minProb, eProb){
                                            return(max(2,qbinom(minProb, x, eProb))) },
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
        append.gdsn(var.geno, g, check=TRUE)
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
        message("End ", profileName, " ",  Sys.time())
    }

    ## Success
    return(0L)
}


#' @title Add information related to a specific study and specific samples
#' into a GDS Sample file
#'
#' @description This function add entries related to 1) a specific study and
#' 2) specific samples into a GDS Sample file.
#' The study information is appended to the GDS Sample file "study.list" node
#' when the node is already present in the file. Otherwise, the node is
#' created and then, the information is added.
#' The sample information for all selected samples is appended to the GDS
#' Sample file "study.annot" node
#' when the node is already present in the file. Otherwise, the node is
#' created and then, the information is added.
#'
#' @param gdsProfile an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param pedProfile a \code{data.frame} with the sample information. The
#' \code{data.frame} must have the columns:
#' "Name.ID", "Case.ID", "Sample.Type", "Diagnosis" and "Source".
#' The unique sample identifier of the \code{data.frame} is the "Name.ID"
#' column and the row names of the \code{data.frame} must be the "Name.ID"
#' values.
#'
#' @param batch a \code{integer} corresponding the batch associated to the
#' study.
#'
#' @param listSamples a \code{vector} of \code{character} string representing
#' the samples (samples identifiers) that are saved into the GDS. All the
#' samples must be present in the 'pdeDF' \code{data.frame}.
#' If \code{NULL}, all samples present in the \code{dfPedProfile} are used.
#'
#' @param studyDF a \code{data.frame} with at least the 3 columns: "study.id",
#' "study.desc" and "study.platform". The three columns are in character
#' string format (no factor).
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function.
#'
#' @return a \code{vector} of \code{character} strings representing the sample
#' identifiers that have been saved in the GDS Sample file.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Create a temporary GDS file in an current directory
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_11.gds")
#'
#' ## Create and open the GDS file
#' tmpGDS  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create a PED data frame with sample information
#' ped1KG <- data.frame(Name.ID=c("1KG_sample_01", "1KG_sample_02"),
#'         Case.ID=c("1KG_sample_01", "1KG_sample_02"),
#'         Sample.Type=rep("Reference", 2), Diagnosis=rep("Reference", 2),
#'         Source=rep("IGSR", 2), stringsAsFactors=FALSE)
#'
#' ## Create a Study data frame with information about the study
#' ## All samples are associated to the same study
#' studyInfo <- data.frame(study.id="Ref.1KG",
#'         study.desc="Unrelated samples from 1000 Genomes",
#'         study.platform="GRCh38 1000 genotypes",
#'         stringsAsFactors=FALSE)
#'
#' ## Add the sample information to the GDS Sample file
#' ## The information for all samples is added (listSamples=NULL)
#' RAIDS:::addStudyGDSSample(gdsProfile=tmpGDS, pedProfile=ped1KG, batch=1,
#'         listSamples=NULL, studyDF=studyInfo, verbose=FALSE)
#'
#' ## Read study information from GDS Sample file
#' read.gdsn(index.gdsn(node=tmpGDS, path="study.list"))
#'
#' ## Read sample information from GDS Sample file
#' read.gdsn(index.gdsn(node=tmpGDS, path="study.annot"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=tmpGDS)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn append.gdsn
#' @encoding UTF-8
#' @keywords internal
addStudyGDSSample <- function(gdsProfile, pedProfile, batch, listSamples,
                                studyDF, verbose) {

    ## Used only the selected samples (all when listSamples == NULL)
    if(!(is.null(listSamples))) {
        if(length(listSamples) == length(intersect(listSamples,
                                                   pedProfile$Name.ID))) {
            # if we remove the names we should manage the listSamples order
            # something like

            tmp <- order(as.character(listSamples))
            pedProfile <- pedProfile[which(pedProfile$Name.ID %in% listSamples), ]
            pedProfile <- pedProfile[order(pedProfile$Name.ID), ][order(tmp),]

        } else {
            stop("List of samples includes samples not present in ",
                "the \'pedProfile\' data frame. The sample names must be ",
                "present in the \'Name.ID\' column. The row names should",
                " be assigned the \'Name.ID\'.")
        }
    } else {
        listSamples <- pedProfile$Name.ID
    }

    ## Create the study data frame that is going to be saved
    df <- data.frame(study.id=studyDF$study.id, study.desc=studyDF$study.desc,
        study.platform=studyDF$study.platform, stringsAsFactors=FALSE)

    ## Append study information to "study.list" when node already present
    ## Otherwise, create node and add study information into it
    if(! "study.list" %in% ls.gdsn(gdsProfile)) {
        ## Create study node and add study information into GDS Sample file
        add.gdsn(gdsProfile, "study.list", df)

        ## Create data frame containing sample information and add it to GDS
        study.annot <- data.frame(data.id=pedProfile[, "Name.ID"],
            case.id=pedProfile[, "Case.ID"],
            sample.type=pedProfile[, "Sample.Type"],
            diagnosis=pedProfile[, "Diagnosis"],
            source=pedProfile[, "Source"],
            study.id=rep(studyDF$study.id, nrow(pedProfile)),
            stringsAsFactors=FALSE)
        add.gdsn(gdsProfile, "study.annot", study.annot)

        if(verbose) { message("study.annot DONE ", Sys.time()) }
    } else{
        ## Append study information to existing node
        append.gdsn(index.gdsn(gdsProfile, "study.list/study.id"),
                        df$study.id, check=TRUE)
        append.gdsn(index.gdsn(gdsProfile, "study.list/study.desc"),
                        df$study.desc, check=TRUE)
        append.gdsn(index.gdsn(gdsProfile, "study.list/study.platform"),
                        df$study.platform, check=TRUE)

        ## Create data frame containing sample information
        study.annot <- data.frame(data.id=pedProfile[, "Name.ID"],
            case.id=pedProfile[, "Case.ID"],
            sample.type=pedProfile[, "Sample.Type"],
            diagnosis=pedProfile[, "Diagnosis"],
            source=pedProfile[, "Source"],
            study.id=rep(studyDF$study.id, nrow(pedProfile)),
            stringsAsFactors=FALSE)

        ## Append sample information to existing node
        append.gdsn(index.gdsn(gdsProfile, "study.annot/data.id"),
                        study.annot$data.id, check=TRUE)
        append.gdsn(index.gdsn(gdsProfile, "study.annot/case.id"),
                        study.annot$case.id, check=TRUE)
        append.gdsn(index.gdsn(gdsProfile, "study.annot/sample.type"),
                        study.annot$sample.type, check=TRUE)
        append.gdsn(index.gdsn(gdsProfile, "study.annot/diagnosis"),
                        study.annot$diagnosis, check=TRUE)
        append.gdsn(index.gdsn(gdsProfile, "study.annot/source"),
                        study.annot$source, check=TRUE)
        append.gdsn(index.gdsn(gdsProfile, "study.annot/study.id"),
                        study.annot$study.id, check=TRUE)

        if(verbose) { message("study.annot DONE ", Sys.time()) }
    }

    ## Return the vector of profile IDs that have been added to Profile GDS file
    return(pedProfile[,"Name.ID"])
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
#' \describe{
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
#' are: "corr", "r", "dprime" and "composite".
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
#' ## Required
#' library(SNPRelate)
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
runLDPruning <- function(gds, method,
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
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param dfPedReference a \code{data.frame} with the information about
#' the sample(s).
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
#' ## Required library
#' library(gdsfmt)
#'
#' ## Create a temporary GDS file in an test directory
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_03.gds")
#'
#' ## Create and open the GDS file
#' tmpGDS <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create "sample.id" node (the node must be present)
#' add.gdsn(node=tmpGDS, name="sample.id", val=c("sample_01",
#'         "sample_02"))
#'
#' ## Create "sample.annot" node (the node must be present)
#' add.gdsn(node=tmpGDS, name="sample.annot", val=data.frame(
#'             Name.ID=c("sample_01", "sample_02"),
#'             sex=c(1,1),  # 1:Male  2: Female
#'             pop.group=c("ACB", "ACB"),
#'             superPop=c("AFR", "AFR"),
#'             batch=c(1, 1),
#'             stringsAsFactors=FALSE))
#'
#' sync.gds(gdsfile=tmpGDS)
#'
#' ## Create a data.frame with information about samples
#' sample_info <- data.frame(Name.ID=c("sample_04", "sample_05",
#'                                 "sample_06"),
#'                         sex=c(1,2,1),  # 1:Male  2: Female
#'                         pop.group=c("ACB", "ACB", "ACB"),
#'                         superPop=c("AFR", "AFR", "AFR"),
#'                         stringsAsFactors=FALSE)
#'
#' ## The row names must be the sample identifiers
#' rownames(sample_info) <- sample_info$Name.ID
#'
#' ## Add information about 2 samples to the GDS file
#' RAIDS:::appendGDSRefSample(gdsReference=tmpGDS,
#'         dfPedReference=sample_info,
#'         batch=2, listSamples=c("sample_04", "sample_06"), verbose=FALSE)
#'
#' ## Read sample identifier list
#' ## Only "sample_04" and "sample_06" should have been added
#' read.gdsn(index.gdsn(node=tmpGDS, path="sample.id"))
#'
#' ## Read sample information from GDS file
#' ## Only "sample_04" and "sample_06" should have been added
#' read.gdsn(index.gdsn(node=tmpGDS, path="sample.annot"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=tmpGDS)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn append.gdsn
#' @encoding UTF-8
#' @keywords internal
appendGDSRefSample <- function(gdsReference, dfPedReference, batch=1,
                                listSamples=NULL, verbose=TRUE) {

    ## Only keep selected samples
    if(!(is.null(listSamples))){
        dfPedReference <- dfPedReference[listSamples,]
    }

    ## Append sample identifiers to the "sample.id" node
    sampleGDS <- index.gdsn(gdsReference, "sample.id")
    append.gdsn(sampleGDS, val=dfPedReference$Name.ID, check=TRUE)

    ## Create the data.frame with the sample information
    samp.annot <- data.frame(sex = dfPedReference[, "sex"],
                                pop.group=dfPedReference[, "pop.group"],
                                superPop=dfPedReference[, "superPop"],
                                batch=rep(batch, nrow(dfPedReference)),
                                stringsAsFactors=FALSE)

    if(verbose) { message("Annot") }

    ## Append data.frame to "sample.annot" node
    curAnnot <- index.gdsn(gdsReference, "sample.annot/sex")
    append.gdsn(curAnnot, samp.annot$sex, check=TRUE)
    curAnnot <- index.gdsn(gdsReference, "sample.annot/pop.group")
    append.gdsn(curAnnot, samp.annot$pop.group, check=TRUE)
    curAnnot <- index.gdsn(gdsReference, "sample.annot/superPop")
    append.gdsn(curAnnot, samp.annot$superPop, check=TRUE)
    curAnnot <- index.gdsn(gdsReference, "sample.annot/batch")
    append.gdsn(curAnnot, samp.annot$batch, check=TRUE)

    if(verbose) { message("Annot done") }

    return(0L)
}


#' @title Add the pruned.study entry related to the SNV dataset in the
#' Profile GDS file
#'
#' @description This function adds the names of the SNVs into the node called
#' "pruned.study" in GDS
#' Sample file. If a "pruned.study" entry is already present, the entry is
#' deleted and a new entry is created.
#'
#' @param gdsProfile an object of class \link[gdsfmt]{gds.class} (a GDS file),
#' the opened Profile GDS file.
#'
#' @param pruned a \code{vector} of \code{character} string representing the
#' name of the SNVs.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Create a temporary GDS file in an test directory
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_1.gds")
#'
#' ## Create and open the GDS file
#' tmpGDS  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Vector of low allelic fraction
#' study <- c("s19222", 's19588', 's19988', 's20588', 's23598')
#'
#' ## Add segments to the GDS file
#' RAIDS:::addGDSStudyPruning(gdsProfile=tmpGDS, pruned=study)
#'
#' ## Read lap information from GDS file
#' read.gdsn(index.gdsn(node=tmpGDS, path="pruned.study"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=tmpGDS)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @keywords internal
addGDSStudyPruning <- function(gdsProfile, pruned) {

    ## Delete the pruned.study entry if present in the Profile GDS file
    if("pruned.study" %in% ls.gdsn(gdsProfile)) {
            delete.gdsn(index.gdsn(node=gdsProfile, "pruned.study"))
    }

    ## Create the pruned.study node in the Profile GDS file
    varPruned <- add.gdsn(node=gdsProfile, name="pruned.study", val=pruned)

    # Write the data cached in memory to the Profile GDS file
    sync.gds(gdsfile=gdsProfile)

    return(0L)
}


#' @title Add information related to low allelic fraction associated to
#' the SNV dataset for a specific sample into a GDS file
#'
#' @description The function adds the information related to low allelic
#' fraction
#' associated to the SNV dataset for a specific sample into a
#' GDS file, more specifically, in the "lap" node. The "lap" node must
#' already be present in the GDS file.
#'
#' @param gdsProfile an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), a GDS file.
#'
#' @param snpLap a \code{vector} of \code{numeric} value representing the
#' low allelic fraction for each SNV present in the SNV dataset. The
#' values should be between \code{0} and \code{0.50}. The
#' length of the \code{vector} should correspond to the number of SNVs
#' present in the "snp.id" entry of the GDS sample file.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Create a temporary GDS file
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP.gds")
#'
#' ## Create and open the GDS file
#' gdsFile  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create a "lap" node
#' add.gdsn(node=gdsFile, name="lap", val=rep(10L, 12))
#' sync.gds(gdsFile)
#'
#' ## Vector of low allelic fraction
#' lap <- c(0.1, 0.23, 0.34, 0.00, 0.12, 0.11, 0.33, 0.55)
#'
#' ## Add segments to the GDS file
#' RAIDS:::addUpdateLap(gdsProfile=gdsFile, snpLap=lap)
#'
#' ## Read lap information from GDS file
#' read.gdsn(index.gdsn(node=gdsFile, path="lap"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=gdsFile)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @keywords internal
addUpdateLap <- function(gdsProfile, snpLap) {

    snpLap2 <- write.gdsn(index.gdsn(gdsProfile, "lap"), snpLap)

    sync.gds(gdsProfile)

    return(0L)
}


#' @title Extract the block identifiers for a list of SNVs
#'
#' @description The function uses the GDS Reference Annotation file to extract
#' the unique block identifiers for a list of SNVs. The block type that is
#' going to be used to extract the information has to be provided by the
#' user.
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the opened Reference SNV Annotation GDS file.
#'
#' @param snpIndex a \code{vectcor} of \code{integer} representing the
#' indexes of the SNVs of interest.
#'
#' @param blockTypeID a \code{character} string corresponding to the block
#' type used to extract the block identifiers. The block type must be
#' present in the GDS Reference Annotation file.
#'
#' @return a \code{vector} of \code{numeric} corresponding to the
#' block identifiers for the SNVs of interest.
#'
#' @examples
#'
#' # Required library
#' library(gdsfmt)
#'
#' ## Path to the demo 1KG Annotation GDS file located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' path1KG <- file.path(dataDir, "tests")
#' fileAnnotGDS <- file.path(path1KG, "ex1_good_small_1KG_Annot.gds")
#'
#' gdsRefAnnotation <- openfn.gds(fileAnnotGDS)
#'
#' ## The indexes for the SNVs of interest
#' snpIndex <- c(1,3,5,6,9)
#'
#' ## Extract the block identifiers for the SNVs represented by their indexes
#' ## for the block created using the genes from Hsapiens Ensembl v86
#' RAIDS:::getBlockIDs(gdsRefAnnot=gdsRefAnnotation, snpIndex=snpIndex,
#'                         blockTypeID="GeneS.Ensembl.Hsapiens.v86")
#'
#' closefn.gds(gdsRefAnnotation)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @keywords internal
getBlockIDs <- function(gdsRefAnnot, snpIndex, blockTypeID) {

    block.annot <- read.gdsn(index.gdsn(gdsRefAnnot, "block.annot"))
    pos <- which(block.annot$block.id == blockTypeID)

    if(length(pos) != 1) {
        stop("The following block type is not found in the ",
                "GDS Annotation file: \'", blockTypeID, "\'")
    }

    b <- read.gdsn(index.gdsn(gdsRefAnnot, "block"), start=c(1, pos),
                    count = c(-1, 1))[snpIndex]

    return(b)
}


#' @title Add information related to segments associated to the SNV
#' dataset for a specific sample into a GDS file
#'
#' @description The function adds the information related to segments
#' associated to the SNV dataset for a specific sample into a
#' GDS file, more specifically, in the "segment" node. If the "segment" node
#' already exists, the previous information is erased.
#'
#' @param gdsProfile an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), a GDS Sample file.
#'
#' @param snpSeg a \code{vector} of \code{integer} representing the segment
#' identifiers associated to each SNV selected for the specific sample. The
#' length of the \code{vector} should correspond to the number of SNVs
#' present in the "snp.id" entry of the GDS sample file.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Temporary GDS file
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Vector of segment identifiers
#' segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L)
#'
#' ## Add segments to the GDS file
#' RAIDS:::addUpdateSegment(gdsProfile=GDS_file_tmp, snpSeg=segments)
#'
#' ## Read segments information from GDS file
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="segment"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @keywords internal
addUpdateSegment <- function(gdsProfile, snpSeg) {

    if("segment" %in% ls.gdsn(gdsProfile)) {
        snpLap <- write.gdsn(index.gdsn(gdsProfile, "segment"), snpSeg)
    } else{
        snpLap <- add.gdsn(gdsProfile, "segment", snpSeg, storage="uint32")
    }

    sync.gds(gdsProfile)

    ## Successful
    return(0L)
}


#' @title Append sample names into a GDS file
#'
#' @description This function append the sample identifiers into the
#' "samples.id" node of a GDS file.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param listSamples a \code{vector} of \code{character} string representing
#' the sample identifiers to be added to GDS file.
#'
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Temporary GDS file in current directory
#' gdsFilePath <- file.path(tempdir(), "GDS_TEMP_04.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create "sample.id" node (the node must be present)
#' add.gdsn(node=GDS_file_tmp, name="sample.id", val=c("sample_01",
#'         "sample_02"))
#'
#' sync.gds(gdsfile=GDS_file_tmp)
#'
#' ## Add information about 2 samples to the GDS file
#' RAIDS:::appendGDSSampleOnly(gds=GDS_file_tmp,
#'         listSamples=c("sample_03", "sample_04"))
#'
#' ## Read sample identifier list
#' ## Only "sample_03" and "sample_04" should have been added
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="sample.id"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn append.gdsn
#' @encoding UTF-8
#' @keywords internal
appendGDSSampleOnly <- function(gds, listSamples) {

    sampleGDS <- index.gdsn(gds, "sample.id")

    append.gdsn(sampleGDS, val=listSamples, check=TRUE)

    return(0L)
}

#' @title Append information associated to ld blocks, as indexes, into the
#' Population Reference SNV Annotation GDS file
#'
#' @description The function appends the information about the ld blocks into
#' the Population Reference SNV Annotation GDS file. The information is
#' extracted from the parameter listBlock.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class}
#' (GDS file), an opened Reference Annotation GDS file.
#'
#' @param listBlock a \code{array} of integer
#' representing the linkage disequilibrium block for
#' each SNV in the in the same order than the variant
#' in Population reference dataset.
#'
#' @param blockName a \code{character} string representing the id of the block.
#' The blockName should not exist in \'gdsRefAnnotFile\'.
#'
#' @param blockDesc a \code{character} string representing the description of
#' the block.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library for GDS
#' library(gdsfmt)
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#  ## Temporary file
#' fileAnnotGDS <- file.path(tempdir(), "ex1_good_small_1KG_Ann_GDS.gds")
#'
#'
#' file.copy(file.path(dataDir, "tests",
#'     "ex1_NoBlockGene.1KG_Annot_GDS.gds"), fileAnnotGDS)
#'
#'
#' fileReferenceGDS  <- file.path(dataDir, "tests",
#'     "ex1_good_small_1KG.gds")
#'  \donttest{
#'     gdsRef <- openfn.gds(fileReferenceGDS)
#'     listBlock <- read.gdsn(index.gdsn(gdsRef, "snp.position"))
#'     listBlock <- rep(-1, length(listBlock))
#'     closefn.gds(gdsRef)
#'     gdsAnnot1KG <- openfn.gds(fileAnnotGDS, readonly=FALSE)
#'     ## Append information associated to blocks
#'     RAIDS:::addGDS1KGLDBlock(gds=gdsAnnot1KG,
#'         listBlock=listBlock,
#'         blockName="blockEmpty",
#'         blockDesc="Example")
#'
#'     closefn.gds(gdsAnnot1KG)
#'
#'     gdsAnnot1KG <- openfn.gds(fileAnnotGDS)
#'     print(gdsAnnot1KG)
#'
#'     closefn.gds(gdsAnnot1KG)
#' }
#'
#' ## Remove temporary file
#' unlink(fileAnnotGDS, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn ls.gdsn compression.gdsn
#' @importFrom gdsfmt append.gdsn sync.gds
#' @encoding UTF-8
#' @keywords internal
addGDS1KGLDBlock <- function(gds, listBlock, blockName, blockDesc) {

    blockAnnot <- data.frame(block.id=blockName,
                             block.desc=blockDesc,
                             stringsAsFactors=FALSE)

    if(! ("block.annot" %in% ls.gdsn(gds))) {
        varBlockAnnot <- add.gdsn(gds, "block.annot", blockAnnot)
    }else {
        curAnnot <- index.gdsn(gds, "block.annot/block.id")
        append.gdsn(curAnnot,blockAnnot$block.id)
        curAnnot <- index.gdsn(gds, "block.annot/block.desc")
        append.gdsn(curAnnot, blockAnnot$block.desc)
    }

    varBlock <- NULL
    if(! ("block" %in% ls.gdsn(gds))){
        varBlock <- add.gdsn(gds, "block",
                             valdim=c(length(listBlock), 1),
                             listBlock, storage="int32",
                             compress = "LZ4_RA")
        readmode.gdsn(varBlock)

    }else {
        if(is.null(varBlock)) {
            varBlock <- index.gdsn(gds, "block")
            varBlock <- compression.gdsn(varBlock, "")
        }
        append.gdsn(varBlock, listBlock)
        varBlock <- compression.gdsn(varBlock, "LZ4_RA")
    }

    sync.gds(gds)

    return(0L)
}

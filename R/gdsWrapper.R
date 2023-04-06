#' @title Append sample names into a GDS file
#'
#' @description This function append the sample identifiers into the
#' "samples.id" node of a GDS file.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param listSample a \code{vector} of \code{character} string representing
#' the sample identifiers to be added to GDS file.
#'
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP_04.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create "sample.id" node (the node must be present)
#' add.gdsn(node=GDS_file_tmp, name="sample.id", val=c("sample_01",
#'     "sample_02"))
#'
#' sync.gds(gdsfile=GDS_file_tmp)
#'
#' ## Add information about 2 samples to the GDS file
#' RAIDS:::appendGDSSampleOnly(gds=GDS_file_tmp,
#'     listSamples=c("sample_03", "sample_04"))
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
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened GDS file.
#'
#' @param pedDF a \code{data.frame} with the sample information. The
#' \code{data.frame} must have the columns:
#' "Name.ID", "Case.ID", "Sample.Type", "Diagnosis" and "Source".
#' The unique sample identifier of the \code{data.frame} is the "Name.ID"
#' column and the row names of the \code{data.frame} must be the "Name.ID"
#' values.
#'
#' @param batch a \code{integer} corresponding the batch associated to the
#' study. Default: \code{1}.
#'
#' @param listSamples a \code{vector} of \code{character} string representing
#' the samples (samples identifiers) that are saved into the GDS. All the
#' samples must be present in the 'pdeDF' \code{data.frame}.
#' If \code{NULL}, all samples present in the \code{pedDF} are used.
#' Default: \code{NULL}.
#'
#' @param studyDF a \code{data.frame} with at least the 3 columns: "study.id",
#' "study.desc" and "study.platform". The three columns are in character
#' string format (no factor).
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{TRUE}.
#'
#' @return a \code{vector} of \code{character} strings representing the sample
#' identifiers that have been saved in the GDS Sample file.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP_11.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create a PED data frame with sample information
#' ped1KG <- data.frame(Name.ID=c("1KG_sample_01", "1KG_sample_02"),
#'     Case.ID=c("1KG_sample_01", "1KG_sample_02"),
#'     Sample.Type=rep("Reference", 2), Diagnosis=rep("Reference", 2),
#'     Source=rep("IGSR", 2), stringsAsFactors=FALSE)
#'
#' ## Create a Study data frame with information about the study
#' ## All samples are associated to the same study
#' studyInfo <- data.frame(study.id="Ref.1KG",
#'     study.desc="Unrelated samples from 1000 Genomes",
#'     study.platform="GRCh38 1000 genotypes",
#'     stringsAsFactors=FALSE)
#'
#' ## Add the sample information to the GDS Sample file
#' ## The information for all samples is added (listSamples=NULL)
#' RAIDS:::addStudyGDSSample(gds=GDS_file_tmp, pedDF=ped1KG, batch=1,
#'     listSamples=NULL, studyDF=studyInfo, verbose=FALSE)
#'
#' ## Read study information from GDS Sample file
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="study.list"))
#'
#' ## Read sample information from GDS Sample file
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="study.annot"))
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
addStudyGDSSample <- function(gds, pedDF, batch=1, listSamples=NULL, studyDF,
                                verbose=TRUE) {

    ## This validation is not necessary as the function is internal
    if(sum(c("study.id", "study.desc", "study.platform") %in%
                colnames(studyDF)) != 3 ) {
        stop("The \'studyDF\' data frame is incomplete. ",
                "One or more mandatory column is missing.\n")
    }

    ## Used only the selected samples (all when listSamples == NULL)
    if(!(is.null(listSamples))) {
        if(length(listSamples) == length(intersect(listSamples,
                                                    rownames(pedDF)))) {
            pedDF <- pedDF[listSamples,]
        } else {
            stop("List of samples includes samples not present in ",
                    "the \'pedDF\' data frame. The sample names must be ",
                    "present in the \'Name.ID\' column. The row names should",
                    " be assigned the \'Name.ID\'.")
        }
    } else {
        listSamples <- pedDF$Name.ID
    }

    ## Create the study data frame that is going to be saved
    df <- data.frame(study.id=studyDF$study.id,
                        study.desc=studyDF$study.desc,
                        study.platform=studyDF$study.platform,
                        stringsAsFactors=FALSE)

    ## Append study information to "study.list" when node already present
    ## Otherwise, create node and add study information into it
    if(! "study.list" %in% ls.gdsn(gds)) {

        ## Create study node and add study information into GDS Sample file
        add.gdsn(gds, "study.list", df)

        ## Create data frame containing sample information
        study.annot <- data.frame(data.id=pedDF[, "Name.ID"],
                                case.id=pedDF[, "Case.ID"],
                                sample.type=pedDF[, "Sample.Type"],
                                diagnosis=pedDF[, "Diagnosis"],
                                source=pedDF[, "Source"],
                                study.id=rep(studyDF$study.id, nrow(pedDF)),
                                stringsAsFactors=FALSE)

        ## Create node and add sample information
        add.gdsn(gds, "study.annot", study.annot)

        if(verbose) { message("study.annot DONE ", Sys.time()) }
    } else{
        ## Append study information to existing node
        append.gdsn(index.gdsn(gds, "study.list/study.id"),
                        df$study.id, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.list/study.desc"),
                        df$study.desc, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.list/study.platform"),
                        df$study.platform, check=TRUE)

        ## Create data frame containing sample information
        study.annot <- data.frame(data.id=pedDF[, "Name.ID"],
                                    case.id=pedDF[, "Case.ID"],
                                    sample.type=pedDF[, "Sample.Type"],
                                    diagnosis=pedDF[, "Diagnosis"],
                                    source=pedDF[, "Source"],
                                    study.id=rep(studyDF$study.id, nrow(pedDF)),
                                    stringsAsFactors=FALSE)

        ## Append sample information to existing node
        append.gdsn(index.gdsn(gds, "study.annot/data.id"),
                        study.annot$data.id, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.annot/case.id"),
                        study.annot$case.id, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.annot/sample.type"),
                        study.annot$sample.type, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.annot/diagnosis"),
                        study.annot$diagnosis, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.annot/source"),
                        study.annot$source, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.annot/study.id"),
                        study.annot$study.id, check=TRUE)

        if(verbose) { message("study.annot DONE ", Sys.time()) }
    }

    ## Return the vector of sample identifiers that have been added
    ## to the GDS Sample file
    return(pedDF[,"Name.ID"])
}


#' @title Add information related to SNVs into a Reference GDS file
#'
#' @description the function adds the SNV information into a Reference
#' GDS file.
#'
#' @param an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Reference GDS file.
#'
#' @param fileFREQ a \code{character} string representing the path and file
#' name of the RDS file with the filtered SNP information.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required package
#' library(withr)
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## Temporary Reference GDS file
#' file1KG <- local_file(file.path(data.dir, "1KG_TEMP_002.gds"))
#' filenewGDS <- createfn.gds(file1KG)
#'
#' ## The RDS file containing the filtered SNP information
#' fileFilerterSNVs <- file.path(data.dir, "mapSNVSelected_Demo.rds")
#'
#' ## Add SNV information to Reference GDS
#' RAIDS:::generateGDSSNPinfo(gds=filenewGDS, fileFREQ=fileFilerterSNVs,
#'     verbose=TRUE)
#'
#' ## Close GDS file (important)
#' closefn.gds(filenewGDS)
#'
#' ## Remove temporary files
#' deferred_run()
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
generateGDSSNPinfo <- function(gds, fileFREQ, verbose) {

    mapSNVSel <- readRDS(file=fileFREQ)

    if(verbose) { message("Read mapSNVSel DONE ", Sys.time()) }

    add.gdsn(node=gds, name="snp.id", paste0("s",seq_len(nrow(mapSNVSel))))

    if(verbose) { message("SNP part snp.id DONE ", Sys.time()) }

    add.gdsn(node=gds, name="snp.chromosome",
            as.integer(gsub("chr", "", mapSNVSel$CHROM)), storage = "uint16")

    if(verbose) { message("SNP part snp.chromosome DONE ", Sys.time()) }

    add.gdsn(node=gds, name="snp.position", as.integer(mapSNVSel$POS),
                storage="int32")

    if(verbose) { message("SNP part snp.position DONE ", Sys.time()) }

    add.gdsn(node=gds, name="snp.allele",
                paste0(mapSNVSel$REF, "/", mapSNVSel$ALT))

    if(verbose) { message("SNP part 1 DONE ", Sys.time()) }

    add.gdsn(node=gds, name="snp.AF", as.numeric(mapSNVSel$AF),
                    storage="packedreal24")

    if(verbose) { message("SNP part AF DONE ", Sys.time()) }

    add.gdsn(node=gds, name="snp.EAS_AF", val=as.numeric(mapSNVSel$EAS_AF),
                    storage="packedreal24")
    add.gdsn(node=gds, name="snp.EUR_AF", val=as.numeric(mapSNVSel$EUR_AF),
                    storage="packedreal24")
    add.gdsn(node=gds, name="snp.AFR_AF", val=as.numeric(mapSNVSel$AFR_AF),
                    storage="packedreal24")
    add.gdsn(node=gds, name="snp.AMR_AF", val=as.numeric(mapSNVSel$AMR_AF),
                    storage="packedreal24")
    add.gdsn(node=gds, name="snp.SAS_AF", val=as.numeric(mapSNVSel$SAS_AF),
                    storage="packedreal24")

    ## Return OL when successful
    return(0L)
}



#' @title Add information related to profile genotype into a Reference GDS file
#'
#' @description This function adds the genotype fields with the associated
#' information into the Reference GDS file for the selected profiles.
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened Reference GDS file.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the 1K genotyping files for each sample are located. The name of the
#' genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#'
#' @param fileLSNP a \code{character} string representing the path and file
#' name of the RDS file that contains the indexes of the retained SNPs. The
#' file must exist. The file must be a RDS file.
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to profiles (must be the profile identifiers) that will be
#' retained and added to the Reference GDS file.
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn write.gdsn
#' @importFrom utils read.csv2
#' @encoding UTF-8
#' @keywords internal
generateGDSgenotype <- function(gds, pathGeno, fileLSNP, listSamples,
                                    verbose) {

    # File with the description of the SNP keep
    listMat1k <- dir(pathGeno, pattern=".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMat1k)

    listSNP <- readRDS(fileLSNP)

    for(i in seq_len(length(listSamples))) {
        pos <- which(listSample1k == listSamples[i])

        if(verbose) { message(listSamples[i]) }

        if( length(pos) == 1) {
            matSample <- read.csv2(file.path(pathGeno, listMat1k[pos]),
                                        row.names=NULL)
            matSample <- matSample[listSNP,, drop=FALSE]
            if(i == 1) {
                var.geno <- add.gdsn(gds, "genotype",
                                valdim=c(nrow(matSample),
                                            length(listSamples)),
                                            storage="bit2")
            }

            # Easier to read
            matSample[matSample[,1] == "0|0", 1] <- 0
            matSample[matSample[,1] == "0|1" | matSample[,1] == "1|0", 1] <- 1
            matSample[matSample[,1] == "1|1", 1] <- 2

            g <- as.matrix(matSample)[,1]

            write.gdsn(var.geno, g, start=c(1, i), count=c(-1, 1))

            rm(matSample)

            if(verbose) { message(listMat1k[pos], " ", i) }
        }else{
            stop("Missing samples genotype in ", listSamples[i])
        }
    }

    return(0L)
}

#' @title This function append the field genotype in the gds file
#'
#' @description TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param pathGeno TODO a PATH to a directory with the a file for each
#' samples with the genotype.
#'
#' @param fileLSNP TODO list of SNP to keep in the file genotype
#'
#' @param listSamples  a \code{array} with the sample to keep
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running. Default: \code{FALSE}.
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils read.csv2
#' @encoding UTF-8
#' @keywords internal
appendGDSgenotype <- function(gds, listSample, pathGeno, fileLSNP,
                                verbose=FALSE) {

    # File with the description of the SNP keep
    listMat1k <- dir(pathGeno, pattern = ".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMat1k)

    listSNP <- readRDS(file=fileLSNP)
    geno.var <- index.gdsn(gds, "genotype")
    g <- read.gdsn(node=geno.var, start=c(1, 1), count=c(1,-1))
    nbSample <- length(g)
    if(verbose) { message(nbSample) }
    for(i in seq_len(length(listSample))) {
        pos <- which(listSample1k == listSample[i])
        if( length(pos) == 1) {
            matSample <- read.csv2(file.path(pathGeno, listMat1k[pos]),
                                        row.names = NULL)
            matSample <- matSample[listSNP,, drop=FALSE]


            # Not faster but harder to read
            # matSample[,1] <- rowSums(t(matrix(as.numeric(unlist(strsplit(
            #                      matSample[,1], "\\|"))),nr=2)))
            # Easier to read
            matSample[matSample[,1] == "0|0",1] <- 0
            matSample[matSample[,1] == "0|1" | matSample[,1] == "1|0",1] <- 1
            matSample[matSample[,1] == "1|1",1] <- 2

            g <- as.matrix(matSample)[,1]
            append.gdsn(geno.var,g, check=TRUE)

            rm(matSample)
            if(verbose) { message(listMat1k[pos], " ", i) }
        }else {
            stop("Missing 1k samples ", listSample[i])
        }
    }

    return(0L)
}


#' @title create a file tfam file for plink from the gds file
#'
#' @description TODO
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file.
#'
#' @param listSample a \code{array} with the sample to keep TODO
#'
#' @param pedOUT TODO a PATH and file name to the output file
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils write.table
#' @encoding UTF-8
#' @keywords export
gds2tfam <- function(gds, listSample, pedOUT) {

    sampleGDS <- index.gdsn(gds, "sample.id")
    sampleId <-read.gdsn(sampleGDS)
    listS <- which(sampleId %in% listSample)

    sampleGDS <- index.gdsn(gds, "sample.annot")
    sampleANNO <-read.gdsn(sampleGDS)

    pedFile <- data.frame(famId=paste0("F", seq_len(length(listSample))),
                            id=sampleId[listS],
                            fa=rep("0",length(listSample)),
                            mo=rep("0",length(listSample)),
                            sex=sampleANNO$sex[listS],
                            pheno=rep(1,length(listSample)),
                            stringsAsFactors=FALSE)

    write.table(pedFile, pedOUT,
                    quote=FALSE, sep="\t",
                    row.names=FALSE,
                    col.names=FALSE)

}

#' @title create a file tfam file for plink from the gds file
#'
#' @description TODO
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1 KG GDS file.
#'
#' @param listSample  a \code{array} with the sample to keep
#'
#' @param sampleANNO a \code{data.frame} with at least column sex and the name
#' must be sample.id
#'
#' @param pedOUT TODO a PATH and file name to the output file
#'
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils write.table
#' @encoding UTF-8
#' @keywords export
gds2tfamSample <- function(gds, listSample, sampleANNO, pedOUT) {

    sampleGDS <- index.gdsn(gds, "sample.id")
    sampleId <-read.gdsn(sampleGDS)
    listS <- which(sampleId %in% listSample)

    pedFile <- data.frame(famId=paste0("F", seq_len(length(listSample))),
                            id=sampleId[listS],
                            fa=rep("0",length(listSample)),
                            mo=rep("0",length(listSample)),
                            sex=sampleANNO[sampleId[listS], "sex"],
                            pheno=rep(1,length(listSample)),
                            stringsAsFactors=FALSE)

    write.table(pedFile, pedOUT, quote=FALSE, sep="\t",
                    row.names=FALSE, col.names=FALSE)

}


#' @title create a file tped file for plink from the gds file
#'
#' @description TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param listSample  a \code{array} with the sample to keep
#'
#' @param listSNP  a \code{array} with the snp.id to keep
#'
#' @param pedOUT TODO a PATH and file name to the output file
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @importFrom utils write.table
#' @encoding UTF-8
#' @keywords export
gds2tped <- function(gds, listSample, listSNP, pedOUT) {

    sampleGDS <- index.gdsn(gds, "sample.id")
    sampleId <-read.gdsn(sampleGDS)
    listS <- which(sampleId %in% listSample)

    snpGDS <- index.gdsn(gds, "snp.id")
    snpId <- read.gdsn(snpGDS)
    listKeep <- which(snpId %in% listSNP)
    snpId <- snpId[listKeep]

    snpGDS <- index.gdsn(gds, "snp.chromosome")
    snpChr <- read.gdsn(snpGDS)
    snpChr <- snpChr[listKeep]

    snpGDS <- index.gdsn(gds, "snp.position")
    snpPos <- read.gdsn(snpGDS)
    snpPos <- snpPos[listKeep]

    tped <- list()
    tped[[1]] <- snpChr
    tped[[2]] <- snpId
    tped[[3]] <- rep(0,length(snpId))
    tped[[4]] <- snpPos
    k<-4
    geno.var <- index.gdsn(gds, "genotype")
    for(i in listS){
        k <- k + 1

        tmp <- read.gdsn(geno.var, start=c(1, i), count=c(-1,1))[listKeep]

        # 0 1 0 1 0 1
        tped[[k]] <- (tmp == 2) + 1
        k <- k + 1
        tped[[k]] <- (tmp > 0) + 1

    }

    write.table(tped, pedOUT, quote=FALSE, sep="\t", row.names=FALSE,
                        col.names=FALSE)

}


#' @title Add the pruned.study entry related to the SNV dataset in the
#' Profile GDS file
#'
#' @description This function adds the names of the SNVs into the node called
#' "pruned.study" in GDS
#' Sample file. If a "pruned.study" entry is already present, the entry is
#' deleted and a new entry is created.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' opened Profile GDS file.
#'
#' @param pruned a \code{vector} of \code{character} string representing the
#' name of the SNVs.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' #' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP_1.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Vector of low allelic fraction
#' study <- c("s19222", 's19588', 's19988', 's20588', 's23598')
#'
#' ## Add segments to the GDS file
#' RAIDS:::addGDSStudyPruning(gds=GDS_file_tmp, pruned=study)
#'
#' ## Read lap information from GDS file
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="pruned.study"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @keywords internal
addGDSStudyPruning <- function(gds, pruned) {

    ## Delete the pruned.study entry if present in the GDS Sample file
    if("pruned.study" %in% ls.gdsn(gds)) {
        delete.gdsn(index.gdsn(node=gds, "pruned.study"))
    }

    ## Create the pruned.study node in the GDS Sample file
    var.Pruned <- add.gdsn(node=gds, name="pruned.study", val=pruned)

    # Write the data cached in memory to the GDS Sample file
    sync.gds(gdsfile=gds)

    return(0L)
}

#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened for the sample
#'
#' @param listBlock TODO
#'
#' @param blockName TODO
#'
#' @param blockDesc TODO
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn ls.gdsn compression.gdsn
#' @importFrom gdsfmt append.gdsn sync.gds
#' @encoding UTF-8
#' @keywords internal
addGDS1KGLDBlock <- function(gds, listBlock, blockName, blockDesc) {

    block.annot <- data.frame(block.id=blockName,
                                block.desc=blockDesc,
                                stringsAsFactors=FALSE)

    if(! ("block.annot" %in% ls.gdsn(gds))) {
        var.block.annot <- add.gdsn(gds, "block.annot", block.annot)
    }else {
        curAnnot <- index.gdsn(gds, "block.annot/block.id")
        append.gdsn(curAnnot,block.annot$block.id)
        curAnnot <- index.gdsn(gds, "block.annot/block.desc")
        append.gdsn(curAnnot, block.annot$block.desc)
    }

    var.block <- NULL
    if(! ("block" %in% ls.gdsn(gds))){
        var.block <- add.gdsn(gds, "block",
                                valdim=c(length(listBlock), 1),
                                listBlock, storage="int32",
                                compress = "LZ4_RA")
        readmode.gdsn(var.block)

    }else {
        if(is.null(var.block)) {
            var.block <- index.gdsn(gds, "block")
            var.block <- compression.gdsn(var.block, "")
        }
        append.gdsn(var.block, listBlock)
        var.block <- compression.gdsn(var.block, "LZ4_RA")
    }

    sync.gds(gds)

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
#' @param gds an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), a GDS file.
#'
#' @param snp.lap a \code{vector} of \code{numeric} value representing the
#' low allelic fraction for each SNV present in the SNV dataset. The
#' values should be between \code{0} and \code{0.50}. The
#' length of the \code{vector} should correspond to the number of SNVs
#' present in the "snp.id" entry of the GDS sample file.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Create a "lap" node
#' lap_initial <- rep(0.0, 8)
#' add.gdsn(node=GDS_file_tmp, name="lap", val=rep(10L, 12))
#' sync.gds(GDS_file_tmp)
#'
#' ## Vector of low allelic fraction
#' lap <- c(0.1, 0.23, 0.34, 0.00, 0.12, 0.11, 0.33, 0.55)
#'
#' ## Add segments to the GDS file
#' RAIDS:::addUpdateLap(gds=GDS_file_tmp, snp.lap=lap)
#'
#' ## Read lap information from GDS file
#' read.gdsn(index.gdsn(node=GDS_file_tmp, path="lap"))
#'
#' ## Close GDS file
#' closefn.gds(gdsfile=GDS_file_tmp)
#'
#' ## Delete the temporary GDS file
#' unlink(x=gdsFilePath, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @keywords internal
addUpdateLap <- function(gds, snp.lap) {

    snpLap <- write.gdsn(index.gdsn(gds, "lap"), snp.lap)

    sync.gds(gds)

    return(0L)
}


#' @title Add information related to segments associated to the SNV
#' dataset for a specific sample into a GDS file
#'
#' @description The function adds the information related to segments
#' associated to the SNV dataset for a specific sample into a
#' GDS file, more specifically, in the "segment" node. If the "segment" node
#' already exists, the previous information is erased.
#'
#' @param gds an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), a GDS Sample file.
#'
#' @param snp.seg a \code{vector} of \code{integer} representing the segment
#' identifiers associated to each SNV selected for the specific sample. The
#' length of the \code{vector} should correspond to the number of SNVs
#' present in the "snp.id" entry of the GDS sample file.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Create a temporary GDS file in an test directory
#' data.dir <- system.file("extdata/tests", package="RAIDS")
#' gdsFilePath <- file.path(data.dir, "GDS_TEMP.gds")
#'
#' ## Create and open the GDS file
#' GDS_file_tmp  <- createfn.gds(filename=gdsFilePath)
#'
#' ## Vector of segment identifiers
#' segments <- c(1L, 1L, 1L, 2L, 2L, 3L, 3L)
#'
#' ## Add segments to the GDS file
#' RAIDS:::addUpdateSegment(gds=GDS_file_tmp, snp.seg=segments)
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
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @keywords internal
addUpdateSegment <- function(gds, snp.seg) {

    if("segment" %in% ls.gdsn(gds)) {
        snpLap <- write.gdsn(index.gdsn(gds, "segment"), snp.seg)
    } else{
        snpLap <- add.gdsn(gds, "segment", snp.seg, storage="uint32")
    }

    sync.gds(gds)

    ## Successful
    return(0L)
}


#' @title Get the block number for each SNV in snp.index
#'
#' @description TODO
#'
#' @param gdsRefAnnot an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the1 1KG SNV Annotation GDS file. RNA specific
#' Default: \code{NULL}.
#'
#' @param snp.index TODO
#'
#' @param block.id a \code{character} string corresponding to the block
#' identifier in \code{gdsRefAnnot}.  RNA specific
#' Default: \code{NULL}
#'
#' @return TODO a \code{vector} of \code{numeric} corresponding to the
#' block number
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @keywords internal
get.Gene.Block <- function(gdsRefAnnot, snp.index, block.id) {

    block.annot <- read.gdsn(index.gdsn(gdsRefAnnot, "block.annot"))
    pos <- which(block.annot$block.id == block.id)

    if(length(pos) != 1) {
        stop("Try to get Gene.Block with block.id problematic ", block.id)
    }

    b <- read.gdsn(index.gdsn(gdsRefAnnot, "block"), start=c(1,pos),
                    count = c(-1,1))[snp.index]

    return(b)
}


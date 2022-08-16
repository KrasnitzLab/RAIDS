#' @title Initialization of the section related to the sample
#' information in the \code{gds} file.
#'
#' @description This function initializes the section related to the sample
#' information in the \code{gds} file. The information is extracted from
#' the \code{data.frame} \code{pedDF} passed to the function.
#'
#' @param gds a \code{gds}.
#'
#' @param pedDF a \code{data.frame} containing the information related to the
#' samples. It must have those columns: "sample.id", "Name.ID", "sex",
#' "pop.group", "superPop" and "batch". The unique id of pedDF
#' is Name.ID and the row.name is Name.ID too.
#'
#' @param listSamples a \code{array} with the sample from pedDF to keep
#'
#' @return the a \code{array} with the sample from pedDF keept
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
generateGDSSample <- function(gds, pedDF, listSamples=NULL){

    if(!(is.null(listSamples))){
        pedDF <- pedDF[listSamples,]
    }

    add.gdsn(gds, "sample.id", pedDF[, "Name.ID"])

    ## Create a data.frame containing the information form the samples
    samp.annot <- data.frame(sex=pedDF[, "sex"],
                                pop.group=pedDF[, "pop.group"],
                                superPop=pedDF[, "superPop"],
                                batch=pedDF[, "batch"],
                                stringsAsFactors=FALSE)

    ## Add the data.frame to the gds object
    add.gdsn(gds, "sample.annot", samp.annot)

    return(pedDF[, "sample.id"])

}

#' @title The function add an array sample.ref to the gds file.It define base
#' on a list of unrelated samples.
#'
#' @description This function create the field sample.ref which is 1 when
#' the samples are a reference and 0 otherwise. The sample.ref is fill based
#' on the file filePart$unrels
#' from  in GENESIS TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param filePart a \code{list} from the function pcairPartition in GENESIS
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
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

    return(0L)
}


#' @title This function append the fields related to the samples. If the
#' samples are part of a study you must uses the addStudyGDSSample
#'
#' @description This function append the fields related to the samples.
#' The fields append are sample.id and the \code{data.frame} sample.annot.
#' If the samples are in the section study the field related to the
#' study must be fill.
#'
#' @param gds a \code{gds}.
#'
#' @param pedDF a \code{data.frame} with the sample info. Must have the column
#' sample.id, Name.ID, sex, pop.group, superPop and batch. The unique id
#' of pedDF is Name.ID and the row.name is Name.ID too.
#'
#' @param listSamples a \code{array} with the sample from pedDF$Name.ID to keep
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{TRUE}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn append.gdsn
#' @encoding UTF-8
#' @keywords internal
appendGDSSample <- function(gds, pedDF, batch=1, listSamples=NULL,
                                verbose=TRUE) {

    if(!(is.null(listSamples))){
        pedDF <- pedDF[listSamples,]
    }

    sampleGDS <- index.gdsn(gds, "sample.id")

    append.gdsn(sampleGDS, val=pedDF$Name.ID, check=TRUE)

    samp.annot <- data.frame(sex = pedDF[, "sex"],
                                pop.group=pedDF[, "pop.group"],
                                superPop=pedDF[, "superPop"],
                                batch=rep(batch, nrow(pedDF)),
                                stringsAsFactors=FALSE)

    if(verbose) { message("Annot") }

    curAnnot <- index.gdsn(gds, "sample.annot/sex")
    append.gdsn(curAnnot,samp.annot$sex, check=TRUE)
    curAnnot <- index.gdsn(gds, "sample.annot/pop.group")
    append.gdsn(curAnnot, samp.annot$pop.group, check=TRUE)
    curAnnot <- index.gdsn(gds, "sample.annot/superPop")
    append.gdsn(curAnnot, samp.annot$superPop, check=TRUE)
    curAnnot <- index.gdsn(gds, "sample.annot/batch")
    append.gdsn(curAnnot, samp.annot$batch, check=TRUE)

    if(verbose) { message("Annot done") }

    return(0L)
}

#' @title This function append the fields sample.id.
#'
#' @description This function append the fields samples.id.
#' The fields append are sample.id with the listSample
#'
#' @param gds a \code{gds}.
#'
#' @param listSample a \code{array} of sample.id to add.
#'
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
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



#' @title This function create the gds file fields related to the study and
#' the sample in it.
#'
#' @description TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param pedDF a \code{data.frame} with the sample info. Must have the column
#' sample.id, Name.ID, sex, pop.group, superPop and batch. The unique id of
#' pedDF is Name.ID and the row.name is Name.ID too.
#'
#' @param batch a \code{integer} corresponding
#'
#' @param listSamples a \code{array} with the sample from pedDF to keep
#'
#' @param studyDF a \code{data.frame} with at least the column study.id,
#' study.desc and study.platform
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{TRUE}.
#'
#' @return TODO
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn append.gdsn
#' @encoding UTF-8
#' @keywords internal
addStudyGDSSample <- function(gds, pedDF, batch=1, listSamples=NULL, studyDF,
                                verbose=TRUE) {

    if(sum(c("study.id", "study.desc", "study.platform") %in%
                colnames(studyDF)) != 3 ) {
        stop("studyDF incomplete in addStudyGDSSample\n")
    }

    if(!(is.null(listSamples))) {
        if(length(listSamples) == length(intersect(listSamples,
                                                    rownames(pedDF)))) {
            pedDF <- pedDF[listSamples,]
        } else {
            stop("List of samples not include in the ped or ped ",
                    "don't have rownames equal to Name.ID\n")
        }
    } else {
        listSamples <- pedDF$Name.ID
    }

    # sampleGDS <- index.gdsn(gds, "sample.id")

    # samplePres <- read.gdsn(sampleGDS)
    # print(paste0("appendGDSSample ", Sys.time()))
    # appendGDSSample(gds, pedDF, batch=batch, listSamples=listSamples)
    # print(paste0("appendGDSSample DONE ", Sys.time()))

    # add.gdsn(gds, "study.offset", length(samplePres))

    df <- data.frame(study.id=studyDF$study.id,
                        study.desc=studyDF$study.desc,
                        study.platform=studyDF$study.platform,
                        stringsAsFactors=FALSE)

    if(! "study.list" %in% ls.gdsn(gds)){
        add.gdsn(gds, "study.list", df)

        study.annot <- data.frame(data.id=pedDF[, "Name.ID"],
                                    case.id=pedDF[, "Case.ID"],
                                    sample.type=pedDF[, "Sample.Type"],
                                    diagnosis=pedDF[, "Diagnosis"],
                                    source=pedDF[, "Source"],
                                    study.id=rep(studyDF$study.id, nrow(pedDF)),
                                    stringsAsFactors=FALSE)

        add.gdsn(gds, "study.annot", study.annot)
        if(verbose) { message("study.annot DONE ", Sys.time()) }
    } else{
        #append.gdsn(gds, "study.list", df)
        append.gdsn(index.gdsn(gds, "study.list/study.id"),
                        df$study.id, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.list/study.desc"),
                        df$study.desc, check=TRUE)
        append.gdsn(index.gdsn(gds, "study.list/study.platform"),
                        df$study.platform, check=TRUE)

        study.annot <- data.frame(data.id=pedDF[, "Name.ID"],
                                  case.id=pedDF[, "Case.ID"],
                                  sample.type=pedDF[, "Sample.Type"],
                                  diagnosis=pedDF[, "Diagnosis"],
                                  source=pedDF[, "Source"],
                                  study.id=rep(studyDF$study.id, nrow(pedDF)),
                                  stringsAsFactors=FALSE)

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

    return(pedDF[,"Name.ID"])
}


#' @title This function create the fields related to the snp
#'
#' @description TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param fileFREQ a \code{character} string with the path and the file
#' naame to a RDS file containing the frequency information
#' TODO describe the file
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{TRUE}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
generateGDSSNPinfo <- function(gds, fileFREQ, verbose=TRUE) {

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



#' @title This function create the field genotype in the gds file
#'
#' @description TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param PATHGENO TODO a PATH to a directory with the a file for each samples
#' with the genotype.
#'
#' @param fileLSNP TODO list of SNP to keep in the file genotype
#'
#' @param listSamples a \code{array} with the sample to keep
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running. Default: \code{FALSE}.
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
generateGDSgenotype <- function(gds, PATHGENO, fileLSNP, listSamples,
                                    verbose=FALSE) {

    # File with the description of the SNP keep
    listMat1k <- dir(PATHGENO, pattern=".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMat1k)

    listSNP <- readRDS(fileLSNP)

    for(i in seq_len(length(listSamples))) {
        pos <- which(listSample1k == listSamples[i])

        if(verbose) { message(listSamples[i]) }

        if( length(pos) == 1) {
            matSample <- read.csv2(file.path(PATHGENO, listMat1k[pos]),
                                        row.names=NULL)
            matSample <- matSample[listSNP,, drop=FALSE]
            if(i == 1) {
                var.geno <- add.gdsn(gds, "genotype",
                                valdim=c(nrow(matSample),
                                            length(listSamples)),
                                            storage="bit2")
            }

            # Not faster but harder to read
            # matSample[,1] <- rowSums(t(matrix(as.numeric(unlist(strsplit(
            #                        matSample[,1], "\\|"))),nr=2)))
            # Easier to read
            matSample[matSample[,1] == "0|0",1] <- 0
            matSample[matSample[,1] == "0|1" | matSample[,1] == "1|0",1] <- 1
            matSample[matSample[,1] == "1|1",1] <- 2

            g <- as.matrix(matSample)[,1]

            write.gdsn(var.geno, g, start=c(1, i), count=c(-1,1))

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
#' @param PATHGENO TODO a PATH to a directory with the a file for each
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
appendGDSgenotype <- function(gds, listSample, PATHGENO, fileLSNP,
                                verbose=FALSE) {

    # File with the description of the SNP keep
    listMat1k <- dir(PATHGENO, pattern = ".+.csv.bz2")
    listSample1k <- gsub(".csv.bz2", "", listMat1k)

    listSNP <- readRDS(file=fileLSNP)
    geno.var <- index.gdsn(gds, "genotype")
    g <- read.gdsn(node=geno.var, start=c(1, 1), count=c(1,-1))
    nbSample <- length(g)
    if(verbose) { message(nbSample) }
    for(i in seq_len(length(listSample))) {
        pos <- which(listSample1k == listSample[i])
        if( length(pos) == 1) {
            matSample <- read.csv2(file.path(PATHGENO, listMat1k[pos]),
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

#' @title This function append the field genotype in the gds file
#'
#' @description TODO
#'
#' @param gds a \code{gds} object.
#'
#' @param matG  a \code{matrix} with the genotype
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
#' @param PATHGENO TODO a PATH to a directory with the a file for each
#' samples with the genotype.
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to the sample identifiers that will have a GDS Sample file created. The
#' sample identifiers must be present in the "Name.ID" column of the
#' \code{data.frame} passed to the \code{pedDF} parameter.
#'
#' @param listPos  a \code{array}
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
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running. Default: \code{FALSE}.
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
generateGDS1KGgenotypeFromSNPPileup <- function(PATHGENO,
                                                listSamples,listPos, offset,
                                                minCov=10, minProb=0.999,
                                                seqError=0.001,
                                                pedStudy,
                                                batch,
                                                studyDF,
                                                PATHGDSSAMPLE=NULL,
                                                verbose=FALSE) {

    # File with the description of the SNP keep
    listMat <- dir(PATHGENO, pattern = ".+.txt.gz")
    listSampleFile <- gsub(".txt.gz", "", listMat)

    g <- as.matrix(rep(-1, nrow(listPos)))

    for(i in seq_len(length(listSamples))) {
        pos <- which(listSampleFile == listSamples[i])

        if(verbose) { message(listSamples[i], " ", Sys.time(), " ", i) }

        if(length(pos) == 1){

            matSample <- read.csv(file.path(PATHGENO, listMat[pos]))


            matSample[, "Chromosome"] <- as.integer(gsub("chr", "",
                                                    matSample[, "Chromosome"]))
            matSample[, "Position"] <- matSample[, "Position"] + offset
            matSample[, "count"] <- rowSums(matSample[, c("File1R", "File1A",
                                                    "File1E", "File1D")])

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

            z <- cbind(c(listPos$snp.chromosome,
                            matSample$Chromosome,
                            matSample$Chromosome),
                        c(listPos$snp.position,
                            matSample$Position,
                            matSample$Position),
                        c(rep(1,nrow(listPos)),
                            rep(0,nrow(matSample)),
                            rep(2,nrow(matSample))),
                        c(rep(0,nrow(listPos)),
                            matSample[, "File1R"],
                            -1 * matSample[, "File1R"]),
                        c(rep(0,nrow(listPos)),
                            matSample[, "File1A"],
                            -1 * matSample[, "File1A"]),
                        c(rep(0,nrow(listPos)),
                            matSample[, "count"],
                            -1 * matSample[, "count"]))
            rm(matSample)
            z <- z[order(z[,1], z[,2], z[,3]),]

            matAll <- data.frame(Chromosome=z[z[,3]==1, 1],
                                    Position=z[z[,3]==1, 2],
                                    File1R=cumsum(z[,4])[z[,3]==1],
                                    File1A=cumsum(z[,5])[z[,3]==1],
                                    count=cumsum(z[,6])[z[,3]==1])
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

                # id <- add.gdsn(gdsSample, "sampleStudy",
                #                 listSamples[i])
            }

            if(! "Ref.count" %in% ls.gdsn(gdsSample)){
                var.Ref <- add.gdsn(gdsSample, "Ref.count",
                                        matAll$File1R,
                                        valdim=c( nrow(listPos), 1),
                                        storage="sp.int16")
                var.Alt <- add.gdsn(gdsSample, "Alt.count",
                                        matAll$File1A,
                                        valdim=c( nrow(listPos), 1),
                                        storage="sp.int16")
                var.Count <- add.gdsn(gdsSample, "Total.count",
                                        matAll$count,
                                        valdim=c( nrow(listPos), 1),
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
                                                batch=batch,
                                                listSamples=c(listSamples[i]),
                                                studyDF=studyDF)

                #closefn.gds(gdsSample)



            listCount <- table(matAll$count[matAll$count >= minCov])
            cutOffA <-
                data.frame(count=unlist(vapply(as.integer(names(listCount)),
                FUN=function(x, minProb, eProb){
                    return(max(2,qbinom(minProb, x,eProb))) },
                FUN.VALUE=numeric(1), minProb=minProb, eProb=2 * seqError )),
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
                                     valdim=c(length(g),
                                              1),
                                     g,
                                     storage="bit2", compress = "LZMA_RA.fast")
                readmode.gdsn(var.geno)
            }


            rm(g)
            closefn.gds(gdsfile=gdsSample)
            if(verbose) { message(listMat[pos], " ", i, " ", Sys.time()) }

        }else{
            stop("Missing samples ", listSamples[i])
        }
    }

    return(0L)
}



#' @title create a file tfam file for plink from the gds file
#'
#' @description TODO
#'
#' @param gds an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the 1KG GDS file.
#'
#' @param listSample  a \code{array} with the sample to keep TODO
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

    #tped <- do.call(cbind, tped)

    write.table(tped, pedOUT,
                    quote=FALSE, sep="\t",
                    row.names=FALSE,
                    col.names=FALSE)

}


#' @title Function just wrap snpgdsIBDKING
#'
#' @description TODO
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param sampleId  a \code{vector} of \code{character} strings representing
#' the samples to keep for the analysis. If \code{NULL}, all samples are used.
#' Default: \code{NULL}.
#'
#' @param snp.id  a \code{vector} of \code{character} strings representing
#' the SNPs to keep for the analysis. If \code{NULL}, all SNPs are used.
#' Default: \code{NULL}.
#'
#' @param maf  a single \code{numeric} representing the threshold for the minor
#' allele frequency. Only the SNPs with ">= maf" are retained.
#' Default: \code{0.05}.
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
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom SNPRelate snpgdsIBDKING
#' @encoding UTF-8
#' @keywords internal
runIBDKING <- function(gds, sampleId=NULL, snp.id=NULL, maf=0.05) {

    ## Calculate IBD coefficients by KING method of moment
    ibd.robust <- snpgdsIBDKING(gds, sample.id=sampleId,
                                    snp.id=snp.id,
                                    maf=maf, type="KING-robust")
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
#' @param listKeep a \code{vector} of SNP identifiers specifying selected
#' SNPs; if \code{NULL}, all SNPs are used in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}} function. Default: \code{NULL}.
#'
#' @param slide.max.bp.v a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done in the \code{\link[SNPRelate]{snpgdsLDpruning}}()
#' function.
#' Default: \code{500000L}.
#'
#' @param ld.threshold.v a single \code{numeric} value that represents the LD
#' threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#' Default: \code{sqrt(0.1)}.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param verbose.v a \code{logical} indicating if information is shown during
#' the process in the \code{\link[SNPRelate]{snpgdsLDpruning}}() function.
#' Default: \code{FALSE}.
#'
#' @return a \code{list} of SNP IDs stratified by chromosomes as generated
#' by \code{\link[SNPRelate]{snpgdsLDpruning}} function.
#'
#' @details
#'
#' The SNP pruning is based on linkage disequilibrium (LD) and is done by the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function in the
#' \href{=https://bioconductor.org/packages/SNPRelate/}{SNPRelate} package.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#'
#' @importFrom SNPRelate snpgdsOpen snpgdsLDpruning
#' @importFrom gdsfmt closefn.gds
#' @encoding UTF-8
#' @keywords internal
runLDPruning <- function(gds, method=c("corr", "r", "dprime", "composite"),
                            listSamples=NULL,
                            listKeep=NULL,
                            slide.max.bp.v = 500000L,
                            ld.threshold.v=sqrt(0.1),
                            np=1L, verbose.v=FALSE) {

    ## Call SNP LD pruning
    snpset <- snpgdsLDpruning(gds, method="corr",
                                sample.id=listSamples,
                                snp.id=listKeep,
                                slide.max.bp=slide.max.bp.v,
                                ld.threshold=ld.threshold.v,
                                num.thread=np,
                                verbose=verbose.v)
    return(snpset)
}


#' @title Add the pruned.study entry in the GDS Sample file
#'
#' @description This function adds the pruned.study entry in the GDS
#' Sample file.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' GDS Sample file.
#'
#' @param pruned TODO
#'
#' @param sample.id a \code{character} string corresponding to
#' the sample identifier.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn index.gdsn delete.gdsn sync.gds ls.gdsn
#' @encoding UTF-8
#' @keywords internal
addGDSStudyPruning <- function(gds, pruned, sample.id) {

    ## Delete the pruned.study entry if present in the GDS Sample file
    if("pruned.study" %in% ls.gdsn(gds)) {
        delete.gdsn(index.gdsn(node=gds, "pruned.study"))
    }

    ## Create the pruned.study node in the GDS Sample file
    var.Pruned <- add.gdsn(gds, "pruned.study", pruned)

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
#' @importFrom gdsfmt add.gdsn index.gdsn ls.gdsn compression.gdsn append.gdsn sync.gds
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



#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened for the sample
#'
#' @param snp.lap TODO
#'
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
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


#' @title TODO
#'
#' @description TODO
#'
#' @param gds an object of class \code{gds} opened for the sample
#'
#' @param snp.seg TODO
#'
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' # TODO
#' gds <- "Demo GDS TODO"
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

    return(0L)
}




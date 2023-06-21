#' @title Extract the pruned SNVs in a reference data set (1KG) by chromosome
#' and/or allelic frequency
#'
#' @description The function extracts the pruned SNVs in a reference data
#' set (1KG) by chromosome and/or allelic frequency. The pruning is done
#' through the linkage disequilibrium analysis.
#'
#' @param gdsReference an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, an
#' opened SNP GDS file.
#'
#' @param method a \code{character} string that represents the method that
#' will be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime" and "composite". Default: \code{"corr"}.
#'
#' @param listSamples a \code{character} string that represents the method that
#' will be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime" and "composite". Default: \code{"corr"}.
#'
#' @param slideWindowMaxBP a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done in the \code{runLDPruning}
#' function. Default: \code{5e5}.
#'
#' @param thresholdLD a single \code{numeric} value that represents the LD
#' threshold used in the \code{runLDPruning} function.
#' Default: \code{sqrt(0.1)}.
#'
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#'
#' @param verbose a \code{logical} specifying if the function must provide
#' more information about the process. Default: \code{FALSE}.
#'
#' @param chr a \code{character} string representing the chromosome where the
#' selected SNVs should belong. Only one chromosome can be handled. If
#' \code{NULL}, the chromosome is not used as a filtering criterion.
#' Default: \code{NULL}.
#'
#' @param minAF a single positive \code{numeric} representing the
#' minimum allelic frequency used to select the SNVs. If \code{NULL}, the
#' allelic frequency is not used as a filtering criterion.
#' Default: \code{NULL}.
#'
#' @param outPrefix a \code{character} string that represents the prefix of the
#' RDS files that will be generated. Default: \code{"pruned_1KG"}.
#'
#' @param keepObj a \code{logical} specifying if the function must save the
#' the processed information into a RDS object. Default: \code{FALSE}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn
#' @encoding UTF-8
#' @keywords internal
pruning1KGbyChr <- function(gdsReference, method="corr", listSamples=NULL,
                    slideWindowMaxBP=5e5, thresholdLD=sqrt(0.1),
                    np=1L, verbose=FALSE, chr=NULL,
                    minAF=NULL, outPrefix="pruned_1KG", keepObj=FALSE) {

    filePruned <- file.path(paste0(outPrefix, ".rds"))
    fileObj <- file.path(paste0(outPrefix, "Obj.rds"))
    snpGDS <- index.gdsn(gdsReference, "snp.id")
    listKeep <- NULL
    if(is.null(minAF)){
        if(!is.null(chr)){
            snpGDS <- index.gdsn(gdsReference, "snp.id")
            snpID <- read.gdsn(snpGDS)

            chrGDS <- index.gdsn(gdsReference, "snp.chromosome")
            snpCHR <- read.gdsn(chrGDS)

            listKeep <- snpID[which(snpCHR == chr)]
        }
    } else{
        snpGDS <- index.gdsn(gdsReference, "snp.id")
        snpID <- read.gdsn(snpGDS)
        afGDS <- index.gdsn(gdsReference, "snp.AF")
        snpAF <- read.gdsn(afGDS)

        if(is.null(chr)){
            listKeep <- snpID[which(snpAF >= minAF & snpAF <= 1-minAF)]
        } else{
            chrGDS <- index.gdsn(gdsReference, "snp.chromosome")
            snpCHR <- read.gdsn(chrGDS)

            listKeep <- snpID[which(snpCHR == chr & snpAF >= minAF &
                                        snpAF <= 1-minAF)]
        }
    }

    snpset <- runLDPruning(gds=gdsReference, method=method,
                    listSamples=listSamples, listKeep=listKeep,
                    slideWindowMaxBP=slideWindowMaxBP,
                    thresholdLD=thresholdLD, np=np, verbose=verbose)

    pruned <- unlist(snpset, use.names=FALSE)
    saveRDS(pruned, filePruned)
    if(keepObj){
        saveRDS(snpset, fileObj)
    }

    return(0L)
}


#' @title Generate two indexes based on gene annotation for gdsAnnot1KG block
#'
#' @description TODO
#'
#' @param gdsReference an object of class
#' \link[gdsfmt]{gds.class} (a GDS file), the opened 1KG GDS file (reference).
#'
#' @param winSize a single positive \code{integer} representing the
#' size of the window to use to group the SNVs when the SNVs are in a
#' non-coding region. Default: \code{10000}.
#'
#' @param EnsDb An object of class \code{EnsDb} with the Ensembl genome
#' annotation. By default, the \code{EnsDb.Hsapiens.v86} class has been used.
#'
#' @return  a \code{data.frame} with those columns:
#' \itemize{
#' \item{chr} {}
#' \item{pos} {}
#' \item{snp.allele} {}
#' \item{Exon} {}
#' \item{GName} {}
#' \item{Gene} {}
#' \item{GeneS} {}
#' }
#' "chr", "pos", "snp.allele", "Exon", "GName", "Gene", "GeneS"
#' Example for GName and the two indexes "Gene", "GeneS"
#' GName Gene GeneS
#' 470                                 ENSG00000230021   17  3820
#' 471                                 ENSG00000230021   17  3820
#' 472                 ENSG00000230021:ENSG00000228794   17  3825
#' 473                 ENSG00000230021:ENSG00000228794   17  3825
#' 481 ENSG00000230021:ENSG00000228794:ENSG00000225880   17  3826
#' 482 ENSG00000230021:ENSG00000228794:ENSG00000225880   17  3826
#' 483 ENSG00000230021:ENSG00000228794:ENSG00000225880   17  3826
#' 492                 ENSG00000230021:ENSG00000228794   17  3825
#' 493                 ENSG00000230021:ENSG00000228794   17  3825
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alex Krasnitz
#' @importFrom S4Vectors Rle
#' @importFrom BSgenome strand
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom IRanges IRanges
#' @importFrom AnnotationDbi select
#' @importFrom ensembldb exonsBy toSAF genes
#' @importFrom AnnotationFilter GeneIdFilter
#' @encoding UTF-8
#' @keywords internal
generateGeneBlock <- function(gdsReference, winSize=10000, EnsDb) {

    edb <- EnsDb
    listEnsId <- unique(names(genes(edb)))

    cols <- c("GENEID", "SYMBOL", "GENENAME", "GENESEQSTART",
                    "GENESEQEND", "SEQNAME")

    annot <- select(edb, keys=listEnsId, columns=cols, keytype="GENEID")
    annot <- annot[which(annot$SEQNAME %in% c(seq_len(22), "X")),]

    # All the genes
    grGene <- GRanges(seqnames=annot$SEQNAME,
        ranges=IRanges(annot$GENESEQSTART, end=annot$GENESEQEND),
        strand=Rle(strand(rep("+", nrow(annot)))),
        mcols=annot[,c("GENEID", "GENENAME")])

    # Data frame of the all genes
    dfGenneAll <- as.data.frame(grGene)

    # group the overlapping gene
    grGeneReduce <- reduce(grGene)
    # data.frame version of grGeneReduce
    dfGene <- as.data.frame(grGeneReduce)

    # All exons
    allExon <- exonsBy(edb, by="gene", filter=GeneIdFilter(listEnsId))
    # Transforming the GRangesList into a data.frame in SAF format
    dfExon <- toSAF(allExon)
    # remove the duplicates
    dfExon <- unique(dfExon)
    # Group the overlaps
    exonReduce <- reduce(allExon)
    # Transforming the GRangesList into a data.frame in SAF format
    dfExonReduce <- toSAF(exonReduce)
    listMat <- list()

    matFreqAll <- data.frame(chr=read.gdsn(index.gdsn(gdsReference,
                                                        "snp.chromosome")),
                pos=read.gdsn(index.gdsn(gdsReference, "snp.position")),
                snp.allele=read.gdsn(index.gdsn(gdsReference, "snp.allele")),
                stringsAsFactors=FALSE)
    offsetGene <- 0
    offsetGeneS <- 0
    offsetGene.O <- 0

    for(chr in seq_len(22))
    {
        dfExonChr <- dfExonReduce[which(dfExonReduce$Chr == chr),]
        dfGenneAllChr <- dfGenneAll[which(dfGenneAll$seqnames == chr),]
        dfGeneChr <- dfGene[which(dfGene$seqnames == chr),]

        # colnames(matFreq) <- c("chr", "pos", "ref", "alt", "af", "EAS_AF",
        #                        "EUR_AF","AFR_AF", "AMR_AF", "SAS_AF")
        message(system.time({
            # SNV in the GDS
            matFreq <- matFreqAll[which(matFreqAll$chr == chr),]
            # create two vector (one for the exon and one for the gene) of char
            # with 1 entry for each SNV in the GDS
            # I will keep the name of the gene and exon at this position
            listSNVExons <- character(nrow(matFreq))
            listSNVGenes <- character(nrow(matFreq))

            listPos <- seq_len(nrow(matFreq))
            listPos <- listPos[order(matFreq$pos)]
            # Create an index to accelerate the process
            startIndex <- seq(1, nrow(matFreq), 1000)
            # Add if the last entry is not the last position
            # is not the nb row of matFreq add the the last
            #position
            if(startIndex[length(startIndex)] < nrow(matFreq)){
                startIndex <- c(startIndex, nrow(matFreq))
            }
            # For gene in the chr
            # slow but acceptable
            #    user  system elapsed
            #    26.116   0.074  26.201
            # see blockAnnotation.R for slower alternatives
            for (genePos in seq_len(nrow(dfGenneAllChr))) {
                # the gene is where SNV exists
                if (dfGenneAllChr$end[genePos] >= matFreq$pos[listPos[1]] &
                dfGenneAllChr$start[genePos] <= matFreq$pos[nrow(matFreq)]) {
                    # In which partitions from the index the gene is located
                    vStart <- max(c(which(matFreq$pos[startIndex] <=
                                            dfGenneAllChr$start[genePos]), 1))
                    vEnd <- min(c(which(matFreq$pos[startIndex] >=
                                            dfGenneAllChr$end[genePos]),
                                            length(startIndex)))
                    # List of SNV in the gene
                    listP <- which(matFreq$pos[listPos[startIndex[vStart]:startIndex[vEnd]]] >= dfGenneAllChr$start[genePos] &
                        matFreq$pos[listPos[startIndex[vStart]:startIndex[vEnd]]] <= dfGenneAllChr$end[genePos])

                    # if SNV in the gene
                    if (length(listP) > 0) {
                        # listPos in the gene
                        listP <-
                            listPos[startIndex[vStart]:startIndex[vEnd]][listP]

                        # Add the name of the gene of SNVs
                        listSNVGenes[listP] <- paste0(listSNVGenes[listP], ":",
                                        dfGenneAllChr$mcols.GENEID[genePos])

                        # Allow run on all without check if the SNV have
                        # already gene name
                        listSNVGenes[listP] <- gsub("^:", "",
                                                        listSNVGenes[listP])

                        # Exon of the gene
                        dfExon <- dfExonChr[which(dfExonChr$GeneID ==
                                    dfGenneAllChr$mcols.GENEID[genePos]),]
                        k <- 1

                        listE <- list()
                        for (pos in listP) {
                            if(length(which(dfExon$Start <= matFreq$pos[pos] &
                                    dfExon$End >= matFreq$pos[pos])) > 0) {
                                listE[[k]] <- pos
                                k <- k + 1
                            }
                        }

                        if (length(listE) > 0) {
                            listE <- do.call(c, listE)
                            listSNVExons[listE] <- paste0(listSNVExons[listE],
                                    ":", dfGenneAllChr$mcols.GENEID[genePos])
                            listSNVExons[listE] <- gsub("^:", "",
                                                    listSNVExons[listE])
                        }
                    }
                }
            }
        }))


        # add the column Exon with the list of gene with an exon with the SNV
        matFreq$Exon <- listSNVExons
        # add the column GName with the list of gene with which include the SNV
        matFreq$GName <- listSNVGenes

        # dfGeneChr are reduced (merge all the overlap interval)
        z <- cbind(c(dfGeneChr$start, dfGeneChr$end, as.integer(matFreq$pos)),
                c(seq_len(nrow(dfGeneChr)), -1 * seq_len(nrow(dfGeneChr)),
                    rep(0, nrow(matFreq))))
        z <- z[order(z[,1], -1 * z[,2]),]

        # group by interval which in overlap a gene
        matFreq$Gene[listPos] <- cumsum(z[,2])[z[,2] == 0]
        matFreq$Gene[matFreq$Gene > 0] <- matFreq$Gene[matFreq$Gene > 0] +
            offsetGene
        offsetGene <- max(offsetGene, max(matFreq$Gene))

        listD <- which(matFreq$Gene > 0)

        tmp <- paste0(matFreq[listD, "GName"], "_", matFreq[listD, "Gene"])
        listO <- order(tmp)


        # Create an index for each gene different if
        # two gene overlap the order don't have meaning.
        # gene ex: ENSG00000238009:ENSG00000239945_6 and ENSG00000238009_6
        # have a different number.
        # Note the order is base on the name not on position
        # Ex:
        #         GeneN                         indexNew
        # 135                 ENSG00000230021  4089
        # 136                 ENSG00000230021  4089
        # 148 ENSG00000230021:ENSG00000237973  4094
        # 149 ENSG00000230021:ENSG00000237973  4094
        # 159 ENSG00000229344:ENSG00000230021  4036
        # 160 ENSG00000229344:ENSG00000230021  4036
        # 161 ENSG00000230021:ENSG00000248527  4095
        # 162 ENSG00000198744:ENSG00000230021  3168
        # 163                 ENSG00000230021  4089
        # 164                 ENSG00000230021  4089
        # 165                 ENSG00000230021  4089
        #
        indexNew <- cumsum(!(duplicated(tmp[listO])))

        matFreq$GeneS <- rep(0, nrow(matFreq))
        matFreq$GeneS[listD][listO] <- indexNew + offsetGeneS
        offsetGeneS <- max(offsetGeneS, max(matFreq$GeneS))

        matFreq$GeneS[matFreq$GeneS < 0] <- 0
        matFreq$GeneS[matFreq$Gene < 0] <- 0
        listOrph <- which(matFreq$GeneS == 0)
        flag <- TRUE
        v <- offsetGene.O - 1
        i <- 1
        curZone <- "GeneS"
        curZone1 <- "Gene"
        winSize <- 10000

        if(length(listOrph) > 0){
            # Very slow can do better
            # but just run 1 time so less priority
            #    user  system elapsed
            # 517.595   7.035 524.658
            #    user  system elapsed
            # 558.526   2.274 561.043
            while(flag){
                #use the index
                vStart <- min(c(which(matFreq$pos[startIndex] >
                                    (matFreq[listOrph[i], "pos"] + winSize)),
                                        length(startIndex)))

                preList <- listOrph[i]:startIndex[vStart]
                listWin <- which(matFreq[preList, "pos"] >
                                    (matFreq[listOrph[i], "pos"] + winSize) |
                                    (matFreq[preList, "pos"] >
                                        matFreq[listOrph[i], "pos"] &
                                        matFreq[preList,"GeneS"] > 0))

                j <- ifelse(length(listWin) > 0, preList[listWin[1]] - 1,
                                listOrph[i])

                matFreq[listOrph[i]:j, curZone] <- v
                matFreq[listOrph[i]:j, curZone1] <- v
                v <- v - 1
                i <- which(listOrph == j) + 1
                flag <- ifelse(i <= length(listOrph), TRUE, FALSE)
            }
            offsetGene.O <- min(offsetGene.O, min(matFreq$Gene))
        }

        listMat[[chr]] <- matFreq

        # save the matrix for each chr
        # create the space at the begining
    }

    matGene.Block <- do.call(rbind, listMat)
    rm(listMat)
    return(matGene.Block)
}



#' @title Validate input parameters for prepPed1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{prepPed1KG}} function.
#'
#' @param filePed a \code{character} string representing the path and
#' file name of the pedigree file (PED file) that contains the information
#' related to the profiles present in the 1KG GDS file. The PED file must
#' exist.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the 1KG genotyping files for each profile are located. Only the profiles
#' with associated genotyping files are retained in the creation of the final
#' \code{data.frame}. The name of the genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file
#' (PED file).
#'
#' @param batch a\code{integer} that uniquely identifies the source of the
#' pedigree information. The 1KG is usually \code{0L}.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo pedigree file
#' pedDemoFile <- file.path(dataDir, "PedigreeDemo.ped")
#'
#' ## The validation should be successful
#' RAIDS:::validatePrepPed1KG(filePed=pedDemoFile,
#'      pathGeno=dataDir, batch=1)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validatePrepPed1KG <- function(filePed, pathGeno, batch) {

    ## Validate that the batch is an integer
    if (! isSingleNumber(batch)) {
        stop("The batch must be an integer.")
    }

    ## Validate that the pedigree file exists
    if (! file.exists(filePed)) {
        stop("The file \'", filePed, "\' does not exist." )
    }

    ## Validate that the path for the genotyping files exists
    if (!(is.character(pathGeno) && dir.exists(pathGeno))) {
        stop("The \'pathGeno\' parameter must be a character string ",
                "representing an existing directory.")
    }

    return(0L)
}


#' @title Validate input parameters for generateGDS1KG() function
#'
#' @description This function validates the input parameters for the
#' \code{\link{generateGDS1KG}} function.
#'
#' @param pathGeno a \code{character} string representing the path where
#' the 1K genotyping files for each sample are located. The name of the
#' genotyping files must correspond to
#' the individual identification (Individual.ID) in the pedigree file.
#'
#' @param filePedRDS a \code{character} string representing the path and file
#' name of the RDS file that contains the pedigree information. The file must
#' exist. The file must be a RDS file.
#'
#' @param fileSNVIndex a \code{character} string representing the path and file
#' name of the RDS file that contains the indexes of the retained SNPs. The
#' file must exist. The file must be a RDS file.
#'
#' @param fileSNVSelected a \code{character} string representing the path and
#' file name of the RDS file that contains the filtered SNP information. The
#' file must exist. The file must be a RDS file.
#'
#' @param fileNameGDS a \code{character} string representing the path and file
#' name of the GDS file that will be created. The GDS file will contain the
#' SNP information, the genotyping information and the pedigree information
#' from 1000 Genomes.
#' The extension of the file must be '.gds'.
#'
#' @param listSamples a \code{vector} of \code{character} string corresponding
#' to samples (must be the sample.ids) that will be retained and added to the
#' GDS file. When \code{NULL}, all the samples are retained.
#'
#' @param verbose a \code{logical} indicating if the function must print
#' messages when running.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The RDS file containing the pedigree information
#' pedigreeFile <- file.path(dataDir, "PedigreeDemo.rds")
#'
#' ## The RDS file containing the indexes of the retained SNPs
#' snpIndexFile <- file.path(dataDir, "listSNPIndexes_Demo.rds")
#'
#' ## The RDS file containing the filtered SNP information
#' filterSNVFile <- file.path(dataDir, "mapSNVSelected_Demo.rds")
#'
#' ## Temporary GDS file containing 1KG information
#' gdsFile <- file.path(dataDir, "1KG_TEMP.gds")
#'
#' ## The validation should be successful
#' RAIDS:::validateGenerateGDS1KG(pathGeno=dataDir, filePedRDS=pedigreeFile,
#'     fileSNVIndex=snpIndexFile, fileSNVSelected=filterSNVFile,
#'     fileNameGDS=gdsFile, listSamples=NULL, verbose=FALSE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateGenerateGDS1KG <-  function(pathGeno, filePedRDS, fileSNVIndex,
                                        fileSNVSelected, fileNameGDS,
                                        listSamples, verbose) {

    validateCharacterString(pathGeno, "pathGeno")

    ## Validate that the path for the genotyping files exists
    if (! file.exists(pathGeno)) {
        stop("The path \'", pathGeno, "\' does not exist." )
    }

    validateCharacterString(filePedRDS, "filePedRDS")

    ## Validate that the pedigree file exists
    if (! file.exists(filePedRDS)) {
        stop("The file \'", filePedRDS, "\' does not exist." )
    }

    validateCharacterString(fileSNVIndex, "fileSNVIndex")

    ## Validate that the SNP indexes file exists
    if (! file.exists(fileSNVIndex)) {
        stop("The file \'", fileSNVIndex, "\' does not exist." )
    }

    validateCharacterString(fileSNVSelected, "fileSNVSelected")

    ## Validate that the SNP information file exists
    if (! file.exists(fileSNVSelected)) {
        stop("The file \'", fileSNVSelected, "\' does not exist." )
    }

    validateCharacterString(fileNameGDS, "fileNameGDS")

    ## Validate that the listSamples is null or character string
    if (!(is.null(listSamples) || is.character(listSamples))) {
        stop("The \'listSamples\' must be NULL or a vector of character ",
                    "strings." )
    }

    ## The verbose parameter must be a logical
    validateLogical(logical=verbose, "verbose")

    return(0L)
}



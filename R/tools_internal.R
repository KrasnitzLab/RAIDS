#' @title Validate that the input parameter is a GDS object
#'
#' @description This function validates that the input parameter inherits
#' the \link[gdsfmt]{gds.class} class.
#'
#' @param gds an object of class \link[gdsfmt]{gds.class} (a GDS file), the
#' 1 KG GDS file.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## The 1KG GDS file (opened)
#' gds1KG <- openfn.gds(file.path(dataDir,
#'                 "PopulationReferenceDemo.gds"), readonly=TRUE)
#'
#' ## The validation should be successful
#' RAIDS:::validateGDSClass(gds=gds1KG, name="gds")
#'
#' ## All GDS file must be closed
#' closefn.gds(gdsfile=gds1KG)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateGDSClass <- function(gds, name) {

    ## The gds must be an object of class "gds.class"
    if (!inherits(gds, "gds.class")) {
        stop("The \'", name ,"\' must be an object of class \'gds.class\'.")
    }

    return(0L)
}


#' @title Validate that the input parameter is a logical
#'
#' @description This function validates that the input parameter is
#' a logical. If the parameter is not a logicial, the function generates
#' an error with a specific message.
#'
#' @param logical a \code{logical} that will be validated.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## The validation should be successful
#' RAIDS:::validateLogical(logical=TRUE, name="test")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateLogical <- function(logical, name) {

    ## The parameter must be a logical
    if (!is.logical(logical)) {
        stop("The \'", name, "\' parameter must be a logical (TRUE or FALSE).")
    }

    return(0L)
}


#' @title Validate that the input parameter is a single positive numeric
#' between zero and one (included)
#'
#' @description This function validates that the input parameter is
#' a single numeric between zero and one (included). If the parameter is not,
#' the function generates an error with a specific message.
#'
#' @param value a single positive \code{numeric} that will be validated.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## The validation should be successful
#' RAIDS:::validateSingleRatio(value=0.02, name="test")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom S4Vectors isSingleNumber
#' @encoding UTF-8
#' @keywords internal
validateSingleRatio <- function(value, name) {

    ## The parameter must be a single positive numeric between 0 and 1
    if (!(isSingleNumber(value) && (value >= 0.0) && (value <= 1.0))) {
        stop("The \'", name, "\' must be a single numeric positive ",
                "value between 0 and 1.")
    }

    return(0L)
}


#' @title Validate that the input parameter is a character sting
#'
#' @description This function validates that the input parameter is
#' a character string (vector of 1 entry). If the parameter is not a
#' character string, the function generates
#' an error with a specific message.
#'
#' @param value a \code{character} string that will be validated.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## The validation should be successful
#' RAIDS:::validateCharacterString(value="hi", name="test")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validateCharacterString <- function(value, name) {

    ## The parameter value must be a character string
    if(!(is.character(value) && length(value) == 1)) {
        stop("The \'", name, "\' parameter must be a character string.")
    }

    return(0L)
}


#' @title Validate that the input parameter is a vector of positive
#' \code{numeric}
#'
#' @description This function validates that the input parameter is
#' a \code{vector} of positive \code{numeric} values (vector of 1 entry or
#' more). All values have to be positive (>\code{0}).
#' If the parameter is not respecting the validation, the function generates
#' an error with a specific message.
#'
#' @param value a \code{vector} of \code{numeric} that will be validated.
#'
#' @param name a \code{character} string that represents the name of the
#' parameter that is tested.
#'
#' @return The function returns \code{0L} when successful.
#'
#' @examples
#'
#' ## The validation should be successful
#' RAIDS:::validatePositiveIntegerVector(value=c(1, 2 ,3, 5), name="test")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
validatePositiveIntegerVector <- function(value, name) {

    ## The parameter value must be a character string
    if(!(is.numeric(value) && is.vector(value) && all(value > 0))) {
        stop("The \'", name, "\' parameter must be a vector of positive ",
                    "integers.")
    }

    return(0L)
}

#' @title Read a SNP-pileup file
#'
#' @description The function reads a generic SNP pileup file and
#' returns a data frame
#' containing the information about the read counts for the SNVs present in
#' the file.
#'
#' @param fileName a \code{character} string representing the name, including
#' the path, of a text file containing the SNV read counts as generated by
#' snp-pileup software. The text file must be comma separated.
#' The text file must contain those columns: Chromosome, Position, Ref, Alt,
#' File1R, File1A, File1E and File1D.
#'
#' @param offset a \code{integer} representing the offset to be added to the
#' position of the SNVs. The value of offset
#' is added to the position present in the file. Default: \code{0L}.
#'
#' @return the a \code{data.frame} containing at least:
#' \describe{
#' \item{Chromosome}{ a \code{numeric} representing the name of
#' the chromosome}
#' \item{Position}{ a \code{numeric} representing the position on the
#' chromosome}
#' \item{Ref}{ a \code{character} string representing the reference nucleotide}
#' \item{Alt}{ a \code{character} string representing the alternative
#' nucleotide}
#' \item{File1R}{ a \code{numeric} representing the count for
#' the reference nucleotide}
#' \item{File1A}{ a \code{numeric} representing the count for the
#' alternative nucleotide}
#' \item{File1E}{a \code{numeric} representing the count for the
#' errors}
#' \item{File1D}{a \code{numeric} representing the count for the
#' deletions}
#' \item{count}{ a \code{numeric} representing the total count}
#' }
#'
#' @examples
#'
#' ## Directory where demo SNP-pileup file
#' dataDir <- system.file("extdata/example/snpPileup", package="RAIDS")
#'
#' ## The SNP-pileup file
#' snpPileupFile <- file.path(dataDir, "ex1.txt.gz")
#'
#' info <- RAIDS:::readSNVPileupFile(fileName=snpPileupFile)
#' head(info)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
readSNVPileupFile <- function(fileName, offset = 0L) {

    matSample <- read.csv(fileName)

    # Check if the mandatory column are there
    if (!(all(c("Chromosome", "Position", "Ref", "Alt", "File1R",
                "File1A", "File1E", "File1D") %in% colnames(matSample)))) {
        stop("The SNP-pileup file must contain all those columns: ",
                "\'Chromosome\', \'Position\', \'Ref\', \'Alt\', \'File1R\', ",
                "\'File1A\', \'File1E\', \'File1D\'.")
    }

    matSample[, "Chromosome"] <- as.integer(gsub("chr", "",
                                            matSample[, "Chromosome"]))
    matSample[, "Position"] <- matSample[, "Position"] + offset
    matSample[, "count"] <- rowSums(matSample[, c("File1R", "File1A",
                                            "File1E", "File1D")])

    return(matSample)
}


#' @title Read a generic SNP pileup file
#'
#' @description The function reads a generic SNP pileup file and
#' returns a data frame
#' containing the information about the read counts for the SNVs present in
#' the file.
#'
#' @param fileName a \code{character} string representing the name, including
#' the path, of a text file containing the SNV read counts. The text file must
#' be comma separated. The text file must
#' contain those columns: Chromosome, Position, Ref, Alt, Count,
#' File1R and File1A.
#'
#' @param offset a \code{integer} representing the offset to be added to the
#' position of the SNVs. The value of offset
#' is added to the position present in the file. Default: \code{0L}.
#'
#' @return a \code{data.frame} containing at least:
#' \describe{
#' \item{Chromosome}{ a \code{numeric} representing the name of
#' the chromosome}
#' \item{Position}{ a \code{numeric} representing the position on the
#' chromosome}
#' \item{Ref}{ a \code{character} string representing the reference nucleotide}
#' \item{Alt}{ a \code{character} string representing the alternative
#' nucleotide}
#' \item{File1R}{ a \code{numeric} representing the count for
#' the reference nucleotide}
#' \item{File1A}{ a \code{numeric} representing the count for the
#' alternative nucleotide}
#' \item{count}{ a \code{numeric} representing the total count}
#' }
#'
#' @examples
#'
#'
#' ## Directory where demo SNP-pileup file
#' dataDir <- system.file("extdata/example/snpPileup", package="RAIDS")
#'
#' ## The SNP-pileup file
#' snpPileupFile <- file.path(dataDir, "ex1.generic.txt.gz")
#'
#' info <- RAIDS:::readSNVFileGeneric(fileName=snpPileupFile)
#' head(info)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt add.gdsn
#' @encoding UTF-8
#' @keywords internal
readSNVFileGeneric <- function(fileName, offset = 0L) {

    matSample <- read.csv(fileName)

    # Check if the mandatory column are there
    if (!(all(c("Chromosome", "Position", "Ref", "Alt", "File1R",
                    "File1A", "Count") %in% colnames(matSample)))) {
        stop("The generic SNP pileup file must contain all those columns: ",
                "\'Chromosome\', \'Position\', \'Ref\', \'Alt\', \'File1R\', ",
                "\'File1A\', \'Count\'.")
    }

    matSample[, "Chromosome"] <- as.integer(gsub("chr", "",
                                        matSample[, "Chromosome"]))
    matSample[, "Position"] <- matSample[, "Position"] + offset
    colnames(matSample)[colnames(matSample) == "Count"] <- "count"

    return(matSample)
}

#' @title Read a VCF file with the genotypes use for the ancestry call
#'
#' @description The function reads VCF file and returns a data frame
#' containing the information about the read counts for the SNVs present in
#' the file.
#'
#' @param fileName a \code{character} string representing the name, including
#' the path, of a VCF file containing the SNV read counts.
#' The VCF must contain those genotype fields: GT, AD, DP.
#'
#' @param profileName a \code{character} with Name.ID for the genotype name.
#' Default: \code{NULL}.
#'
#' @param offset a \code{integer} representing the offset to be added to the
#' position of the SNVs. The value of offset
#' is added to the position present in the file. Default: \code{0L}.
#'
#' @return a \code{data.frame} containing at least:
#' \describe{
#' \item{Chromosome}{ a \code{numeric} representing the name of
#' the chromosome}
#' \item{Position}{ a \code{numeric} representing the position on the
#' chromosome}
#' \item{Ref}{ a \code{character} string representing the reference nucleotide}
#' \item{Alt}{ a \code{character} string representing the alternative
#' nucleotide}
#' \item{File1R}{ a \code{numeric} representing the count for
#' the reference nucleotide}
#' \item{File1A}{ a \code{numeric} representing the count for the
#' alternative nucleotide}
#' \item{count}{ a \code{numeric} representing the total count}
#' }
#'
#' @examples
#'
#'
#' ## Directory where demo SNP-pileup file
#' dataDir <- system.file("extdata/example/snpPileup", package="RAIDS")
#'
#' ## The SNP-pileup file
#' snpPileupFile <- file.path(dataDir, "ex1.vcf.gz")
#'
#' info <- RAIDS:::readSNVVCF(fileName=snpPileupFile)
#' head(info)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom VariantAnnotation readVcf geno
#' @importFrom MatrixGenerics rowRanges
#' @importFrom GenomicRanges seqnames start width
#' @encoding UTF-8
#' @keywords internal
readSNVVCF <- function(fileName, profileName=NULL, offset=0L) {

    vcf <- readVcf(fileName)

    gtCur <- geno(vcf)
    genoPos <- 1
    if(! is.null(profileName)){
        listVcfSample <- colnames(gtCur$GT)
        genoPos <- which(listVcfSample == profileName)
    }

    colChr <- as.character(seqnames(vcf))

    infoPos <- rowRanges(vcf)
    start <- start(infoPos)
    widthRef <- width(infoPos$REF)
    refCur <- as.character(infoPos$REF)
    listKeep <- seq_len(length(refCur))
    listKeep <- listKeep[which(widthRef == 1)]

    countV <- as.integer(gtCur$DP)
    countA <- gtCur$AD

    idVCF <- row.names(gtCur$GT)
    tmp <- matrix(unlist(strsplit(idVCF, ":")), nrow=2)[2,]
    alleleChar <- matrix(unlist(strsplit(tmp, "_")), nrow=2)[2,]
    rm(tmp)

    matCur <- lapply(listKeep, FUN=function(x, countA, alleleChar){

            listAlt <- strsplit(alleleChar[x], "\\/")[[1]]

            keep <- ifelse(listAlt[2] %in% c("A", "C", "G", "T"), TRUE, FALSE)

            res <- data.frame(Alt=as.character(listAlt[2]),
                        File1R=countA[x,genoPos][[1]][1], File1A=countA[x,genoPos][[1]][2],
                        keep=keep)

            return(res)
        },
        countA=countA,
        alleleChar=alleleChar)

    matCur <- do.call(rbind, matCur)
    listTmp <- which(matCur$keep)
    listKeep <- listKeep[listTmp]
    matSample <- data.frame(Chromosome=as.integer(gsub("chr", "",
                                                        colChr))[listKeep],
                    Position=start[listKeep] + offset,
                    Ref=refCur[listKeep], Alt=matCur$Alt[listTmp],
                    File1R=matCur$File1R[listTmp],
                    File1A=matCur$File1A[listTmp],
                    count=unname(gtCur$DP[listKeep,genoPos]),
                    stringsAsFactors=FALSE)

    return(matSample)
}




#' @title Filtering the read counts for a specific nucleotide
#'
#' @description The function returns the read counts for the specific
#' nucleotide or zero when read counts are not available.
#'
#' @param nucleotide a \code{vector} of  a \code{character} strings
#' representing the nucleotides (ex: A, C, G or T).
#'
#'
#' @param count a \code{vector} of \code{numeric} representing the counts for
#' each nucleotide listed in \code{nucleotide} parameter.
#'
#' @param curNucleo a \code{character} strings representing the nucleotide
#' that will be retained (ex: A, C, G or T).
#'
#' @return a \code{numeric} representing the counts for the selected
#' nucleotide. The default value is \code{0}.
#'
#' @examples
#'
#' ## Nucleotides vector
#' nuc <- c("A", "G", "C", "T")
#'
#' ## Count vector
#' cnt <- c(100, 200, 4, 32)
#'
#' ## Return the count for the nucleotide "G"
#' RAIDS:::extractNucleotide(nucleotide=nuc, count=cnt, curNucleo="G")
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
extractNucleotide <- function(nucleotide, count, curNucleo) {

    tmp <- which(nucleotide == curNucleo)
    res <- 0
    if(length(tmp) == 1) res <- count[tmp]
    return(res)

}


#' @title Extract SNV information from pileup file for a selected chromosome
#'
#' @description The function reads pileup file and
#' returns a \code{data.frame}
#' containing the information about the read counts for the SNVs present in
#' the selected chromosome.
#'
#' @param chr a \code{character} string representing the name of the
#' chromosome to keep
#'
#' @param resPileup a \code{data.frame} as generated by the \code{pileup}
#' function from \code{Rsamtools} package
#'
#' @param varDf a \code{list} containing a \code{data.frame} representing
#' the positions to keep for each chromosome.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#'
#' @return a \code{data.frame} containing at least:
#' \describe{
#' \item{seqnames}{ a \code{character} representing the name of the chromosome}
#' \item{pos}{ a \code{numeric} representing the position on the chromosome}
#' \item{REF}{ a \code{character} string representing the reference nucleotide}
#' \item{ALT}{ a \code{character} string representing the alternative
#' nucleotide}
#' \item{A}{ a \code{numeric} representing the count for the A nucleotide}
#' \item{C}{ a \code{numeric} representing the count for the C nucleotide}
#' \item{G}{ a \code{numeric} representing the count for the G nucleotide}
#' \item{T}{ a \code{numeric} representing the count for the T nucleotide}
#' \item{count}{ a \code{numeric} representing the total count}
#' }
#'
#' @examples
#'
#' ## Demo pileup result data.frame
#' resDemo <- data.frame(seqnames=rep("chr14", 10),
#'     pos=c(19069583, 19069584, 19069586, 19069588, 19069589, 19069590,
#'             19069591, 19069592, 19069609, 19069760),
#'     strand=c(rep("+", 5), rep("-", 5)),
#'     nucleotide=c("T", "G", "G", "C", "A", "A", "C", "T", "T", "G"),
#'     count=c(5, 3, 2, 4, 1, 2, 1, 8, 7, 4))
#' resDemo$seqnames <- factor(resDemo$seqnames)
#' resDemo$strand <- factor(resDemo$strand)
#' resDemo$nucleotide <- factor(resDemo$nucleotide)
#'
#' ## Position to keep in a data.frame format
#' varInfo <- list("chr14"=data.frame(chr=c("chr14", "chr14"),
#'     start=c(19069584, 19069609), REF=c("A", "G"), ALT=c("T", "A")))
#'
#' ## Extract information from pileup for selected positions
#' RAIDS:::processPileupChrBin(chr="chr14", resPileup=resDemo, varDf=varInfo,
#'     verbose=FALSE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom dplyr group_by summarize %>%
#' @importFrom rlang .data
#' @encoding UTF-8
#' @keywords internal
processPileupChrBin <- function(chr, resPileup, varDf, verbose) {

    resCur <- NULL

    ## Assign FALSE to verbose by default
    if (is.null(verbose)) {
        verbose <- FALSE
    }

    if (chr %in% names(varDf)) {

        keep <- which(resPileup$seqnames == chr)

        if (length(keep) > 0) {
            # restrict the resPileup to the chromosome chr
            snpO <- resPileup[keep,]
            rm(keep)

            # restrict. varDf to the region overlapping snpO
            vcfCur <- varDf[[chr]][varDf[[chr]]$start >= min(snpO$pos) &
                            varDf[[chr]]$start <= max(snpO$pos),, drop=FALSE]

            if (nrow(vcfCur) > 0) {

                # Get the positions to keep in resPileup (snpO)
                tmpTime <- system.time( {z <- cbind( c(vcfCur$start,
                                                snpO$pos, vcfCur$start),
                    c(rep(-1, nrow(vcfCur)), rep(0, nrow(snpO)),
                        rep(1, nrow(vcfCur))),
                    c(seq_len(nrow(vcfCur)), seq_len(nrow(snpO)),
                        seq_len(nrow(vcfCur))))
                    z <- z[order(z[,1]),]
                    listKeep <- which(cumsum(z[,2]) < 0 & z[,2]==0)} )

                if (verbose) {
                    message("processPileupChrBin selected pos user ",
                        chr, " ",
                        round(tmpTime[1],3), " system ", round(tmpTime[2],3),
                        " elapsed ", round(tmpTime[3],3)) }

                # summarize by position with the for base
                tmpTime <- system.time(resCur <- as.data.frame(snpO[z[listKeep, 3],] %>%
                    group_by(.data$pos) %>%
                        summarize(seqnames=.data$seqnames[1],
                        A=extractNucleotide(.data$nucleotide,.data$count, "A"),
                        C=extractNucleotide(.data$nucleotide,.data$count, "C"),
                        G=extractNucleotide(.data$nucleotide,.data$count, "G"),
                        T=extractNucleotide(.data$nucleotide,.data$count, "T"),
                        count=sum(.data$count))))

                if (verbose) {
                    message("processPileupChrBin extracted nucleotides user ",
                        round(tmpTime[1],3), " system ", round(tmpTime[2],3),
                        " elapsed ", round(tmpTime[3],3)) }
                if(length(listKeep) > 0){
                    # Add the reference allele and the alternative allele
                    tmpTime <- system.time({z <- cbind(c(resCur$pos, vcfCur$start, resCur$pos),
                        c(rep(-1, nrow(resCur)), rep(0, nrow(vcfCur)), rep(1, nrow(resCur))),
                        c(seq_len(nrow(resCur)), seq_len(nrow(vcfCur)), seq_len(nrow(resCur))))
                        z <- z[order(z[,1]),]
                        listKeep <- which(cumsum(z[,2]) < 0 & z[,2]==0)
                        resCur$REF <- rep("N", nrow(resCur))
                        resCur$ALT <- rep("N", nrow(resCur))
                        resCur$REF[z[listKeep-1,3]] <- vcfCur[z[listKeep,3], "REF"]
                        resCur$ALT[z[listKeep-1,3]] <- vcfCur[z[listKeep,3], "ALT"]
                        resCur <- resCur[,c("seqnames", "pos", "REF", "ALT", "A",
                                                "C", "G", "T", "count")]} )

                    if(verbose) {
                        message("processPileupChrBin add ref and alt allele user ",
                            round(tmpTime[1],3), " system ", round(tmpTime[2],3),
                            " elapsed ", round(tmpTime[3],3)) }
                } else{
                    resCur <- NULL
                }
            }
        }
    }
    if(verbose) {message("readSNVBAM processPileupChrBin before return ",
                         Sys.time())}
    return(resCur)
}

#' @title Read a VCF file with the genotypes use for the ancestry call
#'
#' @description The function reads VCF file and
#' returns a data frame
#' containing the information about the read counts for the SNVs present in
#' the file.
#'
#' @param fileName a \code{character} string representing the name, including
#' the path, of a BAM file with the index file in the same directory
#'
#'
#' @param paramSNVBAM a \code{list} containing the parameters passed to the
#' BamFile() function. Default: \code{list(ScanBamParam=NULL, PileupParam=NULL,
#' yieldSize=10000000)}.
#'
#' @param varSelected a \code{data.frame} representing the position to keep
#'
#' @param offset a \code{integer} representing the offset to be added to the
#' position of the SNVs. The value of offset
#' is added to the position present in the file. Default: \code{0L}.
#'
#' @param verbose a \code{logical} indicating if messages should be printed
#' to show how the different steps in the function. Default: \code{FALSE}.
#'
#' @return a \code{data.frame} containing at least:
#' \describe{
#' \item{Chromosome}{ a \code{numeric} representing the name of
#' the chromosome}
#' \item{Position}{ a \code{numeric} representing the position on the
#' chromosome}
#' \item{Ref}{ a \code{character} string representing the reference nucleotide}
#' \item{Alt}{ a \code{character} string representing the alternative
#' nucleotide}
#' \item{File1R}{ a \code{numeric} representing the count for
#' the reference nucleotide}
#' \item{File1A}{ a \code{numeric} representing the count for the
#' alternative nucleotide}
#' \item{count}{ a \code{numeric} representing the total count}
#' }
#'
#' @examples
#'
#'
#' ## Required library for this example to run correctly
#' if (requireNamespace("Rsamtools", quietly=TRUE)) {
#'     ## Demo bam
#'     fl <- system.file("extdata", "no_which_buffered_pileup.bam",
#'             package="Rsamtools", mustWork=TRUE)
#'
#'     RAIDS:::readSNVBAM(fl, varSelected=data.frame(chr=c(1,1),
#'                 start=c(3,5), REF=c("A", "A"), ALT=c("C", "C")))
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom Rsamtools BamFile ScanBamParam PileupParam pileup
#' @importFrom MatrixGenerics rowRanges
#' @importFrom GenomicRanges seqnames start width
#' @encoding UTF-8
#' @keywords internal
readSNVBAM <- function(fileName, varSelected, offset=0L,
                        paramSNVBAM=list(ScanBamParam=NULL, PileupParam=NULL,
                            yieldSize=10000000), verbose=FALSE) {
    # Note the offset is apply to the ref not the sequemce (snp-pileup and vcf)
    varSelected$chr <- paste0("chr", varSelected$chr)
    varSelected$start <- varSelected$start - offset
    listChr <- unique(varSelected$chr)
    # listChr <- names(varSelected)
    varSelected <- lapply(listChr,
                    FUN=function(x, varDf){
                        return(varDf[which(varDf$chr == x),])
                    },
                    varDf=varSelected)

    # names(listPos) <- paste0("chr", listChr)
    # varSelected <- lapply(seq_len(length(varSelected)),
    #                       FUN=function(x, varSelected){
    #                           varSelected[[x]]$chr <- paste0("chr", varSelected[[x]]$chr)
    #                           varSelected[[x]]$start <- varSelected[[x]]$start + offset
    #                           return(varSelected[[x]])
    #                       },
    #                       varSelected=varSelected)
    names(varSelected) <- listChr
    #varSelected$chr <- paste0("chr", varSelected$chr)
    #varSelected$start <- varSelected$start - offset

    myBf <- BamFile(fileName, yieldSize=paramSNVBAM$yieldSize)
    bf <- open(myBf)
    # temporary before create a parameters class
    if(is.null(paramSNVBAM$ScanBamParam)){
        sbp <- ScanBamParam() # which=vcf_granges
    }
    if(is.null(paramSNVBAM$ScanBamParam)){
        pup <- PileupParam(max_depth=5000,
                           min_base_quality=20,
                           min_mapq=15,
                           min_nucleotide_depth=1,
                           min_minor_allele_depth=0,
                           distinguish_strands=FALSE,
                           distinguish_nucleotides=TRUE,
                           ignore_query_Ns=TRUE,
                           include_deletions=FALSE,
                           include_insertions=FALSE)
    }

    i<-1
    res <- list()
    vcfChr <- list()
    k<-0

    repeat {
        if(verbose) {message("readSNVBAM pileup user ",
                             k, " i ", i)}
        tmpTime <- system.time(resPileup <- pileup(bf,
                                                   pileupParam = pup,
                                                   ScanBamParam=sbp))
        if(verbose) {message("readSNVBAM pileup user ",
                             round(tmpTime[1],3),
                             " system ", round(tmpTime[2],3),
                             " elapsed ", round(tmpTime[3],3))}

        listChr <- unique(as.character(resPileup$seqnames))


        #print(paste0("nb Chr ", length(listChr)))
        res[[i]] <- NULL
        if(length(listChr) > 0 ) {
            tmpChr <- grep("chr", listChr)
            if(length(tmpChr) != length(listChr)){
                listChg <- listChr[-1* tmpChr]
                for(i in seq_len(length(listChg))){
                    resPileup$seqnames[resPileup$seqnames == listChg[i]] <-
                            paste0("chr", listChg[i])
                }
                listChr <- unique(as.character(resPileup$seqnames))
            }
            if(sum(!(listChr %in% names(varSelected))) == length(listChr)) {
                # message("End chromosome")
                break
            }
            # print(paste0("Current chr ", listChr))
            if(verbose) { message("readSNVBAM processPileupChrBin start ",
                                 Sys.time()) }
            tmp <- lapply(listChr,
                    FUN=function(x, res, varSelected){
                            return(processPileupChrBin(chr=x, res,
                                varDf=varSelected, verbose=verbose))
                          }, res=resPileup,
                          varSelected=varSelected)
            if(verbose) { message("readSNVBAM processPileupChrBin end ",
                                 Sys.time()) }

            if(length(tmp) > 0) {
                res[[i]] <- do.call(rbind, tmp)
                i <- i + 1
                #message(nrow(res[[i]]), " rows in result data.frame")
            }
        }

        if(nrow(resPileup) == 0L){
            k <- k + 1
            if(k > 20) {
                break
            }
        }else{
            k <- 0
            # i <- i + 1
        }
        if(verbose) { message("readSNVBAM pileup end repeat ", Sys.time()) }
    }
    if(verbose) {message("readSNVBAM pileup Done ", Sys.time())}
    resSNP <- do.call(rbind,res)
    close(bf)

    if(verbose) {message("readSNVBAM pileup user ", k, " i ", i)}

    resSNP$File1R <- rep(0, nrow(resSNP))
    resSNP$File1A <- rep(0, nrow(resSNP))
    for(nuc in c("A", "C", "G", "T")){
        tmp <- which(resSNP$REF == nuc)
        if(length(tmp) > 0){
            resSNP$File1R[tmp] <- resSNP[tmp, nuc]
        }
        tmp <- which(resSNP$ALT == nuc)
        if(length(tmp) > 0){
            resSNP$File1A[tmp] <- resSNP[tmp, nuc]
        }
    }
    resSNP <- resSNP[, c("seqnames", "pos", "REF", "ALT", "File1R", "File1A",
                                "count", "A", "C", "G", "T")]
    colnames(resSNP) <- c("Chromosome", "Position", "Ref", "Alt", "File1R",
                                "File1A", "count", "A", "C", "G", "T")
    resSNP$Position <- resSNP$Position  + offset

    if(verbose) {message("readSNVBAM pileup format Done ", Sys.time())}
    return(resSNP)
}



#' @title The function create a vector of integer representing the linkage
#' disequilibrium block for each SNV in the in the same order
#' than the variant in Population reference dataset.
#'
#' @description The function create a vector of integer representing the linkage
#' disequilibrium block for each SNV in the in the same order
#' than the variant in Population reference dataset.
#'
#' @param fileReferenceGDS  a \code{character} string representing the file
#' name of the Reference GDS file. The file must exist.
#'
#' @param fileBlock a \code{character} string representing the file
#' name of output file det from the plink block command for a chromosome.
#'
#' @return a \code{list} containing 2 entries:
#' \describe{
#' \item{\code{chr}}{ a \code{integer} representing a the chromosome from
#' fileBlock.
#' }
#' \item{\code{block.snp}}{ a \code{array} of integer
#' representing the linkage disequilibrium block for
#' each SNV in the in the same order than the variant
#' in Population reference dataset.
#' }
#' }
#'
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo of Reference GDS file containing reference information
#' fileReferenceGDS <- file.path(dataDir, "PopulationReferenceDemo.gds")
#'
#' ## Demo of of output file det from the plink block
#' ## command for chromosome 1
#' fileLdBlock <- file.path(dataDir, "block.sp.EUR.Ex.chr1.blocks.det")
#'
#' listLdBlock <- RAIDS:::processBlockChr(fileReferenceGDS, fileLdBlock)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @keywords internal
processBlockChr <- function(fileReferenceGDS, fileBlock) {

    if (!(is.character(fileReferenceGDS) && (file.exists(fileReferenceGDS)))) {
        stop("The \'fileReferenceGDS\' must be a character string ",
             "representing the Reference GDS file. The file must exist.")
    }
    if (!(is.character(fileBlock) && (file.exists(fileBlock)))) {
        stop("The \'fileBlock\' must be a character string representing the",
                " file .det from plink block result. The file must exist.")
    }

    gdsReference <- snpgdsOpen(filename=fileReferenceGDS)
    blockChr <- read.delim(fileBlock, sep="")

    listChr <- unique(blockChr$CHR)
    if(length(listChr) != 1){
        stop(paste0("There is not just one CHR in ", fileBlock, "\n"))
    }
    listChr <- as.integer(gsub("chr", "", listChr))
    listSNVChr <- read.gdsn(index.gdsn(gdsReference, "snp.chromosome"))
    listSNVChr <- which(listSNVChr == listChr)
    snpKeep <- read.gdsn(index.gdsn(gdsReference, "snp.position"))[listSNVChr]
    closefn.gds(gdsReference)
    z <- cbind(c(blockChr$BP1, snpKeep, blockChr$BP2+1),
               c(seq_len(nrow(blockChr)),
                 rep(0, length(snpKeep)), -1*seq_len(nrow(blockChr))))

    z <- z[order(z[,1]),]
    blockSnp <- cumsum(z[,2])[z[,2] == 0]

    curStart <- 0
    activeBlock <- 0
    blockState <- 0
    blockInter <- rep(0, length(which(blockSnp == 0)))
    k <- 1
    for(i in seq_len(length(blockSnp))){
        if(blockSnp[i] == 0){
            if(activeBlock == 1){
                if(snpKeep[i] - curStart >= 10000) {
                    blockState <- blockState - 1

                    curStart <- snpKeep[i]
                }
            } else{
                blockState <- blockState - 1
                curStart <- snpKeep[i]
                activeBlock <- 1
            }
            if(blockState == 0){
                blockState <- -1
            }
            blockInter[k] <- blockState
            k <- k + 1
        }else{
            activeBlock <- 0
        }
    }

    blockSnp[blockSnp == 0] <- blockInter
    res <- list(chr=listChr,
                block.snp=blockSnp)
    return(res)
}


#' @title Read a VCF file with the genotypes use for the ancestry call
#'
#' @description The function reads VCF file and returns a data frame
#' containing the information about the read counts for the SNVs present in
#' the file.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{list}
#' \describe{
#' \item{matProfile}{
#' a \code{data.frame} containing at least:
#' \describe{
#' \item{Chromosome}{ a \code{numeric} representing the name of
#' the chromosome}
#' \item{Position}{ a \code{numeric} representing the position on the
#' chromosome}
#' \item{Ref}{ a \code{character} string representing the reference nucleotide}
#' \item{Alt}{ a \code{character} string representing the alternative
#' nucleotide}
#' \item{File1R}{ a \code{numeric} representing the count for
#' the reference nucleotide}
#' \item{File1A}{ a \code{numeric} representing the count for the
#' alternative nucleotide}
#' \item{count}{ a \code{numeric} representing the total count}
#' }}
#' \item{pRAIDS}{
#' a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#' }
#' }
#'
#' @examples
#'
#' ## Current directory
#' dataDir <- file.path(tempdir())
#'
#' fileGDS <- file.path(system.file("extdata/tests", package="RAIDS"), "ex1_good_small_1KG.gds")
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
#' pRAIDS <- RAIDSparam(profileFile=file.path(dataDir, "ex1.txt.gz"),
#'     studyDF=studyDF,
#'     genoSource="snp-pileup",
#'     pedStudy=samplePED,
#'     pathProfileGDS=dataDir,
#'     fileReferenceGDS=fileGDS
#'     )
#' infoCov <- RAIDS:::dispatchParseReadCov(pRAIDS)
#' names(infoCov)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @keywords internal
dispatchParseReadCov <- function(pRAIDS){

    matSample <- NULL
    if(genoSource(pRAIDS) == "snp-pileup") {
        matSample <- readSNVPileupFile(profileFile(pRAIDS), posOffset(pRAIDS))
    } else if(genoSource(pRAIDS) == "generic") {
        matSample <- readSNVFileGeneric(profileFile(pRAIDS), posOffset(pRAIDS))
    } else if(genoSource(pRAIDS) == "VCF") {
        # tmpProfile <- gsub(".vcf.gz", "",listMat[pos])
        matSample <- readSNVVCF(profileFile(pRAIDS),
                                profileName=pedStudy(pRAIDS)$Name.ID[1], posOffset(pRAIDS))
    } else if(genoSource(pRAIDS)== "bam"){
        colnames(listPos(pRAIDS))[seq_len(2)] <- c("chr", "start")
        matSample <- readSNVBAM(fileName=profileFile(pRAIDS),
                                varSelected=listPos(pRAIDS), paramSNVBAM=paramAncestry(pRAIDS),
                                posOffset(pRAIDS),
                                verbose=verbose(pRAIDS))
        # listPos <- do.call(rbind, listPos)
        tmp <- listPos(pRAIDS)
        colnames(tmp)[seq_len(2)] <- c("snp.chromosome", "snp.position")
        listPos(pRAIDS) <- tmp
        rm(tmp)

    }
    return(list(matProfile = matSample,
                pRAIDS=pRAIDS))
}

#' @title run ADMIXTURE for the synthetic data
#'
#' @description TODO
#'
#'
#' @param pathOut  a \code{character} representing path
#' where the directory for the synthetic admixture
#' will be processed
#'
#' @param gr  a \code{integer} representing the group
#' (row matGr from getMatrixPopSynthetic) of
#' synthetic data to run together
#'
#' @param matGr a \code{matrix} containing the sample identifiers and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @param matPIndex a \code{array} of the frequency matrix.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{0L} when successful.
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#'
#' ### runSynADMIXTURE in tools.R
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom BiocParallel MulticoreParam bplapply
#' @encoding UTF-8
#' @keywords internal
prepSynADMIXTURE <- function(pathOut, gr, matGr, matPIndex, pRAIDS){

    pathSyn <- file.path(pathOut, paste0("syn.", gr))
    if(!dir.exists(file.path(pathSyn))){
        dir.create(file.path(pathSyn))
    }
    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))

    admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
    pos <- admixDesc[admixDesc$id == anchorId(pRAIDS), "id.seq"]
    popSelected <- read.gdsn(index.gdsn(gdsReference, "admixture.ag"))
    popSelected <- popSelected[popSelected$id.seq == pos,]

    snpgdsClose(gdsReference)

    matGrDataId <- getMatrixDataId(matGr=matGr, pRAIDS=pRAIDS)

    writeBedRef(fileOut=file.path(pathSyn, paste0(pedStudy(pRAIDS)$Name.ID[1], ".synRef.bed")),
                listRM=matGr[gr,], subset=TRUE, pRAIDS=pRAIDS)
    # writeBedProfile(fileOut=file.path(pathOut, paste0(pRAIDS$pedStudy$Name.ID[1], ".syn.bed")),
    #                 listProfile=matGrDataId[i,], listRM=NULL, profileOnly=TRUE, pRAIDS=pRAIDS)
    writeFamRef(fileOut=file.path(pathSyn, paste0(pedStudy(pRAIDS)$Name.ID[1], ".synRef.fam")),
                listRM=matGr[gr,], subset=TRUE, pRAIDS=pRAIDS)

    writeFamProfile(fileOut=file.path(pathSyn, paste0(pedStudy(pRAIDS)$Name.ID[1], ".syn.fam")),
                    listProfile=matGrDataId[gr,], listRM=NULL, profileOnly=TRUE, pRAIDS)
    matPFile <- ifelse(matPIndex[gr] == 0,
                       file.path(pathOut, "matP", paste0( pedStudy(pRAIDS)$Name.ID[1],
                                                          ".P.in")),
                       file.path(pathOut, "matP", paste0( pedStudy(pRAIDS)$Name.ID[1],
                                                          ".syn.", matPIndex[gr], ".P.in")))
    file.copy(matPFile,
              file.path(pathSyn,
                        paste0(pedStudy(pRAIDS)$Name.ID[1],
                               ".synRef.", nrow(popSelected), ".P.in")))
    file.copy(file.path(pathOut, "profileAdmixture", paste0( pedStudy(pRAIDS)$Name.ID[1],
                                               ".bim")),
              file.path(pathSyn,
                        paste0(pedStudy(pRAIDS)$Name.ID[1],
                               ".synRef.bim")))
    writeBedBimFilesProfileFilter(pathOut=pathSyn,
                                  prefOut=paste0(pedStudy(pRAIDS)$Name.ID[1],
                                                 ".syn"),
                                  fileP=file.path(pathSyn,
                                        paste0(pedStudy(pRAIDS)$Name.ID[1],
                                        ".synRef.", nrow(popSelected), ".P.in")),
                                  listProfile=matGrDataId[gr,],
                                  pRAIDS=pRAIDS)
    return(0L)
}

#' @title run ADMIXTURE for the synthetic data
#'
#' @description TODO
#'
#'
#' @param pathOut  a \code{character} representing path
#' where the directory for the synthetic admixture
#' will be processed
#'
#' @param gr  a \code{integer} representing the group
#' (row matGr from getMatrixPopSynthetic) of
#' synthetic data to run together
#'
#' @param matGr a \code{matrix} containing the sample identifiers and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @param matPIndex a \code{array} of the frequency matrix.
#'
#' @param matQGS a \code{matrix} containing contening the admixture gold standard
#' from the gdsReference
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return \code{0L} when successful.
#'
#'
#' @details TODO
#'
#'
#' @examples
#'
#'
#' ### see runSynADMIXTURE in tools.R
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom BiocParallel MulticoreParam bplapply
#' @encoding UTF-8
#' @keywords internal
runGrSynADMIXTURE <- function(pathOut, gr, matGr, matPIndex, matQGS, pRAIDS) {

    if(verbose(pRAIDS)) print(paste0("Syn Admixture ", gr))

    prepSynADMIXTURE(pathOut=pathOut, gr=gr, matGr=matGr, matPIndex=matPIndex, pRAIDS)

    matQSynRef <- runADMIXTURE(pathOut = file.path(pathOut, paste0("syn.", gr)),
        fileBed = paste0(pedStudy(pRAIDS)$Name.ID[1], ".synRef.bed"),
        seedAdmix=seedADMIXTURE(pRAIDS) + gr, pRAIDS)

    matQSyn <- runADMIXTURE(pathOut = file.path(pathOut, paste0("syn.", gr)),
        fileBed = paste0(pedStudy(pRAIDS)$Name.ID[1], ".syn.bed"),
        seedAdmix=seedADMIXTURE(pRAIDS) + gr, pRAIDS)

    matQ <- data.frame(id=row.names(matQSynRef),
                idSyn=row.names(matQSyn),
            stringsAsFactors = FALSE)
    colnames(matQSynRef) <- paste0(colnames(matQSynRef), ".p")

    row.names(matQ) <- NULL
    row.names(matQSynRef) <- NULL

    matQ <- cbind(matQ, matQSyn, matQSynRef, matQGS[matQ$id,])
    row.names(matQ) <- matQ$id
    return(matQ)

}

#' @title Generate a plink fam file
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new fam file.
#'
#' @param listRM a \code{vector} of \code{character} strings containing the
#' identifiers for the reference samples that need to be removed for the
#' PCA analysis.
#'
#' @param subset a \code{logical} extract the profiles from the reference in
#' listRM (TRUE) or extract from the reference not in listRM (FALSE).
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## see runSynADMIXTURE in tools.R and
#' ## prepSynADMIXTURE in tools_internal.R
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom genio write_fam
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @encoding UTF-8
#' @keywords internal
writeFamRef <- function(fileOut, listRM=NULL, subset=FALSE, pRAIDS){

    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))

    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
    if(subset){
        if(! is.null(listRM)){
            sampleId <- sampleId[sampleId %in% listRM ]
        }
    } else{
        if(! is.null(listRM)){
            sampleId <- sampleId[! sampleId %in% listRM ]
        }
    }

    sujetFam <- data.frame(fam=seq_len(length(sampleId)),
                           id=c(sampleId),
                           pat=0,
                           mat=0,
                           sex=1,
                           pheno=0,
                           stringsAsFactors=FALSE)

    write_fam(fileOut, sujetFam)
    snpgdsClose(gdsReference)
    return(0)
}

#' @title Generate a plink bed file
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new bim file.
#'
#' @param listRM a \code{vector} of \code{character} strings containing the
#' identifiers for the reference samples that need to be removed for the
#' PCA analysis.
#'
#' @param subset a \code{logical} extract the genotypes from the reference in
#' listRM (TRUE) or extract from the reference not in listRM (FALSE).
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## see runSynADMIXTURE in tools.R and
#' ## prepSynADMIXTURE in tools_internal.R
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom genio write_bed
#' @importFrom SNPRelate snpgdsOpen snpgdsClose snpgdsGetGeno
#' @encoding UTF-8
#' @keywords internal
writeBedRef <- function(fileOut, listRM=NULL, subset=FALSE, pRAIDS){
    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                            paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))

    snpId <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[
                        read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
    if(subset){
        if(! is.null(listRM)){
            sampleId <- sampleId[sampleId %in% listRM]
        }
    } else{
        if(! is.null(listRM)){
            sampleId <- sampleId[! sampleId %in% listRM]
        }
    }
    g <- snpgdsGetGeno(gdsReference, sample.id=sampleId, snp.id=snpId,
            snpfirstdim=TRUE)

    write_bed(fileOut, g, verbose=TRUE, append=FALSE)
    snpgdsClose(gdsProfile)
    snpgdsClose(gdsReference)
    return(0)
}


#' @title Generate a matrix of synthetic group per pop
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param matGr a \code{matrix} containing the sample identifiers and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return a \code{matrix} containing the data.id for synthetic data and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @examples
#'
#' ## see runSynADMIXTURE in tools.R and
#' ## prepSynADMIXTURE in tools_internal.R
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @encoding UTF-8
#' @keywords internal
getMatrixDataId <- function(matGr, pRAIDS){
    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                            paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)

    profileSyn <- read.gdsn(index.gdsn(gdsProfile,"study.annot"))
    profileSyn <- profileSyn[profileSyn$study.id ==
                                        studyDFSyn(pRAIDS)$study.id[1],]
    sampleId <- read.gdsn(index.gdsn(gdsProfile,"sample.id"))

    listGr <- lapply(seq_len(nrow(matGr)),
                     FUN=function(x,matGr, profileSyn, sampleId){
                         return(sampleId[ sampleId %in% profileSyn$data.id[profileSyn$case.id %in%
                                                        matGr[x,]]])
                     },
                     matGr=matGr,
                     profileSyn=profileSyn,
                     sampleId=sampleId)
    matGrDataId <- do.call(rbind, listGr)

    snpgdsClose(gdsProfile)

    return(matGrDataId)
}


#' @title Generate a plink bed file TODO
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param pathOut a \code{character} string representing the path where
#' to save the output
#'
#' @param prefOut a \code{character} string representing the prefix of
#' file name of the new bed, bim, and P.in files.
#'
#' @param fileP a \code{character} string representing the path and
#' the file name of the new bim file.
#'
#' @param listProfile a \code{vector} of \code{character} representing the list
#' of profiles from gdsProfile to keep.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## see runSynADMIXTURE in tools.R and
#' ## prepSynADMIXTURE in tools_internal.R
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom genio write_bed write_bim
#' @importFrom methods is
#' @importFrom SNPRelate snpgdsGetGeno snpgdsOpen snpgdsClose
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @keywords internal
writeBedBimFilesProfileFilter <- function(pathOut, prefOut, fileP, listProfile,
    pRAIDS) {

    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                                    paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    snpId <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    g <- NULL

    #studyAn <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    gS <- snpgdsGetGeno(gdsProfile, sample.id=listProfile, snp.id=snpId,
                            snpfirstdim=TRUE)
    nbS <- ncol(gS)
    list2RM <- which(rowSums(is.na(gS))/nbS > 0.25)
    if(length(list2RM) > 0){
        gS <- gS[-1 * list2RM,]
    }

    snpInfo <- data.frame(chr=read.gdsn(index.gdsn(gdsProfile,
                                                "snp.chromosome")),
                          id=read.gdsn(index.gdsn(gdsProfile, "snp.id")),
                          posg=0,
                          pos=read.gdsn(index.gdsn(gdsProfile, "snp.position")),
                          alt="A",
                          ref="B",
                          stringsAsFactors=FALSE)
    if(length(list2RM) > 0){
        snpInfo <- snpInfo[-1 * list2RM,]
    }
    write_bim(file.path(pathOut, paste0(prefOut,
        ".bim")), snpInfo)
    write_bed(file.path(pathOut, paste0(prefOut,
        ".bed")), gS, verbose=TRUE, append=FALSE)

    snpgdsClose(gdsProfile)

    resP <- read.delim(fileP, sep=" ", header=FALSE)
    if(length(list2RM) > 0){
        resP <- resP[-1 * list2RM,]
    }
    # 068ba2ae-288c-446d-8d17-72445ce4f788.syn.1.sv.5.P.in
    write.table(resP, file.path(pathOut, paste0(prefOut ,
        ".", ncol(resP), ".P.in")),
                col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")
    # write.table(resP, file.path(pathOut, paste0(pRAIDS$pedStudy$Name.ID[1],
    #                                     ".syn.", indexS,".sv.5.P.in")),
    #             col.names=FALSE, row.names=FALSE, quote=FALSE, sep=" ")
    return(0)
    }

#' @title Generate a VCF with the information from the SNPs that pass
#' a cut-off threshold
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param gds an object of class
#' \code{\link[SNPRelate:SNPGDSFileClass]{SNPRelate::SNPGDSFileClass}}, a SNP
#' GDS file.
#'
#' @param fileOUT a \code{character} string representing the path and file
#' name of the VCF file that will be created wit the retained SNP information.
#'
#' @param offset a single \code{integer} that is added to the SNP position to
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' Default: \code{0L}.
#'
#' @param freqCutoff a single positive \code{numeric} specifying the cut-off to
#' keep a SNP. If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' ## Path to the demo pedigree file is located in this package
#' data.dir <- system.file("extdata", package="aicsPaper")
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom methods is
#' @importFrom S4Vectors isSingleInteger
#' @encoding UTF-8
#' @export
snvListVCF <- function(gds, fileOUT, offset=0L, freqCutoff=NULL) {

    ## Validate that gds is an object of class SNPGDSFileClass
    if (! is(gds, "SNPGDSFileClass")) {
        stop("The \'gds\' parameter must be an object of ",
             "class \'SNPGDSFileClass\'.")
    }

    ## Validate that offset is a single integer
    if (! isSingleInteger(offset)) {
        stop("The \'offset\' must be a single integer.")
    }

    ## Validate that freqCutoff is a single numeric or NULL
    if (! isSingleNumber(freqCutoff) & ! is.null(freqCutoff)) {
        stop("The \'freqCutoff\' must be a single numeric or NULL.")
    }

    snp.chromosome <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
    snp.position <- read.gdsn(index.gdsn(gds, "snp.position"))
    snp.allele <- read.gdsn(index.gdsn(gds, "snp.allele"))

    allele <- matrix(unlist(strsplit(snp.allele, "\\/")), nrow=2)

    df <- NULL

    if(is.null(freqCutoff)){
        snp.AF = read.gdsn(index.gdsn(gds, "snp.AF"))
        df <- data.frame(CHROM=snp.chromosome,
                            POS=as.integer(snp.position + offset),
                            ID=rep(".", length(snp.chromosome)),
                            REF=allele[1,],
                            ALT=allele[2,],
                            QUAL=rep(".", length(snp.chromosome)),
                            FILTER=rep(".", length(snp.chromosome)),
                            INFO=paste0("AF=", snp.AF),
                            stringsAsFactors=FALSE)
    } else{
        freqDF <- data.frame(snp.AF=read.gdsn(index.gdsn(gds, "snp.AF")),
                        snp.EAS_AF=read.gdsn(index.gdsn(gds, "snp.EAS_AF")),
                        snp.EUR_AF=read.gdsn(index.gdsn(gds, "snp.EUR_AF")),
                        snp.AFR_AF=read.gdsn(index.gdsn(gds, "snp.AFR_AF")),
                        snp.AMR_AF=read.gdsn(index.gdsn(gds, "snp.AMR_AF")),
                        snp.SAS_AF=read.gdsn(index.gdsn(gds, "snp.SAS_AF")))

        listKeep <- which(rowSums(freqDF[,2:6] >= freqCutoff & freqDF[,2:6] <= 1 - freqCutoff) > 0)
        df <- data.frame(CHROM=snp.chromosome[listKeep],
                            POS=as.integer(snp.position[listKeep] + offset),
                            ID=rep(".", length(listKeep)),
                            REF=allele[1,listKeep],
                            ALT=allele[2,listKeep],
                            QUAL=rep(".", length(listKeep)),
                            FILTER=rep(".", length(listKeep)),
                            INFO=paste0("AF=", freqDF$snp.AF[listKeep]),
                            stringsAsFactors=FALSE)
    }

    ## Add the header
    ##fileformat=VCFv4.3
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    cat(paste0('##fileformat=VCFv4.3', "\n"), file = fileOUT)
    cat(paste0('##FILTER=<ID=PASS,Description="All filters passed">',
               "\n"), file = fileOUT, append=TRUE)
    cat(paste0('##INFO=<ID=AF,Number=A,Type=Float,',
                'Description="Estimated allele frequency in the range (0,1)">',
                "\n"), file = fileOUT, append=TRUE)
    cat('#', file = fileOUT, append=TRUE)

    write.table(df, file=fileOUT, sep="\t",
                    append=TRUE, row.names=FALSE,
                    col.names=TRUE, quote=FALSE)

    return(0L)
}

#' @title Merge the pruning files by chromosome in one file
#'
#' @description TODO
#'
#' @param PATHPRUNED TODO
#'
#' @param filePref TODO
#'
#' @param fileOUT TODO
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alexander Krasnitz
#'
#' @keywords internal
groupChrPruning <- function(PATHPRUNED, filePref, fileOUT) {

    prunedList <- list()

    # Read the content of each file (one file per chromosome)
    for(i in seq_len(22)) {
        fileChr <- file.path(PATHPRUNED, paste0(filePref, i, ".rds"))

        if(file.exists(fileChr)) {
            prunedList[[i]] <- readRDS(fileChr)
        } else {
            stop("Problem with the file: ", fileChr)
        }
    }

    ## Merge the content of all files
    pruned <- do.call(c, prunedList)

    ## Save all the information into one file
    saveRDS(pruned, fileChr <- file.path(PATHPRUNED, fileOUT))
}


#' @title Merge the pruning files by chromosome in one file
#'
#' @description TODO
#'
#' @param PATHGENOCHR TODO
#'
#' @param PATHOUT TODO
#'
#' @return The integer \code{0} when successful.
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv2 read.csv2
#' @encoding UTF-8
#' @export
groupChr1KGSNV <- function(PATHGENOCHR, PATHOUT) {

    listFiles <- dir(file.path(PATHGENOCHR, "chr1"), ".+\\.chr1\\.vcf\\.bz2")
    listSamples <- gsub("\\.chr1\\.vcf\\.bz2", "", listFiles)

    for(sampleId in listSamples){
        listGeno <- list()

        for(chr in seq_len(22)){

            geno <- read.csv2(file.path(PATHGENOCHR, paste0("chr", chr),
                                        paste0(sampleId,
                                               ".chr", chr,".vcf.bz2") ),
                              sep="\t",
                              row.names = NULL)

            listGeno[[paste0("chr", chr)]] <- geno
        }
        genoAll <- do.call(rbind, listGeno)
        write.csv2(genoAll, file=bzfile(file.path(PATHOUT,
                            paste0(sampleId, ".csv.bz2"))), row.names=FALSE)
    }

    return(0L)
}


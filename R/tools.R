#' @title Generate a vcf for the SNP with at least freqCutoff for all frequency
#' for at least one super population
#'
#' @description TODO
#'
#' @param PATHPRUNED
#'
#' @param filePref
#'
#' @param fileOUT
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
#'
#' @keywords internal

snvListVCF <- function(gds,
                       fileOUT,
                       offset = 0,
                       freqCutoff = NULL){

    snp.chromosome <- read.gdsn(index.gdsn(gds, "snp.chromosome"))
    snp.position <- read.gdsn(index.gdsn(gds, "snp.position"))
    snp.allele <- read.gdsn(index.gdsn(gds, "snp.allele"))

    allele <- matrix(unlist(strsplit(snp.allele, "\\/")), nr=2)

    df <- NULL

    if(is.null(freqCutoff)){
        snp.AF = read.gdsn(index.gdsn(gds, "snp.AF"))
        df <- data.frame(CHROM = snp.chromosome,
                         POS = as.integer(snp.position + offset),
                         ID = rep(".", length(snp.chromosome)),
                         REF = allele[1,],
                         ALT = allele[2,],
                         QUAL = rep(".", length(snp.chromosome)),
                         FILTER = rep(".", length(snp.chromosome)),
                         INFO = paste0("AF=", snp.AF),
                         stringsAsFactors = FALSE)
    } else{
        freqDF <- data.frame(snp.AF = read.gdsn(index.gdsn(gds, "snp.AF")),
                             snp.EAS_AF = read.gdsn(index.gdsn(gds, "snp.EAS_AF")),
                             snp.EUR_AF = read.gdsn(index.gdsn(gds, "snp.EUR_AF")),
                             snp.AFR_AF = read.gdsn(index.gdsn(gds, "snp.AFR_AF")),
                             snp.AMR_AF = read.gdsn(index.gdsn(gds, "snp.AMR_AF")),
                             snp.SAS_AF = read.gdsn(index.gdsn(gds, "snp.SAS_AF")))

        listKeep <- which(rowSums(freqDF[,2:6] >= freqCutoff & freqDF[,2:6] <= 1 - freqCutoff) > 0)
        df <- data.frame(CHROM = snp.chromosome[listKeep],
                         POS = as.integer(snp.position[listKeep] + offset),
                         ID = rep(".", length(listKeep)),
                         REF = allele[1,listKeep],
                         ALT = allele[2,listKeep],
                         QUAL = rep(".", length(listKeep)),
                         FILTER = rep(".", length(listKeep)),
                         INFO = paste0("AF=", freqDF$snp.AF[listKeep]),
                         stringsAsFactors = FALSE)
    }




    ## Add the header
    ##fileformat=VCFv4.3
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    cat(paste0('##fileformat=VCFv4.3', "\n"), file = fileOUT)
    cat(paste0('##FILTER=<ID=PASS,Description="All filters passed">',
               "\n"),
        file = fileOUT, append = TRUE)
    cat(paste0('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">',
               "\n"),
        file = fileOUT, append = TRUE)
    cat('#',
        file = fileOUT, append = TRUE)

    write.table(df, file = fileOUT, sep = "\t",
                append = TRUE,
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE)



}

#' @title Merge the pruning files by chromosome in one file
#'
#' @description TODO
#'
#' @param PATHPRUNED
#'
#' @param filePref
#'
#' @param fileOUT
#'
#' @return TODO a \code{vector} of \code{numeric}
#'
#' @examples
#'
#' # TODO
#'
#' @author Pascal Belleau, Astrid Desch&ecirc;nes and Alex Krasnitz
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


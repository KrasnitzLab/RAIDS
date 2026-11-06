#' @title Generate a VCF with the information from the SNPs that pass
#' a cut-off threshold
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param gdsReference an object of class \code{\link[gdsfmt]{gds.class}}
#' (a GDS file), the 1KG GDS file.
#'
#' @param fileOut a \code{character} string representing the path and file
#' name of the VCF file that will be created wit the retained SNP information.
#' The file should have the ".vcf" extension.
#'
#' @param offset a single \code{integer} that is added to the SNP position to
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' Default: \code{0L}.
#'
#' @param freqCutoff a single positive \code{numeric} specifying the cut-off to
#' keep a SNP. If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo 1KG Reference GDS file
#' fileGDS <- openfn.gds(file.path(dataDir,
#'                     "PopulationReferenceDemo.gds"))
#'
#' ## Output VCF file that will be created (temporary)
#' vcfFile <- file.path(tempdir(), "Demo_TMP_01.vcf")
#'
#' ## Create a VCF file with the SNV dataset present in the GDS file
#' ## No cutoff on frequency, so all SNVs are saved
#' snvListVCF(gdsReference=fileGDS, fileOut=vcfFile, offset=0L,
#'                     freqCutoff=NULL)
#'
#' ## Close GDS file (IMPORTANT)
#' closefn.gds(fileGDS)
#'
#' ## Remove temporary VCF file
#' unlink(vcfFile, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn ls.gdsn
#' @importFrom methods is
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
snvListVCF <- function(gdsReference, fileOut, offset=0L, freqCutoff=NULL) {

    ## Validate that gdsReference is an object of class gds.class
    if (!inherits(gdsReference, "gds.class")) {
        stop("The \'gdsReference\' must be an object of class \'gds.class\'.")
    }

    ## Validate that offset is a single integer
    if (! isSingleNumber(offset)) {
        stop("The \'offset\' must be a single integer.")
    }

    ## Validate that freqCutoff is a single numeric or NULL
    if (! isSingleNumber(freqCutoff) && ! is.null(freqCutoff)) {
        stop("The \'freqCutoff\' must be a single numeric or NULL.")
    }

    snpChromosome <- read.gdsn(index.gdsn(gdsReference, "snp.chromosome"))
    snpPosition <- read.gdsn(index.gdsn(gdsReference, "snp.position"))
    snpAllele <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))

    allele <- matrix(unlist(strsplit(snpAllele, "\\/")), nrow=2)

    df <- NULL

    if(is.null(freqCutoff)){
        snpAF <- read.gdsn(index.gdsn(gdsReference, "snp.AF"))
        df <- data.frame(CHROM=snpChromosome,
                            POS=as.integer(snpPosition + offset),
                            ID=rep(".", length(snpChromosome)),
                            REF=allele[1,],
                            ALT=allele[2,],
                            QUAL=rep(".", length(snpChromosome)),
                            FILTER=rep(".", length(snpChromosome)),
                            INFO=paste0("AF=", snpAF),
                            stringsAsFactors=FALSE)
    } else {
        if( length(which(paste0("snp.",
                        c("EAS", "EUR", "AFR", "AMR", "SAS"),
                                "_AF") %in% ls.gdsn(gdsReference))) == 5){
            freqDF <- data.frame(
                snp.AF=read.gdsn(index.gdsn(gdsReference, "snp.AF")),
                snp.EAS_AF=read.gdsn(index.gdsn(gdsReference, "snp.EAS_AF")),
                snp.EUR_AF=read.gdsn(index.gdsn(gdsReference, "snp.EUR_AF")),
                snp.AFR_AF=read.gdsn(index.gdsn(gdsReference, "snp.AFR_AF")),
                snp.AMR_AF=read.gdsn(index.gdsn(gdsReference, "snp.AMR_AF")),
                snp.SAS_AF=read.gdsn(index.gdsn(gdsReference, "snp.SAS_AF")))

            listKeep <- which(rowSums(freqDF[,2:6] >= freqCutoff &
                                            freqDF[,2:6] <= 1 - freqCutoff) > 0)
            df <- data.frame(CHROM=snpChromosome[listKeep],
                                POS=as.integer(snpPosition[listKeep] + offset),
                                ID=rep(".", length(listKeep)),
                                REF=allele[1,listKeep],
                                ALT=allele[2,listKeep],
                                QUAL=rep(".", length(listKeep)),
                                FILTER=rep(".", length(listKeep)),
                                INFO=paste0("AF=", freqDF$snp.AF[listKeep]),
                                stringsAsFactors=FALSE)
        } else if("AF.superPop" %in% ls.gdsn(gdsReference)) {

            matAF <- index.gdsn(gdsReference, "AF.superPop")
            snpAF <- read.gdsn(matAF)
            freqDF <- cbind(read.gdsn(index.gdsn(gdsReference, "snp.AF")),snpAF)
            rm(snpAF)
            listKeep <- which(rowSums(freqDF[,2:ncol(freqDF)] >= freqCutoff &
                            freqDF[,2:ncol(freqDF)] <= 1 - freqCutoff) > 0)
            df <- data.frame(CHROM=snpChromosome[listKeep],
                                POS=as.integer(snpPosition[listKeep] + offset),
                                ID=rep(".", length(listKeep)),
                                REF=allele[1,listKeep],
                                ALT=allele[2,listKeep],
                                QUAL=rep(".", length(listKeep)),
                                FILTER=rep(".", length(listKeep)),
                                INFO=paste0("AF=", freqDF[listKeep,1]),
                                stringsAsFactors=FALSE)
        }
    }
    df <- df[df$REF %in% c("A", "C", "G", "T") & df$ALT %in%
                    c("A", "C", "G", "T"),]

    df$CHROM <- paste0("chr", df$CHROM)

    ## Add the header
    ##fileformat=VCFv4.3
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency
    ##                            in the range (0,1)">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    cat(paste0('##fileformat=VCFv4.3', "\n"), file=fileOut)
    cat(paste0('##FILTER=<ID=PASS,Description="All filters passed">',
                "\n"), file=fileOut, append=TRUE)
    cat(paste0('##INFO=<ID=AF,Number=A,Type=Float,',
                'Description="Estimated allele frequency in the range (0,1)">',
                "\n"), file=fileOut, append=TRUE)
    cat('#', file=fileOut, append=TRUE)

    write.table(df, file=fileOut, sep="\t", append=TRUE, row.names=FALSE,
                    col.names=TRUE, quote=FALSE)

    ## Successful
    return(0L)
}

#' @title Generate a VCF with the information from the SNPs that pass
#' a cut-off threshold
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param pRAIDS a TODO.
#'
#' @param fileOut a \code{character} string representing the path and file
#' name of the VCF file that will be created wit the retained SNP information.
#' The file should have the ".vcf" extension.
#'
#' @param freqCutoff a single positive \code{numeric} specifying the cut-off to
#' keep a SNP. If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## TODO
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn ls.gdsn openfn.gds closefn.gds
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
snvListVCFRef <- function(pRAIDS, fileOut, freqCutoff=NULL) {
    # gdsReference, fileOut, offset=0L, freqCutoff=NULL
    ## Validate that gdsReference is an object of class gds.class
    # if (!inherits(gdsReference, "gds.class")) {
    #     stop("The \'gdsReference\' must be an object of class \'gds.class\'.")
    # }
    #
    # ## Validate that offset is a single integer
    # if (! isSingleNumber(offset)) {
    #     stop("The \'offset\' must be a single integer.")
    # }
    #
    # ## Validate that freqCutoff is a single numeric or NULL
    # if (! isSingleNumber(freqCutoff) && ! is.null(freqCutoff)) {
    #     stop("The \'freqCutoff\' must be a single numeric or NULL.")
    # }
    # if reference HGDP1kg no offset

    gdsReference <- openfn.gds(fileReferenceGDS(pRAIDS))

    snpChromosome <- read.gdsn(index.gdsn(gdsReference, "snp.chromosome"))
    snpPosition <- read.gdsn(index.gdsn(gdsReference, "snp.position"))
    snpAllele <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))
    listSNPKepp <- NULL
    allele <- matrix(unlist(strsplit(snpAllele, "\\/")), nrow=2)
    if("snp.KeepDefault" %in% ls.gdsn(gdsReference)){
        listSNPKepp <- read.gdsn(index.gdsn(gdsReference, "snp.KeepDefault"))
        snpChromosome <- snpChromosome[listSNPKepp]
        snpPosition <- snpPosition[listSNPKepp]
        allele <- allele[,listSNPKepp]
    }

    df <- NULL

    if(is.null(freqCutoff)){
        snpAF <- read.gdsn(index.gdsn(gdsReference, "snp.AF"))
        if(! is.null(listSNPKepp)){
            snpAF <- snpAF[listSNPKepp]
        }
        df <- data.frame(CHROM=snpChromosome,
                            POS=as.integer(snpPosition + posOffset(pRAIDS)),
                            ID=rep(".", length(snpChromosome)),
                            REF=allele[1,],
                            ALT=allele[2,],
                            QUAL=rep(".", length(snpChromosome)),
                            FILTER=rep(".", length(snpChromosome)),
                            INFO=paste0("AF=", snpAF),
                            stringsAsFactors=FALSE)
    } else {
        if( length(which(paste0("snp.",
                                c("EAS", "EUR",
                                  "AFR", "AMR",
                                  "SAS"),
                                "_AF") %in% ls.gdsn(gdsReference))) == 5){
            freqDF <- data.frame(
                snp.AF=read.gdsn(index.gdsn(gdsReference, "snp.AF")),
                snp.EAS_AF=read.gdsn(index.gdsn(gdsReference, "snp.EAS_AF")),
                snp.EUR_AF=read.gdsn(index.gdsn(gdsReference, "snp.EUR_AF")),
                snp.AFR_AF=read.gdsn(index.gdsn(gdsReference, "snp.AFR_AF")),
                snp.AMR_AF=read.gdsn(index.gdsn(gdsReference, "snp.AMR_AF")),
                snp.SAS_AF=read.gdsn(index.gdsn(gdsReference, "snp.SAS_AF")))

            listKeep <- which(rowSums(freqDF[,2:6] >= freqCutoff &
                                          freqDF[,2:6] <= 1 - freqCutoff) > 0)
            df <- data.frame(CHROM=snpChromosome[listKeep],
                             POS=as.integer(snpPosition[listKeep] + posOffset(pRAIDS)),
                             ID=rep(".", length(listKeep)),
                             REF=allele[1,listKeep],
                             ALT=allele[2,listKeep],
                             QUAL=rep(".", length(listKeep)),
                             FILTER=rep(".", length(listKeep)),
                             INFO=paste0("AF=", freqDF$snp.AF[listKeep]),
                             stringsAsFactors=FALSE)
        }else if("AF.superPop" %in% ls.gdsn(gdsReference)){

            matAF <- index.gdsn(gdsReference, "AF.superPop")
            snpAF <- read.gdsn(matAF)
            freqDF <- cbind(read.gdsn(index.gdsn(gdsReference, "snp.AF")),snpAF)
            rm(snpAF)
            listKeep <- which(rowSums(freqDF[,2:ncol(freqDF)] >= freqCutoff &
                                freqDF[,2:ncol(freqDF)] <= 1 - freqCutoff) > 0)
            df <- data.frame(CHROM=snpChromosome[listKeep],
                        POS=as.integer(snpPosition[listKeep] + posOffset(pRAIDS)),
                        ID=rep(".", length(listKeep)),
                        REF=allele[1,listKeep],
                        ALT=allele[2,listKeep],
                        QUAL=rep(".", length(listKeep)),
                        FILTER=rep(".", length(listKeep)),
                        INFO=paste0("AF=", freqDF[listKeep,1]),
                        stringsAsFactors=FALSE)
        }
    }
    df <- df[df$REF %in% c("A", "C", "G", "T") & df$ALT %in%
                                c("A", "C", "G", "T"),]
    df$CHROM <- paste0("chr", df$CHROM)

    closefn.gds(gdsReference)

    ## Add the header
    ##fileformat=VCFv4.3
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency
    ##                            in the range (0,1)">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    cat(paste0('##fileformat=VCFv4.3', "\n"), file=fileOut)
    cat(paste0('##FILTER=<ID=PASS,Description="All filters passed">',
               "\n"), file=fileOut, append=TRUE)
    cat(paste0('##INFO=<ID=AF,Number=A,Type=Float,',
               'Description="Estimated allele frequency in the range (0,1)">',
               "\n"), file=fileOut, append=TRUE)
    cat('#', file=fileOut, append=TRUE)

    write.table(df, file=fileOut, sep="\t", append=TRUE, row.names=FALSE,
                col.names=TRUE, quote=FALSE)

    ## Successful
    return(0L)
}

#' @title Merge the genotyping files per chromosome into one file
#'
#' @description This function merge all the genotyping files associated to one
#' specific sample into one file. That merged VCF file will be saved in a
#' specified directory and will have the name of the sample. It will also be
#' compressed (bzip). The function will merge the
#' files for all samples present in the input directory.
#'
#' @param pathGenoChr a \code{character} string representing the path where
#' the genotyping files for each sample and chromosome are located. The path
#' must contains sub-directories (one per chromosome) and the genotyping files
#' must be present in those sub-directories.
#' The path must exists.
#'
#' @param pathOut a \code{character} string representing the path where
#' the merged genotyping files for each sample will be created.
#' The path must exists.
#'
#' @return The integer \code{0L} when successful or \code{FALSE} if not.
#'
#' @examples
#'
#' ## Path to the demo vcf files in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#' pathGenoTar <- file.path(dataDir, "demoGenoChr", "demoGenoChr.tar")
#'
#' ## Path where the chromosomes files will be located
#' pathGeno <- file.path(tempdir(), "tempGeno")
#' dir.create(pathGeno, showWarnings=FALSE)
#'
#' ## Untar the file that contains the VCF files for 3 samples split by
#' ## chromosome (one directory per chromosome)
#' untar(tarfile=pathGenoTar, exdir=pathGeno)
#'
#' ## Path where the output VCF file will be created is
#' ## the same where the split VCF are (pathGeno)
#'
#' ## The files must not exist
#' if (!file.exists(file.path(pathGeno, "NA12003.csv.bz2")) &&
#'         !file.exists(file.path(pathGeno, "NA12004.csv.bz2")) &&
#'         !file.exists(file.path(pathGeno, "NA12005.csv.bz2"))) {
#'
#'         ## Return 0 when successful
#'         ## The files "NA12003.csv.bz2", "NA12004.csv.bz2" and
#'         ## "NA12005.csv.bz2" should not be present in the current directory
#'         groupChr1KGSNV(pathGenoChr=pathGeno, pathOut=pathGeno)
#'
#'         ## Validate that files have been created
#'         file.exists(file.path(pathGeno, "NA12003.csv.bz2"))
#'         file.exists(file.path(pathGeno, "NA12004.csv.bz2"))
#'         file.exists(file.path(pathGeno, "NA12005.csv.bz2"))
#'
#' }
#'
#' ## Remove temporary directory
#' unlink(pathGeno, recursive=TRUE, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom utils write.csv2 read.csv2
#' @encoding UTF-8
#' @export
groupChr1KGSNV <- function(pathGenoChr, pathOut) {

    ## Validate that the input path for the genotyping files exists
    if (! file.exists(pathGenoChr)) {
        stop("The path \'", pathGenoChr, "\' does not exist.")
    }

    ## Validate that the output path for the genotyping files exists
    if (! file.exists(pathOut)) {
        stop("The path \'", pathOut, "\' does not exist.")
    }

    ## Obtain the comprehensive list of samples
    listFiles <- dir(file.path(pathGenoChr, "chr1"), ".+\\.chr1\\.vcf\\.bz2")
    listSamples <- gsub("\\.chr1\\.vcf\\.bz2", "", listFiles)

    ## Merge files associated to each samples into one csv file
    results <- lapply(X=listSamples, FUN=function(sampleId, pathOut) {

        ## For each chromosome, read genotyping file and append the information
        listGeno <- lapply(seq_len(22), function(chr, sampleId) {
            geno <- read.csv2(file.path(pathGenoChr, paste0("chr", chr),
                        paste0(sampleId, ".chr", chr,".vcf.bz2")),
                            sep="\t", row.names=NULL)
            return(geno)}, sampleId=sampleId)

        genoAll <- do.call(rbind, listGeno)

        ## Save the genotyping information into one file
        write.csv2(genoAll, file=bzfile(file.path(pathOut,
            paste0(sampleId, ".csv.bz2"))), row.names=FALSE)

        return(TRUE)}, pathOut=pathOut)

    ## Successful or not
    return(ifelse(all(unlist(results)), 0L, FALSE))
}


#' @title Generate a VCF with the information from the SNPs that pass
#' a cut-off threshold
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param freqCutoff a single positive \code{numeric} specifying the cut-off to
#' keep a SNP. If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @param chr a \code{integer} between 1 to 22
#' chr1:22 specify SNP retained are on chr.
#' If \code{NULL}, all SNPs are retained. Default: \code{NULL}.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new VCF.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Path to the demo pedigree file is located in this package
#' dataDir <- system.file("extdata", package="RAIDS")
#'
#' ## Demo 1KG Reference GDS file
#' fileGDS <- openfn.gds(file.path(dataDir,
#'                     "PopulationReferenceDemo.gds"))
#'
#' ## Output VCF file that will be created (temporary)
#' vcfFile <- file.path(tempdir(), "Demo_TMP_01.vcf")
#'
#' ## Create a VCF file with the SNV dataset present in the GDS file
#' ## No cutoff on frequency, so all SNVs are saved
#' snvListVCF(gdsReference=fileGDS, fileOut=vcfFile, offset=0L,
#'                     freqCutoff=NULL)
#'
#' ## Close GDS file (IMPORTANT)
#' closefn.gds(fileGDS)
#'
#' ## Remove temporary VCF file
#' unlink(vcfFile, force=TRUE)
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom methods is
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
generateVCF <- function(freqCutoff=NULL, chr=NULL , fileOut, pRAIDS) {
    # Open GDS Reference file
    gdsReference <- snpgdsOpen(filename=fileReferenceGDS(pRAIDS))

    ## Validate that freqCutoff is a single numeric or NULL
    if (! isSingleNumber(freqCutoff) && ! is.null(freqCutoff)) {
        stop("The \'freqCutoff\' must be a single numeric or NULL.")
    }
    if(! is.null(chr)){
        chr <- as.integer(gsub("chr", "", chr))
        if(is.na(chr)){
            stop(paste0("The \'chr\':", chr, " is not an integer or chr1-chr22\n"))
        }
    }

    snpChromosome <- read.gdsn(index.gdsn(gdsReference, "snp.chromosome"))
    snpPosition <- read.gdsn(index.gdsn(gdsReference, "snp.position"))
    snpAllele <- read.gdsn(index.gdsn(gdsReference, "snp.allele"))

    allele <- matrix(unlist(strsplit(snpAllele, "\\/")), nrow=2)
    listSNP <- seq_len(length(snpPosition))
    if(!is.null(chr)){
        listSNP <- listSNP[snpChromosome == chr]
    }
    if(! is.null(freqCutoff)){
        snpAF <- read.gdsn(index.gdsn(gdsReference, "snp.AF"))[listSNP]
        listSNP <- listSNP[snpAF > freqCutoff &  snpAF < 1-freqCutoff]
    }
    closefn.gds(gdsReference)
    snpVCF <- data.frame(CHROM=snpChromosome[listSNP],
                         # not - the offset because if you load the vcf is + offset
                         POS=as.integer(snpPosition[listSNP] - posOffset(pRAIDS)),
                         ID = rep(".", length(listSNP)),
                         REF=allele[1,listSNP],
                         ALT=allele[2,listSNP],
                         QUAL = rep(".", length(listSNP)),
                         FILTER = rep("PASS", length(listSNP)),
                         INFO = rep(".", length(listSNP)),
                         FORMAT = rep("GT", length(listSNP)),
                         stringsAsFactors = FALSE)

    ## Profile GDS file name
    fileProfileGDS <-  validateProfileGDSExist(pathProfile=pathProfileGDS(pRAIDS),
                                               profile=pedStudy(pRAIDS)$Name.ID[1])
    gdsProfile <- openfn.gds(fileProfileGDS)

    listGeno <- list()
    g <- read.gdsn(index.gdsn(gdsProfile, "geno.ref"))[listSNP]
    closefn.gds(gdsProfile)
    geno <- rep("./.", length(listSNP))
    geno[which(g==0)] <- "0/0"
    geno[which(g==1)] <- "0/1"
    geno[which(g==2)] <- "1/1"
    snpVCF$QUAL[which(g == 3)] <- "."

    ## Add the header
    ##fileformat=VCFv4.3
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency
    ##                            in the range (0,1)">
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

    cat(paste0('##fileformat=VCFv4.3', "\n"), file = fileOut)
    cat(paste0('##FILTER=<ID=PASS,Description="All filters passed">',
               "\n"), file = fileOut, append=TRUE)
    cat(paste0('##INFO=<ID=GT,Number=1,Type=String,',
               'Description="Phased Genotype">',
               "\n"), file = fileOut, append=TRUE)
    listColName <- c('#CHROM', "POS", "ID",
                     "REF", "ALT", "QUAL",
                     "FILTER", "INFO","FORMAT",
                     pedStudy(pRAIDS)$Name.ID[1])
    cat(paste0(paste(listColName,
                     collapse = "\t"), "\n"),
        file = fileOut, append=TRUE)

    snpgVCF <- cbind(snpVCF, geno)

    snpgVCF$CHROM <- paste0("chr", snpgVCF$CHROM)
    write.table(x = snpgVCF, file = fileOut,
                quote = FALSE, sep = "\t",
                col.names = FALSE, row.names = FALSE,
                append = TRUE)


    ## Successful
    return(0L)
}

#' @title Generate a plink bim file
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new bim file.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' pathOut <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathOut)){
#'     dir.create(pathOut)
#' }
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#'
#' writeBimPruned(file.path(pathOut, "test.bim"), pRAIDS)
#'
#' dir(file.path(pathOut))
#'
#' unlink(file.path(pathOut, "test.bim"))
#' unlink(file.path(pathOut), recursive=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom genio write_bim
#' @importFrom methods is
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
writeBimPruned <- function(fileOut, pRAIDS){
    fileGDSProofile <- file.path(pathProfileGDS(pRAIDS), paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- openfn.gds(fileGDSProofile)

    snpInfo <- data.frame(chr = read.gdsn(index.gdsn(gdsProfile, "snp.chromosome")),
                          id = read.gdsn(index.gdsn(gdsProfile, "snp.id")),
                          posg = 0,
                          pos = read.gdsn(index.gdsn(gdsProfile, "snp.position")),
                          alt = "A",
                          ref = "B",
                          stringsAsFactors = FALSE)
    closefn.gds(gdsProfile)
    write_bim(fileOut, snpInfo)
    return(0)
}


#' @title Generate a plink bed file TODO
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param fileOut a \code{character} string representing the path and
#' the file name of the new bim file.
#'
#' @param listProfile a \code{vector} of \code{character} representing the list
#' of profile from gdsProfile to keep
#'
#' @param listRM a \code{list} TODO. Default: \code{NULL}.
#'
#' @param profileOnly a \code{boolean} TODO. Default: \code{FALSE}.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' pathOut <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathOut)){
#'     dir.create(pathOut)
#' }
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#'
#' writeBedProfile(fileOut=file.path(pathOut,
#'     paste0(pedStudy(pRAIDS)$Name.ID[1], ".bed")),
#'     listProfile=pedStudy(pRAIDS)$Name.ID[1],
#'     profileOnly = TRUE, pRAIDS=pRAIDS)
#'
#' dir(file.path(pathOut))
#'
#' unlink(file.path(pathOut, paste0(pedStudy(pRAIDS)$Name.ID[1], ".bed")))
#' unlink(file.path(pathOut), recursive=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom genio write_bed
#' @importFrom methods is
#' @importFrom SNPRelate snpgdsGetGeno snpgdsClose snpgdsOpen
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
writeBedProfile <- function(fileOut, listProfile, listRM=NULL,
    profileOnly=FALSE, pRAIDS) {

    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                                    paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    snpId <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
    g <- NULL
    if(! profileOnly){
        gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))

        sampleId <- read.gdsn(index.gdsn(gdsReference,
                    "sample.id"))[
                        read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
        if(! is.null(listRM)){
            sampleId <- sampleId[! sampleId %in% listRM]
        }
        g <- snpgdsGetGeno(gdsReference, sample.id=sampleId, snp.id=snpId,
                                snpfirstdim=TRUE)
        snpgdsClose(gdsReference)
    }
    #studyAn <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    gS <- snpgdsGetGeno(gdsProfile, sample.id=listProfile, snp.id=snpId,
                                snpfirstdim=TRUE)
    if( profileOnly){
        g<- gS
    }else{
        g <- cbind(g,gS)
    }

    write_bed(fileOut, g, verbose=TRUE, append=FALSE)
    snpgdsClose(gdsProfile)

    return(0)
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
#' @param listProfile a \code{vector} of \code{character} representing the list
#' of profile from gdsProfile to keep.
#'
#' @param listRM a TODO. Default: \code{NULL}.
#'
#' @param profileOnly a \code{integer} TODO. Default: \code{FALSE}.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters.
#'
#' @return The integer \code{0L} when successful.
#'
#' @examples
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' pathOut <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathOut)){
#'     dir.create(pathOut)
#' }
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#'
#' writeFamProfile(fileOut=file.path(pathOut,
#'     paste0(pedStudy(pRAIDS)$Name.ID[1], ".fam")),
#'     listProfile=pedStudy(pRAIDS)$Name.ID[1],
#'     profileOnly = TRUE, pRAIDS=pRAIDS)
#'
#' dir(file.path(pathOut))
#'
#' unlink(file.path(pathOut, paste0(pedStudy(pRAIDS)$Name.ID[1], ".fam")))
#' unlink(file.path(pathOut), recursive=TRUE)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom genio write_fam
#' @importFrom methods is
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
writeFamProfile <- function(fileOut, listProfile, listRM=NULL,
    profileOnly=FALSE, pRAIDS){

    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                            paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    sampleId <- NULL
    if(! profileOnly){
        gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))

        sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
        if(! is.null(listRM)){
            sampleId <- sampleId[! sampleId %in% listRM ]
        }
        snpgdsClose(gdsReference)
    }

    sujetFam <- data.frame(fam=seq_len(length(sampleId) + length(listProfile)),
                           id=c(sampleId, listProfile),
                           pat=0,
                           mat=0,
                           sex=1,
                           pheno=0,
                           stringsAsFactors=FALSE)

    write_fam(fileOut, sujetFam)
    snpgdsClose(gdsProfile)

    return(0)
}

#' @title Generate the row for the ADMIXTURE P-matrix
#'
#' @description This function write with P-matrix
#'
#' @param anchor a \code{vector} of \code{character} representing the
#' ancestral group for each profile from the reference
#' where sample.ref == 1. If the admixture proportion is higher than
#' the threshold for one ancestral group, the value is the name
#' of the ancestral group; else "-".
#'
#' @param snpId a \code{vector} of \code{character} representing the list
#' of snp.id from gdsReference to keep.
#'
#' @param matGr a \code{matrix} containing the sample identifiers and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters.
#'
#' @return a \code{list} containing one entry per ancestral group:
#' \describe{
#' \item{\code{ancestral group}}{ a \code{list} contenains 2 entries:
#' \describe{
#' \item{\code{f}}{ a \code{list} contenains number of rows in
#' matGr + 1 entries:
#' \describe{
#' \item{\code{vector}}{ the first is frequency for all the anchor of
#' each snpId}
#' \item{\code{vector}}{ the other is frequency for all the anchor
#' minus one row of matGr of each snpId}
#' }
#' }
#' \item{\code{index}}{ a \code{vector} contenains 2 entriescontenains
#' number of rows in matGr entries. Each value is the index of the row
#' for which the f of this index corresponds.}
#' }
#' }
#' }
#'
#' @examples
#'
#' ## Required library
#' library(SNPRelate)
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop")
#'
#'
#' matGr <- getMatrixPopSynthetic(pRAIDS=pRAIDS)
#' fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
#'                             paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
#'
#'
#' gdsProfile <- snpgdsOpen(fileGDSProfile)
#' pruned <- read.gdsn(index.gdsn(gdsProfile, "pruned.study"))
#' snpgdsClose(gdsProfile)
#' gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))
#'
#' admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
#' pos <- admixDesc[admixDesc$id == anchorId(pRAIDS), "id.seq"]
#'
#' anchor <- read.gdsn(index.gdsn(gdsReference, "anchor.ref"),
#'                     start=c(1, pos), count=c(-1,1))
#' snpgdsClose(gdsReference)
#'
#' rowF <- prepPMatrix( anchor=anchor, snpId = pruned, matGr=matGr, pRAIDS=pRAIDS)
#' names(rowF)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom genio write_bed write_bim
#' @importFrom methods is
#' @importFrom SNPRelate snpgdsGetGeno snpgdsOpen snpgdsClose
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
prepPMatrix <- function( anchor,
        snpId, matGr=NULL, pRAIDS) {
    # c("EURag", "EASag", "AMRag", "AFRag", "SASag")
    tmpTime <- system.time({ag <- unique(anchor)
    ag <- ag[ag != "-"] # Ancestry group with an anchor
    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))
    profileRef <- read.gdsn(index.gdsn(gdsReference, "sample.id"))
    profileRef <- profileRef[which(read.gdsn(index.gdsn(gdsReference, "sample.ref"))==1)]
    rowF <- lapply( ag,
        FUN=function(x,gdsReference,profileRef,
                    pruned, anchor, matGr){
            keepProfiles <- profileRef[which(anchor == x)]
            matGeno <- snpgdsGetGeno(gdsReference,
                    sample.id=keepProfiles,
                    snp.id=pruned, snpfirstdim=FALSE)
            f <- list()
            # f[[1]] frequencies of SNPs based on all pseudo no-admix
            f[[1]]  <- colSums(matGeno)/ (2*dim(matGeno)[1])
            k<-1

            matchinF <- NULL
            if(! is.null(matGr)){
                # vector of integer with  length of nb groups
                # it contains the position in f
                # for the specific group
                matchinF <- integer(nrow(matGr))
                for(i in seq_len(nrow(matGr))){ # each group
                    # By default use f[[1]]
                    matchinF[i] <- 1
                    tmp <- which(keepProfiles %in% matGr[i,])
                    if(length(tmp) > 0){
                        # remove some anchors for the freq
                        k <- k + 1
                        f[[k]]  <- colSums(matGeno[-1 * tmp, ]) /
                            (2*dim(matGeno[ -1 * tmp, ])[1])
                        # The group i use f[[k]]
                        matchinF[i] <- k
                    }
                }
            }
            # f <- colSums(matGeno)/ (2*dim(matGeno)[1])
            return(list(f = f, index = matchinF))
            },
        gdsReference=gdsReference,
        profileRef=profileRef,
        pruned=snpId,
        anchor = anchor,
        matGr = matGr)})
    if(verbose(pRAIDS)) print(paste0("Time to prepare P-matrices ",
        " user ", unname(tmpTime[1]),
        " system ", unname(tmpTime[2]),
        " elapsed ", unname(tmpTime[3])
    ))
    names(rowF) <- ag
    snpgdsClose(gdsReference)
    return(rowF)
}

#' @title Generate the row for the ADMIXTURE P-matrax
#'
#' @description This function write with P-matrix
#'
#' @param pathOut a \code{character} representing the output
#' path
#'
#' @param rowF a \code{list} containing allele frequencies for each
#' population and and for a subset if intersect quasi-no admix subset
#' intersect the synthetic
#'
#' @param pRAIDS a \code{RAIDSparam} an object with all the RAIDS
#' parameters.
#'
#' @return a \code{array} of the frequency matrix.
#'
#' @examples
#'
#' ## see computePmatrix
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn index.gdsn
#' @importFrom genio write_bed write_bim
#' @importFrom methods is
#' @importFrom SNPRelate snpgdsGetGeno snpgdsOpen snpgdsClose
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
writePMatrix <- function( pathOut,rowF, pRAIDS) {
    # c("EURag", "EASag", "AMRag", "AFRag", "SASag")
    #  nb is the number of group with a matrix different of rowF[[x]]$f[[1]]
    nb <- length(rowF[[1]]$index)

    k <- 0
    matPIndex <- integer(nb)
    matP <- lapply(seq_len(length(rowF)),
        FUN=function(x, rowF){
            return(rowF[[x]]$f[[1]])
        },
        rowF=rowF
        )
    matP <- t(do.call(rbind, matP))
    write.table(round(1-matP, 5), file.path(pathOut, paste0( pedStudy(pRAIDS)$Name.ID[1], ".P.in")),row.names = FALSE,col.names = FALSE, sep=" ")

    for(i in seq_len(nb)){
        flag <- FALSE
        # vector of length number ancestry groups
        matI <- integer(length(rowF))
        for(j in seq_len(length(rowF))){
            matI[j] <- rowF[[j]]$index[i]
            if(rowF[[j]]$index[i] != 1){
                flag <- TRUE
            }
        }
        if(flag){
            matP <- lapply(seq_len(length(rowF)),
                FUN=function(x, rowF, matI){
                    return(rowF[[x]]$f[[matI[x]]])
                },
                rowF=rowF,
                matI = matI
            )
            matP <- t(do.call(rbind, matP))
            k<-k+1
            write.table(round(1-matP, 5), file.path(pathOut, paste0( pedStudy(pRAIDS)$Name.ID[1],".syn.", k, ".P.in")),row.names = FALSE,col.names = FALSE, sep=" ")
        }
        matPIndex[i] <- ifelse(! flag,0, k)

    }
    return(matPIndex)
}

#' @title Generate a matrix of synthetic group per pop
#'
#' @description This function extract the SNPs that pass a frequency cut-off
#' in at least one super population
#' from a GDS SNP information file and save the retained SNP information into
#' a VCF file.
#'
#' @param pRAIDS a \code{RAIDSparam} object with all the RAIDS
#' parameters
#'
#' @return a \code{matrix} containing the sample identifiers and where
#' each column is the name of a subcontinental population. The number of
#' row corresponds to the number of samples for each subcontinental population.
#'
#' @examples
#'
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                       study.desc=studyID,
#'                       study.platform="NotDef",
#'                       stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                        Case.ID=c("ex1"),
#'                        Sample.Type=c("type"),
#'                        Diagnosis="NotDef",
#'                        Source=c("NotDef"),
#'                        stringsAsFactors=FALSE)
#' row.names(pedStudy) <- pedStudy$Name.ID
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                      pedStudy=pedStudy,
#'                      pathProfileGDS=pathProfileGDS,
#'                      fileReferenceGDS=fileReferenceGDS)
#' matGr <- getMatrixPopSynthetic(pRAIDS)
#' dim(matGr)
#'
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt read.gdsn
#' @importFrom genio write_bed
#' @importFrom methods is
#' @importFrom S4Vectors isSingleNumber
#' @importFrom utils write.table
#' @encoding UTF-8
#' @export
getMatrixPopSynthetic <- function(pRAIDS) {
    fileGDSProfile <- file.path(pathProfileGDS(pRAIDS),
                                paste0(pedStudy(pRAIDS)$Name.ID[1], ".gds"))
    gdsProfile <- snpgdsOpen(fileGDSProfile)
    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))
    tmp <- which(read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1)
    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[tmp]

    sampleAnnot <- read.gdsn(index.gdsn(gdsReference, "sample.annot"))[tmp,]
    # studySyn <- read.gdsn(index.gdsn(gdsProfile,""))
    # pRAIDS$studyDFSyn$study.id[1]
    profileSyn <- read.gdsn(index.gdsn(gdsProfile, "study.annot"))
    profileSyn <- profileSyn[profileSyn$study.id ==
                                        studyDFSyn(pRAIDS)$study.id[1],]

    tmp <- which(sampleId %in% profileSyn$case.id)
    listPop <- unique(sampleAnnot[tmp, fieldSubPop(pRAIDS)])
    listGr <- lapply(listPop,
                    FUN=function(x,sampleSyn, sampleAnnotSyn, fieldSubPop){
                        return(sampleSyn[sampleAnnotSyn[,fieldSubPop] == x])
                },
                sampleSyn=sampleId[tmp],
                sampleAnnotSyn=sampleAnnot[tmp,],
                fieldSubPop=fieldSubPop(pRAIDS))
    matGr <- do.call(cbind, listGr)
    snpgdsClose(gdsProfile)
    snpgdsClose(gdsReference)

    return(matGr)
}


#' @title Call ADMIXTURE
#'
#' @description TODO
#'
#'
#' @param pathOut  a \code{character} representing path
#' where the file for ADMIXTURE (BED,BIM, FAM, and P.in)
#' are located and where the output of ADMIXTURE
#' will be generated.
#'
#' @param fileBed  a \code{character} represented
#' the BED file name (without path) for a
#'
#'
#' @param seedAdmix  a \code{integer} representing the seed
#' pass to ADMIXTURE.
#' Default: \code{3435L}.
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
#' \dontrun{
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' pathOut <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathOut)){
#'     dir.create(pathOut)
#' }
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#'
#' ## Path to admixture
#' ## pathADMIXTURE <- "PATH/admixture"
#'
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop",
#'                       pathADMIXTURE=pathADMIXTURE)
#' matPIndex <- RAIDS:::computePmatrix(pathOut=pathOut, pRAIDS=pRAIDS)
#' RAIDS:::generateProfilePlinkFiles(pathOut, pRAIDS)
#'
#' profileQ <- runADMIXTURE(pathOut = file.path(pathOut, "profileAdmixture"),
#'     fileBed = paste0(pedStudy(pRAIDS)$Name.ID[1], ".bed"),
#'     seedAdmix=seedADMIXTURE(pRAIDS), pRAIDS)
#'
#'
#' dir(file.path(pathOut,"profileAdmixture"))
#' print(profileQ)
#' unlink(file.path(pathOut), recursive = TRUE)
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen
#' @encoding UTF-8
#' @export
runADMIXTURE <- function( pathOut, fileBed, seedAdmix=3435L, pRAIDS) {

    if(is.null(pathADMIXTURE(pRAIDS)) || ! file.exists(pathADMIXTURE(pRAIDS))){
        stop("Fatal, you need to defined pathADMIXTURE in RAIDSparam to use ADMIXTURE")
    }
    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))

    admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
    pos <- admixDesc[admixDesc$id == anchorId(pRAIDS), "id.seq"]
    popSelected <- read.gdsn(index.gdsn(gdsReference, "admixture.ag"))
    popSelected <- popSelected[popSelected$id.seq == pos,]

    snpgdsClose(gdsReference)

    pathCur <- getwd()
    setwd(pathOut)
    on.exit(setwd(pathCur), add=TRUE)
    cmdAdmix <- c("-s", seedAdmix, "-P",
            fileBed,
            nrow(popSelected))

    # cmdAdmix <- c("-s", seedAdmix, "-P",
    #         pathBed,
    #          "--out", outPref,
    #     5)
    # TODO select the population group

    res <- system2(command = pathADMIXTURE(pRAIDS), args = cmdAdmix,
        stdout=TRUE)
    matQ <- read.delim(gsub(".bed$",
                    paste0(".", nrow(popSelected), ".Q"), fileBed),
                header = FALSE, sep=" ")
    colnames(matQ) <- popSelected[order(popSelected$column), "ag"]
    famF <- read.delim(gsub(".bed$",
                    paste0(".fam"), fileBed),
                header = FALSE, sep="\t")
    rownames(matQ) <- famF[,2]
    setwd(pathCur)

    return(matQ)
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
#' \dontrun{
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoAncestryCall", package="RAIDS")
#'
#' pathOut <- file.path(tempdir(), "out.tmp")
#' if(!dir.exists(pathOut)){
#'     dir.create(pathOut)
#' }
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA"
#'
#' studyDF <- data.frame(study.id=studyID,
#'                              study.desc=studyID,
#'                              study.platform="NotDef",
#'                              stringsAsFactors=FALSE)
#' pathProfileGDS <- file.path(dataDir) # , "ex1.gds"
#' fileReferenceGDS <- system.file("extdata/tests/ex1_good_small_1KG.gds", package="RAIDS")
#' pedStudy <- data.frame(Name.ID=c("ex1"),
#'                              Case.ID=c("ex1"),
#'                              Sample.Type=c("type"),
#'                              Diagnosis="NotDef",
#'                              Source=c("NotDef"),
#'                              stringsAsFactors=FALSE)
#'      row.names(pedStudy) <- pedStudy$Name.ID
#' ## Path to admixture
#' ## pathADMIXTURE <- "PATH/admixture"
#' pRAIDS <- RAIDSparam(studyDF=studyDF,
#'                       pedStudy=pedStudy,
#'                       pathProfileGDS=pathProfileGDS,
#'                       fileReferenceGDS=fileReferenceGDS,
#'                       fieldPopInfAnc="SuperPop",
#'                       pathADMIXTURE=pathADMIXTURE)
#' matPIndex <- RAIDS:::computePmatrix(pathOut=pathOut, pRAIDS=pRAIDS)
#' RAIDS:::generateProfilePlinkFiles(pathOut, pRAIDS)
#'
#' syntheticQ <- runSynADMIXTURE(pathOut=pathOut, matPIndex=matPIndex, pRAIDS=pRAIDS)
#'
#' dir(file.path(pathOut))
#' head(syntheticQ)
#' unlink(file.path(pathOut), recursive = TRUE)
#' }
#'
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @importFrom gdsfmt index.gdsn read.gdsn closefn.gds
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom BiocParallel MulticoreParam bplapply
#' @encoding UTF-8
#' @export
runSynADMIXTURE <- function( pathOut, matPIndex, pRAIDS) {

    gdsReference <- snpgdsOpen(fileReferenceGDS(pRAIDS))

    admixDesc <- read.gdsn(index.gdsn(gdsReference, "admixture.desc"))
    pos <- admixDesc[admixDesc$id == anchorId(pRAIDS), "id.seq"]
    popSelected <- read.gdsn(index.gdsn(gdsReference, "admixture.ag"))
    popSelected <- popSelected[popSelected$id.seq == pos,]
    matQGS <- read.gdsn(index.gdsn(gdsReference, paste0("admixture.", nrow(popSelected))), start=c(1,1,pos), count=c(-1,-1,1))
    sampleId <- read.gdsn(index.gdsn(gdsReference, "sample.id"))[read.gdsn(index.gdsn(gdsReference, "sample.ref")) == 1]
    row.names(matQGS) <- sampleId
    colnames(matQGS) <- paste0(popSelected[order(popSelected$column), "ag"], ".gs")
    snpgdsClose(gdsReference)
    matGr <- getMatrixPopSynthetic(pRAIDS)
    res <- NULL
    if(np(pRAIDS) == 1){
        res <- lapply(seq_len(nrow(matGr)),
            FUN=function(x, pathOut, matGr, matPIndex, matQGS, pRAIDS){
                matQ <- runGrSynADMIXTURE(pathOut=pathOut, gr=x,
                    matGr=matGr, matPIndex=matPIndex,
                    matQGS=matQGS, pRAIDS=pRAIDS)

                return(matQ)
            },
            pathOut=pathOut,
            matGr=matGr,
            matPIndex=matPIndex,
            matQGS=matQGS,
            pRAIDS=pRAIDS)
    }else{
        param <- MulticoreParam(workers = np(pRAIDS))
        res <- bplapply(seq_len(nrow(matGr)),
            FUN=function(x, pathOut, matGr, matPIndex, matQGS, pRAIDS){
                matQ <- runGrSynADMIXTURE(pathOut=pathOut, gr=x,
                    matGr=matGr, matPIndex=matPIndex,
                    matQGS=matQGS, pRAIDS=pRAIDS)

                return(matQ)
            },
        pathOut=pathOut,
        matGr=matGr,
        matPIndex=matPIndex,
        matQGS=matQGS,
        pRAIDS=pRAIDS,
        BPPARAM = param)
    }
    res <- do.call(rbind, res)
    return(res)
}

#' aicsPaper: Ancestry Inference from Cancer Sequencing
#'
#' The aicsPaper package implements specialized algorithms that enable
#' ancestry inference from various cancer sequencing sources (RNA,
#' Exome and Whole-Genome sequencing).
#'
#' The aicsPaper package also implements simulation algorithm that generates
#' synthetic cancer-derived data.
#'
#' This code and analysis pipeline was designed and developed for the
#' following publication:
#'
#' @references
#' Belleau, P, Deschênes A, Tuveson DA, Krasnitz A. Accurate and robust
#' inference of genetic ancestry from cancer-derived molecular data across
#' genomic platforms. bioRxiv 2022.02.01.478737
#' [link to bioRxiv](https://doi.org/10.1101/2022.02.01.478737)
#'
#'
#' @docType package
#'
#' @name aicsPaper-package
#'
#' @aliases aicsPaper-package aicsPaper
#'
#' @author Pascal Belleau, Astrid Deschênes and
#' Alexander Krasnitz
#'
#' Maintainer:
#' Pascal Belleau <pascal_belleau@hotmail.com>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{prepPed1KG}} {This function extracts the
#'     needed information from the 1000 Genomes pedigree file and formats it
#'     into a \code{data.frame} so in can be used in following steps
#'     of the ancestry inference process.}
#'     \item \code{\link{generateMapSnvSel}} {The function applies a cut-off
#'     filter to the SNP information file to retain only the SNP that have a
#'     frequency superior or equal to the specified cut-off in at least one
#'     super population.}
#'     \item \code{\link{generateGDS1KG}} {This function generates the GDS
#'     file that will contain the information from 1KG.}
#'     \item \code{\link{identifyRelative}} {The function identify patients
#'     that are genetically related in the GDS 1KG file.}
#'     \item \code{\link{addRef2GDS1KG}} {TODO}
#' }
#'
#' @encoding UTF-8
#' @keywords package
NULL

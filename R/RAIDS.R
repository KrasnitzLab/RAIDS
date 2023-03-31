#' RAIDS: Ancestry Inference from Cancer Sequencing
#'
#' The RAIDS package implements specialized algorithms that enable
#' ancestry inference from various cancer sequencing sources (RNA,
#' Exome and Whole-Genome sequencing).
#'
#' The RAIDS package also implements simulation algorithm that generates
#' synthetic cancer-derived data.
#'
#' This code and analysis pipeline was designed and developed for the
#' following publication:
#'
#' @references
#' Pascal Belleau, Astrid Deschênes, Nyasha Chambwe, David A. Tuveson,
#' Alexander Krasnitz; Genetic Ancestry Inference from Cancer-Derived
#' Molecular Data across Genomic and Transcriptomic Platforms.
#' Cancer Res 1 January 2023; 83 (1): 49–58.
#' https://doi.org/10.1158/0008-5472.CAN-22-0682
#' [link to paper](https://doi.org/10.1158/0008-5472.CAN-22-0682)
#'
#'
#' @docType package
#'
#' @name RAIDS-package
#'
#' @aliases RAIDS-package RAIDS
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
#'     file that will contain the information from 1KG. }
#'     \item \code{\link{identifyRelative}} {The function identify patients
#'     that are genetically related in the 1KG GDS file. }
#'     \item \code{\link{addRef2GDS1KG}} { This function adds the information
#'     about the unrelated patients to the 1KG GDS file. }
#'     \item \code{\link{add1KG2SampleGDS}} { TODO }
#'     \item \code{\link{addBlockFromPlink2GDS}} { TODO }
#'     \item \code{\link{addPhase1KG2SampleGDSFromGDS}} { TODO }
#'     \item \code{\link{appendStudy2GDS1KG}} { This function creates the
#'     Sample GDS file(s) for one or multiple specific samples
#'     using the information from a Sample RDS description file and the 1KG
#'     GDS file. }
#'     \item \code{\link{basePCASample}} { TODO }
#'     \item \code{\link{estimateAllelicFraction}} { This function estimates
#'     the allelic fraction of the pruned SNVs for a specific sample and add
#'     the information to the associated GDS Sample file. The allelic fraction
#'     estimation method is adapted to the type of study (DNA or RNA). }
#' }
#'
#' @encoding UTF-8
#' @keywords package
NULL
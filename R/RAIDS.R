#' RAIDS: Accurate Inference of Genetic Ancestry from Cancer Sequences
#'
#' The RAIDS package implements specialized algorithms that enable
#' ancestry inference from various cancer data sources (RNA,
#' Exome and Whole-Genome sequencing).
#'
#' The RAIDS package also implements simulation algorithm that generates
#' synthetic cancer-derived data.
#'
#' This code and analysis pipeline was designed and developed for the
#' following publication:
#'
#' Pascal Belleau, Astrid Deschênes, Nyasha Chambwe, David A. Tuveson,
#' Alexander Krasnitz; Genetic Ancestry Inference from Cancer-Derived
#' Molecular Data across Genomic and Transcriptomic Platforms.
#' Cancer Res 1 January 2023; 83 (1): 49–58.
#' https://doi.org/10.1158/0008-5472.CAN-22-0682
#'
#'
#' @references
#' Pascal Belleau, Astrid Deschênes, Nyasha Chambwe, David A. Tuveson,
#' Alexander Krasnitz; Genetic Ancestry Inference from Cancer-Derived
#' Molecular Data across Genomic and Transcriptomic Platforms.
#' Cancer Res 1 January 2023; 83 (1): 49–58.
#' https://doi.org/10.1158/0008-5472.CAN-22-0682
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
#'     \item \code{\link{runExomeAncestry}} {This function runs most steps
#'     leading to the ancestry inference call on a specific exome profile.}
#'     \item \code{\link{runExomeAncestry}} {This function runs most steps
#'     leading to the ancestry inference call on a specific RNA profile.}
#'     \item \code{\link{estimateAllelicFraction}} { This function estimates
#'     the allelic fraction of the pruned SNVs for a specific sample and add
#'     the information to the associated GDS Sample file. The allelic fraction
#'     estimation method is adapted to the type of study (DNA or RNA). }
#'     \item \code{\link{computeSyntheticROC}} { This function calculate the
#'     AUROC of the inferences for specific values of D and K using the
#'     inferred ancestry results from the synthetic profiles.}
#'     \item \code{\link{generateMapSnvSel}} {The function applies a cut-off
#'     filter to the SNV information file to retain only the SNV that have a
#'     frequency superior or equal to the specified cut-off in at least one
#'     super population.}
#' }
#'
#' @return RAIDS
#' @encoding UTF-8
#' @keywords package
NULL

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


#' A small \code{data.frame} containing the information related to
#' synthetic profiles. The ancestry of the profiles used to generate the
#' synthetic profiles must be present.
#'
#' The object is a \code{data.frame} with 7 columns. The row names of
#' the \code{data.frame} must be the profile unique identifiers.
#'
#' This dataset can be
#' used to test the \code{\link{computeSyntheticROC}} function.
#'
#' @name pedSynthetic
#'
#' @docType data
#'
#' @aliases pedSynthetic
#'
#' @format The \code{data.frame} containing the information about the
#' synthetic profiles.  The row names of
#' the \code{data.frame} correspond to the profile unique identifiers.
#' The \code{data.frame} contains 7 columns:
#' \itemize{
#' \item \code{data.id} {a \code{character} string representing the unique
#' synthetic profile identifier.}
#' \item \code{case.id} {a \code{character} string representing the unique
#' profile identifier that was used to generate the synthetic profile.}
#' \item \code{sample.type} {a \code{character} string representing the type
#' of profile. }
#' \item \code{diagnosis} {a \code{character} string representing the
#' diagnosis of profile that was used to generate the synthetic profile. }
#' \item \code{source} {a \code{character} string representing the
#' source of the synthetic profile. }
#' \item \code{study.id} {a \code{character} string representing the
#' name of the study to which the synthetic profile is associated. }
#' \item \code{superPop} {a \code{character} string representing the
#' super population of the profile that was used to generate the synthetic
#' profile. }
#' }
#'
#' @return  The \code{data.frame} containing the information about the
#' synthetic profiles.  The row names of
#' the \code{data.frame} correspond to the profile unique identifiers.
#' The \code{data.frame} contains 7 columns:
#' \itemize{
#' \item \code{data.id} {a \code{character} string representing the unique
#' synthetic profile identifier.}
#' \item \code{case.id} {a \code{character} string representing the unique
#' profile identifier that was used to generate the synthetic profile.}
#' \item \code{sample.type} {a \code{character} string representing the type
#' of profile. }
#' \item \code{diagnosis} {a \code{character} string representing the
#' diagnosis of profile that was used to generate the synthetic profile. }
#' \item \code{source} {a \code{character} string representing the
#' source of the synthetic profile. }
#' \item \code{study.id} {a \code{character} string representing the
#' name of the study to which the synthetic profile is associated. }
#' \item \code{superPop} {a \code{character} string representing the
#' super population of the profile that was used to generate the synthetic
#' profile. }
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{computeSyntheticROC}} {for calculating the AUROC of
#'     the inferences for specific values of D and K using the inferred
#'     ancestry results from the synthetic profiles}
#' }
#'
#' @usage data(pedSynthetic)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading demo dataset containing pedigree information for synthetic
#' ## profiles
#' data(pedSynthetic)
#'
#' ## Loading demo dataset containing the inferred ancestry results
#' ## for the synthetic data
#' data(matKNNSynthetic)
#'
#' ## Retain one K and one D value
#' matKNN <- matKNNSynthetic[matKNNSynthetic$D == 5 & matKNNSynthetic$K == 4, ]
#'
#' ## Compile statistics from the
#' ## synthetic profiles for fixed values of D and K
#' results <- RAIDS:::computeSyntheticROC(matKNN=matKNN,
#'     matKNNAncestryColumn="SuperPop",
#'     pedCall=pedSynthetic, pedCallAncestryColumn="superPop",
#'     listCall=c("EAS", "EUR", "AFR", "AMR", "SAS"))
#'
#' results$matAUROC.All
#' results$matAUROC.Call
#' results$listROC.Call
#'
#'
NULL


#' A small \code{data.frame} containing the
#' inferred ancestry on the synthetic profiles.
#'
#' The object is a \code{data.frame} with 4 columns.
#'
#' This dataset can be
#' used to test the \code{\link{computeSyntheticROC}} function.
#'
#' @name matKNNSynthetic
#'
#' @docType data
#'
#' @aliases matKNNSynthetic
#'
#' @format The \code{data.frame} containing the information about the
#' synthetic profiles. The \code{data.frame} contains 4 columns:
#' \itemize{
#' \item \code{sample.id} {a \code{character} string representing the unique
#' synthetic profile identifier.}
#' \item \code{D} {a \code{numeric} representing the number of dimensions used
#' to infer the ancestry of the synthetic profile.}
#' \item \code{K} {a \code{numeric} representing the number of neighbors used
#' to infer the ancestry of the synthetic profile.}
#' \item \code{SuperPop} {a \code{character} string representing the
#' inferred ancestry of the synthetic profile for the specific D and K values.}
#' }
#'
#' @return  The \code{data.frame} containing the information about the
#' synthetic profiles. The \code{data.frame} contains 4 columns:
#' \itemize{
#' \item \code{sample.id} {a \code{character} string representing the unique
#' synthetic profile identifier.}
#' \item \code{D} {a \code{numeric} representing the number of dimensions used
#' to infer the ancestry of the synthetic profile.}
#' \item \code{K} {a \code{numeric} representing the number of neighbors used
#' to infer the ancestry of the synthetic profile.}
#' \item \code{SuperPop} {a \code{character} string representing the
#' inferred ancestry of the synthetic profile for the specific D and K values.}
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{computeSyntheticROC}} {for calculating the AUROC of
#'     the inferences for specific values of D and K using the inferred
#'     ancestry results from the synthetic profiles}
#' }
#'
#' @usage data(matKNNSynthetic)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading demo dataset containing pedigree information for synthetic
#' ## profiles
#' data(pedSynthetic)
#'
#' ## Loading demo dataset containing the inferred ancestry results
#' ## for the synthetic data
#' data(matKNNSynthetic)
#'
#' ## Retain one K and one D value
#' matKNN <- matKNNSynthetic[matKNNSynthetic$D == 5 & matKNNSynthetic$K == 4, ]
#'
#' ## Compile statistics from the
#' ## synthetic profiles for fixed values of D and K
#' results <- RAIDS:::computeSyntheticROC(matKNN=matKNN,
#'     matKNNAncestryColumn="SuperPop",
#'     pedCall=pedSynthetic, pedCallAncestryColumn="superPop",
#'     listCall=c("EAS", "EUR", "AFR", "AMR", "SAS"))
#'
#' results$matAUROC.All
#' results$matAUROC.Call
#' results$listROC.Call
#'
#'
NULL


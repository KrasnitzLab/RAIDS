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


#' A small \code{data.frame} containing the
#' SNV information.
#'
#' The object is a \code{data.frame} with 17 columns.
#'
#' This dataset can be
#' used to test the \code{\link{calcAFMLRNA}} and \code{\link{tableBlockAF}}
#' internal functions.
#'
#' @name snpPositionDemo
#'
#' @docType data
#'
#' @aliases snpPositionDemo
#'
#' @format The \code{data.frame} containing the information about the
#' synthetic profiles. The \code{data.frame} contains 4 columns:
#' \itemize{
#' \item \code{cnt.tot} {a \code{integer} representing the number of reads at
#' the SNV position.}
#' \item \code{cnt.ref} {a \code{integer} representing the number of reads
#' corresponding to the reference at the SNV position.}
#' \item \code{cnt.alt} {a \code{integer} representing the number of reads
#' different than the reference at the SNV position.}
#' \item \code{snp.pos} {a \code{integer} representing the position of the
#' SNV on the chromosome.}
#' \item \code{snp.chr} {a \code{integer} representing the chromosome on which
#' the SNV is located.}
#' \item \code{normal.geno} {a \code{integer} representing the genotype
#' (0=wild-type reference; 1=heterozygote; 2=homozygote alternative; 3=unkown).}
#' \item \code{pruned} {a \code{logical} indicated if the SNV is pruned.}
#' \item \code{snp.index} {a \code{integer} representing the index of the
#' SNV in the reference SNV GDS file.}
#' \item \code{keep} {a \code{logical} indicated if the genotype
#' exists for the SNV.}
#' \item \code{hetero} {a \code{logical} indicated if the SNV is heterozygote.}
#' \item \code{homo} {a \code{logical} indicated if the SNV is homozygote.}
#' \item \code{block.id} {a \code{integer} representing the block identifier
#' associated to the current SNV.}
#' \item \code{phase} {a \code{integer} representing the block identifier
#' associated to the current SNV.}
#' \item \code{lap} {a \code{numeric} representing the lower allelic fraction.}
#' \item \code{LOH} {a \code{integer} indicating if the SNV is in an LOH region
#' (0=not LOH, 1=in LOH).}
#' \item \code{imbAR} {a \code{integer} indicating if the SNV is in an
#' imbalanced region (-1=not classified as imbalanced or LOH, 0=in LOH;
#' 1=tested positive for imbalance in at least 1 window).}
#' \item \code{freq} {a \code{numeric} representing the frequency of the
#' variant in the the reference.}
#' }
#'
#' @return The \code{data.frame} containing the information about the
#' synthetic profiles. The \code{data.frame} contains 4 columns:
#' \itemize{
#' \item \code{cnt.tot} {a \code{integer} representing the number of reads at
#' the SNV position.}
#' \item \code{cnt.ref} {a \code{integer} representing the number of reads
#' corresponding to the reference at the SNV position.}
#' \item \code{cnt.alt} {a \code{integer} representing the number of reads
#' different than the reference at the SNV position.}
#' \item \code{snp.pos} {a \code{integer} representing the position of the
#' SNV on the chromosome.}
#' \item \code{snp.chr} {a \code{integer} representing the chromosome on which
#' the SNV is located.}
#' \item \code{normal.geno} {a \code{integer} representing the genotype
#' (0=wild-type reference; 1=heterozygote; 2=homozygote alternative; 3=unkown).}
#' \item \code{pruned} {a \code{logical} indicated if the SNV is pruned.}
#' \item \code{snp.index} {a \code{integer} representing the index of the
#' SNV in the reference SNV GDS file.}
#' \item \code{keep} {a \code{logical} indicated if the genotype
#' exists for the SNV.}
#' \item \code{hetero} {a \code{logical} indicated if the SNV is heterozygote.}
#' \item \code{homo} {a \code{logical} indicated if the SNV is homozygote.}
#' \item \code{block.id} {a \code{integer} representing the block identifier
#' associated to the current SNV.}
#' \item \code{phase} {a \code{integer} representing the block identifier
#' associated to the current SNV.}
#' \item \code{lap} {a \code{numeric} representing the lower allelic fraction.}
#' \item \code{LOH} {a \code{integer} indicating if the SNV is in an LOH region
#' (0=not LOH, 1=in LOH).}
#' \item \code{imbAR} {a \code{integer} indicating if the SNV is in an
#' imbalanced region (-1=not classified as imbalanced or LOH, 0=in LOH;
#' 1=tested positive for imbalance in at least 1 window).}
#' \item \code{freq} {a \code{numeric} representing the frequency of the
#' variant in the the reference.}
#' }
#'
#' @usage data(snpPositionDemo)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Loading demo dataset containing SNV information
#' data(snpPositionDemo)
#'
#' ## Only use a subset of heterozygote SNVs related to one block
#' subset <- snpPositionDemo[which(snpPositionDemo$block.id == 2750 &
#'                   snpPositionDemo$hetero), c("cnt.ref", "cnt.alt", "phase")]
#'
#' ## Compute the log likelihood ratio based on the coverage of
#' ## each allele in a specific block
#' result <- RAIDS:::calcAFMLRNA(subset)
#' head(result)
#'
#'
NULL


#' The PCA results of the demo 1KG reference dataset for demonstration purpose.
#' Beware that the PCA has been run on a very small subset of the
#' 1KG reference dataset
#' and should not be used to call ancestry inference on a real profile.
#'
#' The object is a \code{list}.
#'
#' This object can be
#' used to test the \code{\link{computePCAMultiSynthetic}} function.
#'
#' @name demoPCA1KG
#'
#' @docType data
#'
#' @aliases demoPCA1KG
#'
#' @format The \code{list} containing the PCA results for a small subset of
#' the reference 1KG dataset. The \code{list} contains 2 entries:
#' \itemize{
#' \item{pruned}{ a \code{vector} of SNV identifiers specifying selected SNVs
#' for the PCA analysis.}
#' \item{pca.unrel}{ a \code{snpgdsPCAClass} object containing the eigenvalues
#' as generated by \link[SNPRelate]{snpgdsPCA} function.}
#' }
#'
#' @return The \code{list} containing the PCA results for a small subset of
#' the reference 1KG dataset. The \code{list} contains 2 entries:
#' \itemize{
#' \item{pruned}{ a \code{vector} of SNV identifiers specifying selected SNVs
#' for the PCA analysis.}
#' \item{pca.unrel}{ a \code{snpgdsPCAClass} object containing the eigenvalues
#' as generated by \link[SNPRelate]{snpgdsPCA} function.}
#' }
#'
#' @usage data(demoPCA1KG)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Loading demo PCA on subset of 1KG reference dataset
#' data(demoPCA1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' samplesRM <- c("HG00246", "HG00325", "HG00611", "HG01173", "HG02165",
#'     "HG01112", "HG01615", "HG01968", "HG02658", "HG01850", "HG02013",
#'     "HG02465", "HG02974", "HG03814", "HG03445", "HG03689", "HG03789",
#'     "NA12751", "NA19107", "NA18548", "NA19075", "NA19475", "NA19712",
#'     "NA19731", "NA20528", "NA20908")
#' names(samplesRM) <- c("GBR", "FIN", "CHS","PUR", "CDX", "CLM", "IBS",
#'     "PEL", "PJL", "KHV", "ACB", "GWD", "ESN", "BEB", "MSL", "STU", "ITU",
#'     "CEU", "YRI", "CHB", "JPT", "LWK", "ASW", "MXL", "TSI", "GIH")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#'
#' ## Projects synthetic profiles on demo 1KG PCA
#' results <- computePCAMultiSynthetic(gdsProfile=gdsProfile,
#'     listPCA=demoPCA1KG, sampleRef=samplesRM, studyIDSyn=studyID,
#'     verbose=FALSE)
#'
#' ## The eigenvectors for the synthetic profiles
#' head(results$eigenvector)
#'
#' ## Close Profile GDS file (important)
#' closefn.gds(gdsProfile)
#'
NULL


#' The PCA result of demo synthetic profiles projected on the demo subset
#' 1KG reference PCA.
#'
#' The object is a \code{list}.
#'
#' This object can be
#' used to test the \code{\link{computeKNNRefSynthetic}} function.
#'
#' @name demoPCASyntheticProfiles
#'
#' @docType data
#'
#' @aliases demoPCASyntheticProfiles
#'
#' @format The \code{list} containing the PCA result of demo synthetic
#' profiles projected on the demo subset 1KG reference PCA.
#' The \code{list} contains 3 entries:
#' \itemize{
#' \item{sample.id}{ a \code{character} string representing the unique
#' identifier of the synthetic profiles.}
#' \item{eigenvector.ref} { a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{eigenvector}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current synthetic profiles projected on the demo
#' PCA 1KG reference profiles.}
#' }
#'
#' @return The \code{list} containing the PCA result of demo synthetic
#' profiles projected on the demo subset 1KG reference PCA.
#' The \code{list} contains 3 entries:
#' \itemize{
#' \item{sample.id}{ a \code{character} string representing the unique
#' identifier of the synthetic profiles.}
#' \item{eigenvector.ref} { a \code{matrix} of \code{numeric} containing
#' the eigenvectors for the reference profiles.}
#' \item{eigenvector}{ a \code{matrix} of \code{numeric} containing the
#' eigenvectors for the current synthetic profiles projected on the demo
#' PCA 1KG reference profiles.}
#' }
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{computeKNNRefSynthetic}} {for running a k-nearest
#'     neighbors analysis on a subset of the synthetic data set.}
#' }
#'
#' @usage data(demoPCASyntheticProfiles)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' ## Projects synthetic profiles on 1KG PCA
#' results <- computeKNNRefSynthetic(gdsProfile=gdsProfile,
#'     listEigenvector=demoPCASyntheticProfiles,
#'     listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"), studyIDSyn=studyID,
#'     spRef=demoKnownSuperPop1KG)
#'
#' ## The inferred ancestry for the synthetic profiles for differents values
#' ## of D and K
#' head(results$matKNN)
#'
#' ## Close Profile GDS file (important)
#' closefn.gds(gdsProfile)
#'
NULL




#' The known super population ancestry of the demo 1KG reference profiles.
#'
#' The object is a \code{vector}.
#'
#' This object can be
#' used to test the \code{\link{computeKNNRefSynthetic}} and
#' \code{\link{computePoolSyntheticAncestryGr}} functions.
#'
#' @name demoKnownSuperPop1KG
#'
#' @docType data
#'
#' @aliases demoKnownSuperPop1KG
#'
#' @format The \code{vector} containing the know super population ancestry
#' for the demo 1KG reference profiles.
#'
#' @return The \code{vector} containing the know super population ancestry
#' for the demo 1KG reference profiles.
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{computeKNNRefSynthetic}} {for running a k-nearest
#'     neighbors analysis on a subset of the synthetic data set.}
#'     \item \code{\link{computePoolSyntheticAncestryGr}} { for running a
#'     PCA analysis using 1 synthetic profile from each sub-continental
#'     population.}
#' }
#'
#' @usage data(demoKnownSuperPop1KG)
#'
#' @keywords datasets
#'
#' @examples
#'
#' ## Required library
#' library(gdsfmt)
#'
#' ## Load the demo PCA on the synthetic profiles projected on the
#' ## demo 1KG reference PCA
#' data(demoPCASyntheticProfiles)
#'
#' ## Load the known ancestry for the demo 1KG reference profiles
#' data(demoKnownSuperPop1KG)
#'
#' ## Path to the demo Profile GDS file is located in this package
#' dataDir <- system.file("extdata/demoKNNSynthetic", package="RAIDS")
#'
#' ## Open the Profile GDS file
#' gdsProfile <- snpgdsOpen(file.path(dataDir, "ex1.gds"))
#'
#' # The name of the synthetic study
#' studyID <- "MYDATA.Synthetic"
#'
#' ## Projects synthetic profiles on 1KG PCA
#' results <- computeKNNRefSynthetic(gdsProfile=gdsProfile,
#'     listEigenvector=demoPCASyntheticProfiles,
#'     listCatPop=c("EAS", "EUR", "AFR", "AMR", "SAS"), studyIDSyn=studyID,
#'     spRef=demoKnownSuperPop1KG)
#'
#' ## The inferred ancestry for the synthetic profiles for differents values
#' ## of D and K
#' head(results$matKNN)
#'
#' ## Close Profile GDS file (important)
#' closefn.gds(gdsProfile)
#'
NULL

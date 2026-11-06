#' @title Class Union for data.frame or NULL
#' 
#' @description A virtual class that groups data.frame and \code{NULL} 
#' together.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @name DataFrameOrNULL-class
#' @rdname DataFrameOrNULL-class
#' @exportClass DataFrameOrNULL
setClassUnion("DataFrameOrNULL", members = c("data.frame", "NULL"))


#' @title Class Union for character or NULL
#' 
#' @description A virtual class that groups character and \code{NULL} 
#' together.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @name CharacterOrNULL-class
#' @rdname CharacterOrNULL-class
#' @exportClass CharacterOrNULL
setClassUnion("CharacterOrNULL", members = c("character", "NULL"))


#' @title An S4 class to represent the RAIDS parameters
#'
#' @slot studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factors). Default: 
#' \code{data.frame(study.id="NotDef", study.desc="NotDef", 
#' study.platform="NotDef", stringsAsFactors=FALSE)}.
#' 
#' @slot studyDFSyn a \code{data.frame} containing the information about the
#' synthetic data to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be in \code{character} strings (no factors). Default: 
#' \code{data.frame(study.id="NotDef.Synthetic", 
#' study.desc="NotDef synthetic data", study.platform="Synthetic", 
#' stringsAsFactors=FALSE)}.
#' 
#' @slot pedStudy a \code{data.frame} containing the information TODO 
#' with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", and "Source". All columns must be in
#' \code{character} strings (no factors). All row names should correspond to the
#' values in the "Name.ID" column.
#' Default: \code{data.frame(Name.ID=c("ProfileId"), Case.ID=c("ProfileId"), 
#' Sample.Type=c("type"), Diagnosis="NotDef", Source=c("NotDef"), 
#' stringsAsFactors=FALSE, row.names = c("ProfileId"))}.
#' 
#' @slot studyType a single \code{character} string representing the type of 
#' study. The possible choices are: "LD" and "GeneAware". The type of study 
#' affects the way the estimation of the allelic fraction is done. 
#' Default: \code{"LD"}.
#' 
#' @slot genoSource a single \code{character} string corresponding to the type 
#' of file with the genotype and the allele information of the profile that 
#' will be provided in the 'profileFile' slot. The valid options are: "VCF", 
#' "generic", "snp-pileup", and "bam". The "generic" format CSV file
#' must have at least these columns: 'Chromosome', 'Position', 'Ref', 'Alt', 
#' 'Count', 'File1R', and 'File1A'. The 'Count' is the depth at the 
#' specified position; 'FileR' is the depth of the reference allele, and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally, in the case of a "VCF" file, the file must have at least those 
#' genotype fields: GT, AD, and DP. Default: \code{NULL}.
#' 
#' @slot blockTypeId a single \code{character} string corresponding to 
#' the block type used to extract the block identifiers. The block type must 
#' be present in the GDS Reference Annotation file. 
#' Default: \code{"GeneS.Ensembl.Hsapiens.v86"}.
#' 
#' @slot reference a single \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' Default: \code{"1KGv1.0"}.
#' 
#' @slot genome a \code{character} string with one possible value:
#' 'HG38'. It specifies the genome used. Default: \code{"HG38"}.
#' 
#' @slot chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section. 
#' Default: \code{seqlengths(Hsapiens)[seq_len(25)]}.
#' 
#' @slot paramAncestry a \code{list} of parameters related to ancestry. 
#' The \code{list} should contain those three entries: \code{ScanBamParam}, 
#' \code{PileupParam}, and \code{yieldSize}.  
#' Default: \code{list(ScanBamParam=NULL, PileupParam=NULL, 
#' yieldSize=10000000)}.
#' 
#' @slot profileFile a \code{character} string representing the path to the
#' file with genotype and the allele information of the profile. The format 
#' of the accepted file is determined by the value in the 'genoSource' slot.
#' A profile would have a file with extension "*vcf.gz" when 'genoSource' slot 
#' is "VCF". If 'genoSource' is "generic" or "snp-pileup", then "*.txt.gz". 
#' If 'genoSource' is "bam", then "*.bam" (the file needs to be indexed with a
#' existing corresponding "*.bai" file). Default: \code{NULL}.
#' 
#' @slot profileFileGeno not used in the current version
#' 
#' @slot pathProfileGDS a single \code{character} string representing the 
#' path to the directory where the GDS Profile files will be created. If 
#' specified, the directory must exist. 
#' Default: \code{NULL}.
#' 
#' @slot fileReferenceGDS a single \code{character} string representing 
#' the file name of the Reference GDS file. If specified, the file  
#' must exist. Default: \code{NULL}.
#' 
#' @slot fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Population Reference GDS Annotation file. If specified,
#' the file must exist. Default: \code{NULL}.
#' 
#' @slot inferenceType a single \code{character} string representing the 
#' genotype ancestry inference method. The two possible values:
#' 'PCAknn' and 'haploAdmixture'. Default: \code{"PCAknn"}.
#' 
#' @slot sampleRef a \code{vector} of \code{character} strings that represent 
#' a subset of sample identifiers from the GDS Reference File. Those samples 
#' will be used for the ancestry inference. If \code{NULL}, 
#' all the samples are retained. This is usefull if you want to removed some 
#' specific samples (ex: samples with admixed ancestry).
#' Default: \code{NULL}.
#' 
#' @slot batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore. 
#' Default: \code{1L}.
#' 
#' @slot prefix a single \code{character} string representing at the
#' beginning of the synthetic name id. Default: \code{"1"}.
#' 
#' @slot nbSim a positive \code{integer} representing number of 
#' replicated synthitic data (restricted to 1 for 
#' the moment). Default: \code{1L}.
#' 
#' @slot posOffset a single \code{integer} that is added to the SNP position to
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' Default: \code{-1L}.
#' 
#' @slot minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#' 
#' @slot minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#' 
#' @slot seqError a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error. 
#' Default: \code{0.001}.
#' 
#' @slot seqErrorSyn a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error for synthetic.
#' Default: \code{0.001}.
#' 
#' @slot pRecomb a single positive \code{numeric} between \code{0} and 
#' \code{1} that represents the frequency of phase switching in the 
#' synthetic profiles. Default: \code{0.01}.
#' 
#' @slot np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#' 
#' @slot listPos a \code{data.frame} containing 2 columns named: 
#' "snp.chromosome" and "snp.position". The first column,
#' called "snp.chromosome", contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position", contains the
#' position of the SNV on the chromosome. Default: \code{NULL}.
#' 
#' @slot syntheticRefDF a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The \code{data.frame} must have those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#' Default: \code{NULL}.
#' 
#' @slot pruningMethod a single \code{character} string representing the 
#' method that will be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime", and "composite". Default: \code{"corr"}.
#' 
#' @slot slideWindowMaxBP a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done by the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. Default: \code{500000L}.
#' 
#' @slot thresholdLD a single positive \code{numeric} value that represents 
#' the LD threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. Default: \code{sqrt(0.1)}.
#' 
#' @slot specificSNV a \code{data.frame} containing 2 columns. The first 
#' column, called "snp.chromosome" contains the name of the chromosome where 
#' the SNV is located. The second column, called "snp.position", contains the
#' position of the SNV on the chromosome. Optionally, the column "snvKeep" 
#' contains the SNV index in the reference or -1 if not in the reference. It 
#' is used during the pruning step. Default: \code{NULL}.
#' 
#' @slot genoType a single \code{character} string representing the field 
#' in the GDS profile file that contains all the genotype information. 
#' Default: \code{"geno.ref"}.
#' 
#' @slot phaseType a TODO. Default: \code{"phase.ref"}.
#' 
#' @slot phase a \code{logical} indicating TODO. Default: \code{FALSE}.
#' 
#' @slot PCAmissingRate a positive \code{numeric} representing the maximum 
#' missing rate retained accepted to use SNPs in the PCA analysis done 
#' with the \link[SNPRelate]{snpgdsPCA} function. If \code{NaN}, no 
#' missing threshold. Default: \code{0.025}.
#' 
#' @slot PCAalgorithm a single \code{character} string representing the 
#' algorithm to use with the \link[SNPRelate]{snpgdsPCA} function. The 
#' algorithm must be implemented and available to the 
#' \link[SNPRelate]{snpgdsPCA} function. The current options are:
#' 'exact' or 'randomized'.
#' Default: \code{"exact"}.
#' 
#' @slot eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned. 
#' Default: \code{32L}.
#' 
#' @slot eigenCountSyn a TODO. Default: \code{15L}.
#' 
#' @slot kList a \code{vector} of positive \code{integer} representing the 
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. 
#' Default: \code{seq(2L, 15L, 1L)}.
#' 
#' @slot pcaList a \code{vector} of positive \code{integer} representing 
#' the values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis. 
#' Default: \code{seq(2L, 15L, 1L)}.
#' 
#' @slot fieldPopInRef a single \code{character} string representing the name 
#' of the column that contains the known ancestry for the reference profiles in
#' the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter). Default: \code{"superPop"}.
#' 
#' @slot fieldPopInfAnc a \code{character} string representing the name of 
#' the column in the data frame that contains the inferred 
#' super-population ancestry for the samples. The column should be 
#' present in the data frame. Default: \code{"superPop"}.
#' 
#' @slot fieldSubPop a \code{character} string representing the name of 
#' the column in the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter) that contains the sub-population 
#' information for the samples. The column should be present in the Population 
#' Reference GDS file. Default: \code{"pop.group"}.
#' 
#' @slot anchorId a \code{character} string representing the name of the
#' identifier of the quasi non-admixed profile from the reference. 
#' Default: \code{"sp0.95v0.2"}.
#' 
#' @slot pathADMIXTURE a \code{character} string representing the path pf
#' the software ADMIXTURE and the binary. If it use a function using ADMIXTURE, then
#' if it is null, it will try to find it with \code{Sys.which("admixture")}.
#' Default: \code{NULL}.
#' 
#' @slot seedADMIXTURE a single \code{integer} seed for
#' ADMIXTURE program. Default: \code{432L}.  
#' 
#' @slot cleanADMIXTUREFiles a single \code{logical} indicating if RAIDS
#' removes ADMIXTURE files after processing. Default: \code{TRUE}.  
#' 
#' @slot verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{FALSE}.
#' 
#' 
#' @details
#'
#' The `chrInfo` slot contains the length of the chromosomes. The
#' length of the chromosomes can be obtain through the
#' \code{\link[Seqinfo]{seqlengths}} library.
#'
#' As example, for hg38 genome:
#'
#' ```
#'
#' if (requireNamespace("Seqinfo", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'      chrInfo <- Seqinfo::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#' }
#'
#' ```
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @aliases RAIDSparam-class
#' @name RAIDSparam-class
#' @rdname RAIDSparam-class
#' 
#' @keywords classes
#' @exportClass RAIDSparam
#' @importFrom Seqinfo seqlengths
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 Hsapiens
#' @importFrom stringr str_detect
#' @export
setClass("RAIDSparam",
  slots = c(
    studyDF = "data.frame",
    studyDFSyn = "data.frame",
    pedStudy="data.frame",
    studyType="character",
    genoSource="CharacterOrNULL",
    blockTypeId="character",
    reference="character",
    genome="character",
    chrInfo="integer",
    paramAncestry="list",
    profileFile="CharacterOrNULL",
    profileFileGeno="CharacterOrNULL",
    pathProfileGDS="CharacterOrNULL",
    fileReferenceGDS="CharacterOrNULL",
    fileReferenceAnnotGDS="CharacterOrNULL",
    inferenceType="character",
    sampleRef="CharacterOrNULL",
    batch="integer",
    prefix="character",
    nbSim="integer",
    posOffset="integer",
    minCov="integer",
    minProb="numeric",
    seqError="numeric",
    seqErrorSyn="numeric",
    pRecomb="numeric",
    np="integer",
    listPos="DataFrameOrNULL",
    syntheticRefDF="DataFrameOrNULL",
    pruningMethod="character",
    slideWindowMaxBP="integer",
    thresholdLD="numeric",
    specificSNV="DataFrameOrNULL",
    genoType="character",
    phaseType="character",
    phase="logical",
    PCAmissingRate="numeric",
    PCAalgorithm="character",
    eigenCount="integer",
    eigenCountSyn="integer",
    kList="integer",
    pcaList="integer",
    fieldPopInRef="character",
    fieldPopInfAnc="character",
    fieldSubPop="character",
    anchorId="character",
    pathADMIXTURE="CharacterOrNULL",
    seedADMIXTURE="integer",
    cleanADMIXTUREFiles="logical",
    verbose="logical"
  ),
  prototype = list(
    studyDF=data.frame(study.id="NotDef",
                                study.desc="NotDef",
                                study.platform="NotDef",
                                stringsAsFactors=FALSE),
    studyDFSyn=data.frame(study.id="NotDef.Synthetic",
                        study.desc="NotDef synthetic data",
                        study.platform="Synthetic", stringsAsFactors=FALSE),
    pedStudy=data.frame(Name.ID=c("ProfileId"),
                                Case.ID=c("ProfileId"),
                                Sample.Type=c("type"),
                                Diagnosis="NotDef",
                                Source=c("NotDef"),
                                stringsAsFactors=FALSE, 
                                row.names=c("ProfileId")),
    studyType="LD",
    genoSource=NULL,
    blockTypeId="GeneS.Ensembl.Hsapiens.v86",
    reference="1KGv1.0",
    genome="HG38",
    chrInfo=seqlengths(Hsapiens)[seq_len(25)],
    paramAncestry=list(ScanBamParam=NULL, PileupParam=NULL,
                                yieldSize=10000000),
    profileFile=NULL,
    profileFileGeno=NULL,
    pathProfileGDS=NULL,
    fileReferenceGDS=NULL,
    fileReferenceAnnotGDS=NULL,
    inferenceType="PCAknn",
    sampleRef=NULL,
    batch=1L,
    prefix="1",
    nbSim=1L,
    posOffset=-1L,
    minCov=10L,
    minProb=0.999,
    seqError=0.001,
    seqErrorSyn=0.001,
    pRecomb=0.01,
    np=1L,
    listPos=NULL,
    syntheticRefDF=NULL,
    pruningMethod="corr",
    slideWindowMaxBP=500000L,
    thresholdLD=sqrt(0.1),
    specificSNV=NULL,
    genoType="geno.ref",
    phaseType="phase.ref",
    phase=FALSE,
    PCAmissingRate=0.025,
    PCAalgorithm="exact",
    eigenCount=32L,
    eigenCountSyn=15L,
    kList=seq(2L, 15L, 1L),
    pcaList=seq(2L, 15L, 1L),
    fieldPopInRef="superPop",
    fieldPopInfAnc="superPop",
    fieldSubPop="pop.group",
    anchorId="sp0.95v0.2",
    pathADMIXTURE=NULL,
    seedADMIXTURE=432L,
    cleanADMIXTUREFiles=TRUE,
    verbose=FALSE
  )
)

## Validation
setValidity("RAIDSparam",
    function(object)
    {
        ## Validate the studyDF parameter
        if (!(is.data.frame(object@studyDF) && 
                all(c("study.id", "study.desc", "study.platform") %in% 
                    colnames(object@studyDF)))) {
            return(paste0("'studyDF' slot must be a data.frame with ", 
                "those 3 columns: \"study.id\", \"study.desc\", ", 
                "\"study.platform\"."))
        }
      
        ## Validate the studyDFSyn parameter
        if (!(is.data.frame(object@studyDFSyn) && 
                all(c("study.id", "study.desc", "study.platform") %in% 
                    colnames(object@studyDFSyn)))) {
            return(paste0("'studyDFSyn' slot must be a data.frame ", 
                "with those 3 columns: \"study.id\", \"study.desc\", ", 
                "\"study.platform\"."))
        }

        ## Validate the pedStudy parameter
        if (!(is.data.frame(object@pedStudy) && 
                all(c("Name.ID", "Case.ID", "Sample.Type", "Diagnosis", 
                        "Source") %in% colnames(object@pedStudy)) && 
                all(rownames(object@pedStudy) == object@pedStudy$Name.ID))) {
            return(paste0("'pedStudy' slot must be a data.frame ", 
                "with those 5 columns: \"Name.ID\", \"Case.ID\", ", 
                "\"Sample.Type\", \"Diagnosis\", and \"Source\". All row ", 
                "names should correspond to the Name.ID values."))
        }
      
        ## Validate the studyType parameter
        if (length(object@studyType) != 1 || 
                !object@studyType %in% c("LD", "GeneAware")) {
            return(paste0("'studyType' slot must have one character", 
                " string within those 2 choices: \"LD\" and \"GeneAware\"."))
        }
      
        ## Validate the genoSource parameter 
        if (!is.null(object@genoSource) && !(length(object@genoSource) == 1 && 
            object@genoSource %in% c("VCF", "bam", "generic", "snp-pileup"))) {
            return(paste0("'genoSource' slot must have one character ", 
                "string or NULL. The valid options are: \"VCF\", \"bam\", ", 
                "\"generic\", or \"snp-pileup\"."))
        }
      
        ## Validate thte blockTypeId parameter
        if (length(object@blockTypeId) != 1) {
            return("'blockTypeId' slot must have one character string.")
        }

        ## Validate the reference parameter
        if (length(object@reference) != 1 || 
            !object@reference %in% c("1KGv1.0", "1k_hgdpV0.1")) {
            return(paste0("'reference' slot must have one character string", 
                " within those 2 choices: \"1KGv1.0\", \"1k_hgdpV0.1\"."))
        }

        ## Validate the genome parameter
        if (length(object@genome) != 1 || object@genome != "HG38") {
            return("'genome' slot must be the character string \"HG38\".")
        }

        ## Validate the chrInfo parameter TODO

        ## Validate the paramAncestry parameter 
        if (length(object@paramAncestry) != 3 || 
                    !all(c("ScanBamParam", "PileupParam", "yieldSize") %in% 
                        names(object@paramAncestry))) {
            return(paste0("'paramAncestry' slot must be a list with those ", 
                        "three entries: \"ScanBamParam\", \"PileupParam\", ", 
                        "and \"yieldSize\"."))
        }
      
        ## Validate the profileFile parameter
        if (!is.null(object@profileFile) && length(object@profileFile) != 1) {
            return("'profileFile' slot must have one character string.")
        } 

        ## Validate the profileFile parameter using genoSource information
        if (!is.null(object@profileFile) && length(object@profileFile) == 1) {
            genoSourceVal <- object@genoSource
            if (is.null(genoSourceVal)) {
                return(paste0("'genoSource' slot must be coherent with the ", 
                    "'profileFile' slot. The 'genoSource' slot is currently ", 
                    "NULL."))
            } else if (genoSourceVal == "bam" && 
                    !stringr::str_detect(object@profileFile, ".bam$")) {
                return(paste0("'profileFile' slot must have one character ", 
                    "string representing a file with extension '.bam' ", 
                    "according to 'genoSource' slot."))
            } else if (genoSourceVal == "VCF" && 
                    !stringr::str_detect(object@profileFile, ".vcf.gz$")) {
                return(paste0("'profileFile' slot must have one character ", 
                    "string representing a file with extension '.vcf.gz' ", 
                    "according to 'genoSource' slot."))
            } else if (genoSourceVal %in% c("snp-pileup", "generic") && 
                    !stringr::str_detect(object@profileFile, ".txt.gz$")) {
                return(paste0("'profileFile' slot must have one character ", 
                    "string representing a file with extension '.txt.gz' ", 
                    "according to 'genoSource' slot."))
            }
        } 
        
        ## Validate the profileFile file exists when not null
        if (!is.null(object@profileFile) && !file.exists(object@profileFile)) {
            return(paste0("'profileFile' slot must have one character ", 
                    "string representing an existing file."))
        }
      
        ## Validate the profileFileGeno parameter TODO
      
        ## Validate the pathProfileGDS parameter
        if ((!is.null(object@pathProfileGDS) && 
                    length(object@pathProfileGDS) != 1) || 
            (!is.null(object@pathProfileGDS) && 
                    !dir.exists(object@pathProfileGDS))) {
            return(paste0("'pathProfileGDS' slot must have one character ", 
                    "string representing an existing directory."))
        }
      
        ## Validate the fileReferenceGDS parameter 
        if ((!is.null(object@fileReferenceGDS) && 
                    length(object@fileReferenceGDS) != 1) || 
            (!is.null(object@fileReferenceGDS) && 
                    !file.exists(object@fileReferenceGDS))) {
            return(paste0("'fileReferenceGDS' slot must have one character ", 
                    "string representing an existing file."))
        }
      
        ## Validate the fileReferenceAnnotGDS parameter
        if ((!is.null(object@fileReferenceAnnotGDS) && 
                    length(object@fileReferenceAnnotGDS) != 1) || 
            (!is.null(object@fileReferenceAnnotGDS) && 
                    !file.exists(object@fileReferenceAnnotGDS))) {
            return(paste0("'fileReferenceAnnotGDS' slot must have one ", 
                    "character string representing an existing file."))
        }
      
        ## Validate the inferenceType parameter 
        if (length(object@inferenceType) != 1 || 
                !object@inferenceType %in% c("PCAknn", "haploAdmixture")) {
            return(paste0("'inferenceType' slot must be a single character" , 
                " string. The valid options are: 'PCAknn' and ", 
                "'haploAdmixture'."))
        }

        ## Validate the sampleRef parameter (NOTHING TO DO HERE)
      
        ## Validate the batch parameter
        if (length(object@batch) != 1 || object@batch < 1) {
            return("'batch' slot must have one positive integer.")
        }

        ## Validate the prefix parameter
        if (length(object@prefix) != 1) {
            return("'prefix' slot must have one character string.")
        }

        ## Validate the nbSim parameter
        if ( length(object@nbSim) != 1 || is.na(object@nbSim) || object@nbSim != 1L) { # || object@nbSim < 1L
            return("'nbSim' slot must be 1 for the moment.")# have one positive integer
        }

        ## Validate the posOffset parameter
        if (length(object@posOffset) != 1) {
            return("'posOffset' slot must have one integer.")
        }
        
        ## Validate the minCov parameter
        if (length(object@minCov) != 1 || object@minCov < 1L) {
            return("'minCov' slot must have one positive integer.")
        }

        ## Validate the minProb parameter
        if (length(object@minProb) != 1 || 
                !(object@minProb > 0 && object@minProb < 1)) {
            return(paste0("'minProb' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the seqError parameter
        if (length(object@seqError) != 1 || 
                !(object@seqError > 0 && object@seqError < 1)) {
            return(paste0("'seqError' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the seqErrorSyn parameter
        if (length(object@seqErrorSyn) != 1 || 
                !(object@seqErrorSyn > 0 && object@seqErrorSyn < 1)) {
            return(paste0("'seqErrorSyn' slot must have one positive numeric", 
                        " between 0 and 1."))
        }

        ## Validate the pRecomb parameter
        if (length(object@pRecomb) != 1 || 
                !(object@pRecomb > 0 && object@pRecomb < 1)) {
            return(paste0("'pRecomb' slot must have one positive numeric", 
                        " between 0 and 1."))
        }
      
        ## Validate the np parameter
        if (length(object@np) != 1 || 
                !(object@np > 0)) {
            return("'np' slot must have one positive integer.")
        }

        ## Validate the listPos parameter
        if (!(is.null(object@listPos) || (ncol(object@listPos) == 2 && 
                all(colnames(object@listPos) %in% c("snp.chromosome", 
                "snp.position"))))) {
            return(paste0("'listPos' slot must be NULL or a data.frame with ", 
                "2 columns named \"snp.chromosome\" and \"snp.position\"."))
        }

        ## Validate the syntheticRefDF parameter
        if (!(is.null(object@syntheticRefDF) || 
                ncol(object@syntheticRefDF) >= 3 && 
                    all(c("sample.id", "pop.group", "superPop") %in% 
                            colnames(object@syntheticRefDF)))) {
            return(paste0("'syntheticRefDF' slot must be NULL or a ", 
                "data.frame with 3 columns named \"sample.id\", ", 
                "\"pop.group\", and \"superPop\"."))
        }

        ## Validate the pruningMethod parameter
        if (length(object@pruningMethod) != 1 || 
              !object@pruningMethod %in% c("corr", "r", "dprime",
                  "composite")) {
          return(paste0("'pruningMethod' slot must have one character string ", 
              "within those 4 choices: \"corr\", \"r\", \"dprime\", ", 
              "\"composite\"."))
        }

        ## Validate the slideWindowMaxBP parameter
        if (length(object@slideWindowMaxBP) != 1 || 
                            object@slideWindowMaxBP < 0) {
          return(paste0("'slideWindowMaxBP' slot must have one positive ", 
                "integer value."))
        }
      
        ## Validate the thresholdLD parameter
        if (length(object@thresholdLD) != 1 || object@thresholdLD < 0) {
          return("'thresholdLD' slot must have one positive numeric value.")
        }
      
        ## Validate the specificSNV parameter
        if (!(is.null(object@specificSNV) || 
                (ncol(object@specificSNV) >= 2 && 
                    all(c("snp.chromosome", "snp.position") %in% 
                            colnames(object@specificSNV) )))) {
            return(paste0("'specificSNV' slot must be NULL or a ", 
                "data.frame with 2 columns named \"snp.chromosome\"", 
                " and \"snp.position\"."))
        }

        ## Validate the genoType parameter
        if (length(object@genoType) != 1) {
          return("'genoType' slot must have one character string.")
        }
      
        ## Validate the phaseType parameter
        if (length(object@phaseType) != 1) {
            return("'phaseType' slot must be one character string.")
        }
        
        ## Validate the phase parameter
        if (length(object@phase) != 1) {
          return("'phase' slot must have one logical value.")
        }

        ## Validate the PCAmissingRate parameter
        if (length(object@PCAmissingRate) != 1 || object@PCAmissingRate < 0) {
          return("'PCAmissingRate' slot must have one positive numeric value.")
        }
      
        ## Validate the PCAalgorithm parameter
        if (length(object@PCAalgorithm) != 1 || 
              !object@PCAalgorithm %in% c("exact", "randomized")) {
          return(paste0("'PCAalgorithm' slot must have one character string. ", 
              "The valid options are: \"exact\" or \"randomized\"."))
        }
      
        ## Validate the eigenCount parameter
        if (length(object@eigenCount) != 1) {
          return("'eigenCount' slot must have one integer value.")
        }
      
        ## Validate the eigenCountSyn parameter 
        if (length(object@eigenCountSyn) != 1) {
          return("'eigenCountSyn' slot must have one integer value.")
        }

        ## Validate the kList parameter
        if (!all(object@kList > 0)) {
          return(paste0("'kList' slot must have one or more positive ", 
                    "integer values."))
        }

        ## Validate the pcaList parameter
        if (!all(object@pcaList > 0)) {
          return(paste0("'pcaList' slot must have one or more positive ", 
                    "integer values."))
        }

        ## Validate the fieldPopInRef parameter
        if (length(object@fieldPopInRef) != 1) {
          return("'fieldPopInRef' slot must have one character string.")
        }

        ## Validate the fieldPopInfAnc parameter
        if (length(object@fieldPopInfAnc) != 1) {
          return("'fieldPopInfAnc' slot must have one character string.")
        }

        ## Validate the fieldSubPop parameter
        if (length(object@fieldSubPop) != 1) {
          return("'fieldSubPop' slot must have one character string.")
        }
        
        ## Validate the anchorId parameter
        if (length(object@anchorId) != 1) {
          return("'anchorId' slot must have one character string.")
        }
        
        ## Validate the pathADMIXTURE parameter
        if ((!is.null(object@pathADMIXTURE) && 
                    length(object@pathADMIXTURE) != 1) || 
            (!is.null(object@pathADMIXTURE) && 
                    !file.exists(object@pathADMIXTURE))) {
            return(paste0("'pathADMIXTURE' slot must have one character ", 
                    "string representing an existing file."))
        }

        ## Validate the seedADMIXTURE parameter 
        if (length(object@seedADMIXTURE) != 1) {
          return("'seedADMIXTURE' slot must have one integer value.")
        }
        
        ## Validate the cleanADMIXTUREFiles parameter
        if (length(object@cleanADMIXTUREFiles) != 1) {
          return("'cleanADMIXTUREFiles' slot must have one logical value.")
        }
      
        ## Validate the verbose parameter
        if (length(object@verbose) != 1) {
          return("'verbose' slot must have one logical value.")
        }

        TRUE
    }
)

###########################################################################
## All the getter functions for the RAIDSparam class
###########################################################################

#' Generic function for getting the studyDF slot in a class
#' 
#' @description A generic function for getting the studyDF slot in a 
#' class. The function is implemented for the \code{RAIDSparam} class.
#' 
#' @param x an object.
#' 
#' @return an object.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @export
setGeneric("studyDF", function(x) standardGeneric("studyDF"))

#' Getter for the studyDF slot in a RAIDSparam class
#' 
#' @description A function for getting the studyDF slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{data.frame}.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the studyDF slot for the object
#' studyDF(paramDemo)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @export
setMethod("studyDF", "RAIDSparam", function(x) {
  return(x@studyDF)
})


#' Generic function for getting the studyDFSyn slot in a class
#' 
#' @description A generic function for getting the studyDFSyn slot in a 
#' class. The function is implemented for the \code{RAIDSparam} class.
#' 
#' @param x an object.
#' 
#' @return an object.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @export
setGeneric("studyDFSyn", function(x) standardGeneric("studyDFSyn"))


#' A getter for the studyDFSyn slot in a RAIDSparam class
#' 
#' @description A function for getting the studyDFSyn slot in a 
#' \code{RAIDSparam} class.
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{data.frame}.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the studyDFSyn slot for the object
#' studyDFSyn(paramDemo)
#' 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @export
setMethod("studyDFSyn", "RAIDSparam", function(x) {
  return(x@studyDFSyn)
})


#' Generic function for getting the pedStudy slot in a class
#' 
#' @description A generic function for getting the pedStudy slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the pedStudy slot in the S4 object.
#' 
#' @examples
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(studyDF="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", studyDF="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # studyDF(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pedStudy", function(x) standardGeneric("pedStudy"))


#' A getter for the pedStudy slot in a RAIDSparam class
#' 
#' @description A function for getting the pedStudy slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{data.frame}.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the pedStudy slot for the object
#' pedStudy(paramDemo)
#' 
#' @export
setMethod("pedStudy", "RAIDSparam", function(x) {
  return(x@pedStudy)
})


#' Generic function for getting the studyType slot in a class
#' 
#' @description A generic function for getting the studyType slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the studyType slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(studyType="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", studyType="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # studyType(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("studyType", function(x) standardGeneric("studyType"))


#' A getter for the studyType slot in a RAIDSparam class
#' 
#' @description A function for getting the studyType slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string representing the type 
#' of study.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the studyType slot for the object
#' studyType(paramDemo)
#' 
#' @export
setMethod("studyType", "RAIDSparam", function(x) {
  return(x@studyType)
})


#' Generic function for getting the genoSource slot in a class
#' 
#' @description A generic function for getting the genoSource slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the genoSource slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(genoSource="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", genoSource="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # genoSource(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("genoSource", function(x) standardGeneric("genoSource"))


#' A getter for the genoSource slot in a RAIDSparam class
#' 
#' @description A function for getting the genoSource slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string corresponding to the type 
#' of file with the genotype and the allele information of the profile that 
#' will be provided in the 'profileFile' slot. The valid options are: "VCF", 
#' "generic", "snp-pileup", and "bam". 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the genoSource slot for the object
#' genoSource(paramDemo)
#' 
#' @export
setMethod("genoSource", "RAIDSparam", function(x) {
  return(x@genoSource)
})


#' Generic function for getting the blockTypeId slot in a class
#' 
#' @description A generic function for getting the blockTypeId slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the blockTypeId slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(blockTypeId="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", blockTypeId="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # blockTypeId(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("blockTypeId", function(x) standardGeneric("blockTypeId"))


#' A getter for the blockTypeId slot in a RAIDSparam class
#' 
#' @description A function for getting the blockTypeId slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string corresponding to the block type 
#' used to extract the block identifiers. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the blockTypeId slot for the object
#' blockTypeId(paramDemo)
#' 
#' @export
setMethod("blockTypeId", "RAIDSparam", function(x) {
  return(x@blockTypeId)
})


#' Generic function for getting the reference slot in a class
#' 
#' @description A generic function for getting the reference slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the reference slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(reference="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", reference="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # reference(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("reference", function(x) standardGeneric("reference"))


#' A getter for the reference slot in a RAIDSparam class
#' 
#' @description A function for getting the reference slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string with two possible values: 
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the reference slot for the object
#' reference(paramDemo)
#' 
#' @export
setMethod("reference", "RAIDSparam", function(x) {
  return(x@reference)
})


#' Generic function for getting the genome slot in a class
#' 
#' @description A generic function for getting the genome slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the genome slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(genome="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", genome="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # genome(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("genome", function(x) standardGeneric("genome"))


#' A getter for the genome slot in a RAIDSparam class
#' 
#' @description A function for getting the genome slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return  \code{character} string with one possible value: 'HG38'.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the genome slot for the object
#' genome(paramDemo)
#' 
#' @export
setMethod("genome", "RAIDSparam", function(x) {
  return(x@genome)
})


#' Generic function for getting the chrInfo slot in a class
#' 
#' @description A generic function for getting the chrInfo slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the chrInfo slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(chrInfo="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", chrInfo="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # chrInfo(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("chrInfo", function(x) standardGeneric("chrInfo"))


#' A getter for the chrInfo slot in a RAIDSparam class
#' 
#' @description A function for getting the chrInfo slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{vector} of positive \code{integer} values 
#' representing the length of the chromosomes.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the chrInfo slot for the object
#' chrInfo(paramDemo)
#' 
#' @export
setMethod("chrInfo", "RAIDSparam", function(x) {
  return(x@chrInfo)
})


#' Generic function for getting the paramAncestry slot in a class
#' 
#' @description A generic function for getting the paramAncestry slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the paramAncestry slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(paramAncestry="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", paramAncestry="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # paramAncestry(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("paramAncestry", function(x) standardGeneric("paramAncestry"))


#' A getter for the paramAncestry slot in a RAIDSparam class
#' 
#' @description A function for getting the paramAncestry slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{list} of parameters related to ancestry.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the paramAncestry slot for the object
#' paramAncestry(paramDemo)
#' 
#' @export
setMethod("paramAncestry", "RAIDSparam", function(x) {
  return(x@paramAncestry)
})


#' Generic function for getting the profileFile slot in a class
#' 
#' @description A generic function for getting the profileFile slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the profileFile slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(profileFile="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", profileFile="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # profileFile(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("profileFile", function(x) standardGeneric("profileFile"))


#' A getter for the profileFile slot in a RAIDSparam class
#' 
#' @description A function for getting the profileFile slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return \code{NULL} or a \code{character} string representing the path to 
#' the file with genotype and the allele information of the profile.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the profileFile slot for the object
#' profileFile(paramDemo)
#' 
#' @export
setMethod("profileFile", "RAIDSparam", function(x) {
  return(x@profileFile)
})


#' Generic function for getting the pathProfileGDS slot in a class
#' 
#' @description A generic function for getting the pathProfileGDS slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the pathProfileGDS slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pathProfileGDS="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pathProfileGDS="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pathProfileGDS(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pathProfileGDS", function(x) standardGeneric("pathProfileGDS"))


#' A getter for the pathProfileGDS slot in a RAIDSparam class
#' 
#' @description A function for getting the pathProfileGDS slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return \code{NULL} or a single \code{character} string representing the 
#' path to the directory where the GDS Profile files will be created.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the pathProfileGDS slot for the object
#' pathProfileGDS(paramDemo)
#' 
#' @export
setMethod("pathProfileGDS", "RAIDSparam", function(x) {
  return(x@pathProfileGDS)
})


#' Generic function for getting the fileReferenceGDS slot in a class
#' 
#' @description A generic function for getting the fileReferenceGDS slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the fileReferenceGDS slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fileReferenceGDS="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fileReferenceGDS="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fileReferenceGDS(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fileReferenceGDS", function(x) standardGeneric("fileReferenceGDS"))


#' A getter for the fileReferenceGDS slot in a RAIDSparam class
#' 
#' @description A function for getting the fileReferenceGDS slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return \code{NULL} or a single \code{character} string representing 
#' the file name of the Reference GDS file.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the fileReferenceGDS slot for the object
#' fileReferenceGDS(paramDemo)
#' 
#' @export
setMethod("fileReferenceGDS", "RAIDSparam", function(x) {
  return(x@fileReferenceGDS)
})


#' Generic function for getting the fileReferenceAnnotGDS slot in a class
#' 
#' @description A generic function for getting the fileReferenceAnnotGDS slot 
#' in a S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the fileReferenceAnnotGDS slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fileReferenceAnnotGDS="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fileReferenceAnnotGDS="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fileReferenceAnnotGDS(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fileReferenceAnnotGDS", 
    function(x) standardGeneric("fileReferenceAnnotGDS"))


#' A getter for the fileReferenceAnnotGDS slot in a RAIDSparam class
#' 
#' @description A function for getting the fileReferenceAnnotGDS slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return \code{NULL} or a \code{character} string representing the file 
#' name of the Population Reference GDS Annotation file.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the fileReferenceAnnotGDS slot for the object
#' fileReferenceAnnotGDS(paramDemo)
#' 
#' @export
setMethod("fileReferenceAnnotGDS", "RAIDSparam", function(x) {
  return(x@fileReferenceAnnotGDS)
})



#' Generic function for getting the inferenceType slot in a class
#' 
#' @description A generic function for getting the inferenceType slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the inferenceType slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(inferenceType="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", inferenceType="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # inferenceType(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("inferenceType", 
    function(x) standardGeneric("inferenceType"))


#' A getter for the inferenceType slot in a RAIDSparam class
#' 
#' @description A function for getting the inferenceType slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string representing the genotype 
#' ancestry inference method.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the inferenceType slot for the object
#' inferenceType(paramDemo)
#' 
#' @export
setMethod("inferenceType", "RAIDSparam", function(x) {
  return(x@inferenceType)
})


#' Generic function for getting the sampleRef slot in a class
#' 
#' @description A generic function for getting the sampleRef slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the sampleRef slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(sampleRef="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", sampleRef="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # sampleRef(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("sampleRef", function(x) standardGeneric("sampleRef"))


#' A getter for the inferenceType slot in a RAIDSparam class
#' 
#' @description A function for getting the sampleRef slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return \code{NULL} or a \code{vector} of \code{character} strings that 
#' represent a subset of sample identifiers from the GDS Reference File. 
#' Those samples will be used for the ancestry inference. If \code{NULL}, all 
#' the samples are retained.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the sampleRef slot for the object
#' sampleRef(paramDemo)
#' 
#' @export
setMethod("sampleRef", "RAIDSparam", function(x) {
  return(x@sampleRef)
})


#' Generic function for getting the prefix slot in a class
#' 
#' @description A generic function for getting the prefix slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the prefix slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(prefix="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", prefix="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # prefix(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("prefix", function(x) standardGeneric("prefix"))


#' A getter for the prefix slot in a RAIDSparam class
#' 
#' @description A function for getting the prefix slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string representing TODO.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the prefix slot for the object
#' prefix(paramDemo)
#' 
#' @export
setMethod("prefix", "RAIDSparam", function(x) {
  return(x@prefix)
})


#' Generic function for getting the batch slot in a class
#' 
#' @description A generic function for getting the batch slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the batch slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(batch="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", batch="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # batch(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("batch", function(x) standardGeneric("batch"))


#' A getter for the batch slot in a RAIDSparam class
#' 
#' @description A function for getting the batch slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single positive \code{integer} representing the current 
#' identifier for the batch.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the batch slot for the object
#' batch(paramDemo)
#' 
#' @export
setMethod("batch", "RAIDSparam", function(x) {
  return(x@batch)
})


#' A getter for the prefix slot in a RAIDSparam class
#' 
#' @description A function for getting the prefix slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string representing TODO.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the prefix slot for the object
#' prefix(paramDemo)
#' 
#' @export
setMethod("prefix", "RAIDSparam", function(x) {
  return(x@prefix)
})


#' Generic function for getting the nbSim slot in a class
#' 
#' @description A generic function for getting the nbSim slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the nbSim slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(nbSim="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", nbSim="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # nbSim(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("nbSim", function(x) standardGeneric("nbSim"))


#' A getter for the nbSim slot in a RAIDSparam class
#' 
#' @description A function for getting the nbSim slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a positive \code{integer} representing TODO.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the nbSim slot for the object
#' nbSim(paramDemo)
#' 
#' @export
setMethod("nbSim", "RAIDSparam", function(x) {
  return(x@nbSim)
})


#' A getter for the prefix slot in a RAIDSparam class
#' 
#' @description A function for getting the prefix slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string representing TODO.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the prefix slot for the object
#' prefix(paramDemo)
#' 
#' @export
setMethod("prefix", "RAIDSparam", function(x) {
  return(x@prefix)
})


#' Generic function for getting the minCov slot in a class
#' 
#' @description A generic function for getting the minCov slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the minCov slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(minCov="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", minCov="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # minCov(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("minCov", function(x) standardGeneric("minCov"))


#' A getter for the minCov slot in a RAIDSparam class
#' 
#' @description A function for getting the minCov slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single positive \code{integer} representing the minimum 
#' required coverage. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the minCov slot for the object
#' minCov(paramDemo)
#' 
#' @export
setMethod("minCov", "RAIDSparam", function(x) {
  return(x@minCov)
})


#' Generic function for getting the posOffset slot in a class
#' 
#' @description A generic function for getting the offposOffsetset slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the posOffset slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(posOffset="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", posOffset="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # posOffset(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("posOffset", function(x) standardGeneric("posOffset"))


#' A getter for the posOffset slot in a RAIDSparam class
#' 
#' @description A function for getting the posOffset slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{integer} that is added to the SNP position to 
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the posOffset slot for the object
#' posOffset(paramDemo)
#' 
#' @export
setMethod("posOffset", "RAIDSparam", function(x) {
  return(x@posOffset)
})


#' Generic function for getting the minProb slot in a class
#' 
#' @description A generic function for getting the minProb slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the minProb slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(minProb="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", minProb="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # minProb(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("minProb", function(x) standardGeneric("minProb"))


#' A getter for the minProb slot in a RAIDSparam class
#' 
#' @description A function for getting the minProb slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{numeric} between \code{0} and \code{1} 
#' representing the probability that the calculated genotype call is correct.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the minProb slot for the object
#' minProb(paramDemo)
#' 
#' @export
setMethod("minProb", "RAIDSparam", function(x) {
  return(x@minProb)
})


#' Generic function for getting the seqError slot in a class
#' 
#' @description A generic function for getting the seqError slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the seqError slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(seqError="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", seqError="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # seqError(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("seqError", function(x) standardGeneric("seqError"))


#' A getter for the seqError slot in a RAIDSparam class
#' 
#' @description A function for getting the seqError slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{numeric} between \code{0} and \code{1} 
#' representing the probability of sequencing error.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the seqError slot for the object
#' seqError(paramDemo)
#' 
#' @export
setMethod("seqError", "RAIDSparam", function(x) {
  return(x@seqError)
})


#' Generic function for getting the seqErrorSyn slot in a class
#' 
#' @description A generic function for getting the seqErrorSyn slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the seqErrorSyn slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(seqErrorSyn="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", seqErrorSyn="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # seqErrorSyn(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("seqErrorSyn", function(x) standardGeneric("seqErrorSyn"))


#' A getter for the seqErrorSyn slot in a RAIDSparam class
#' 
#' @description A function for getting the seqErrorSyn slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{numeric} between \code{0} and \code{1} 
#' representing the probability of sequencing error for synthetic.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the seqErrorSyn slot for the object
#' seqErrorSyn(paramDemo)
#' 
#' @export
setMethod("seqErrorSyn", "RAIDSparam", function(x) {
  return(x@seqErrorSyn)
})


#' Generic function for getting the pRecomb slot in a class
#' 
#' @description A generic function for getting the pRecomb slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the pRecomb slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pRecomb="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pRecomb="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pRecomb(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pRecomb", function(x) standardGeneric("pRecomb"))


#' A getter for the pRecomb slot in a RAIDSparam class
#' 
#' @description A function for getting the pRecomb slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single positive \code{numeric} between \code{0} and \code{1} 
#' that represents the frequency of phase switching in the synthetic profiles.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the pRecomb slot for the object
#' pRecomb(paramDemo)
#' 
#' @export
setMethod("pRecomb", "RAIDSparam", function(x) {
  return(x@pRecomb)
})


#' Generic function for getting the np slot in a class
#' 
#' @description A generic function for getting the np slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the np slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(np="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", np="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # np(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("np", function(x) standardGeneric("np"))


#' A getter for the np slot in a RAIDSparam class
#' 
#' @description A function for getting the np slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single positive \code{integer} specifying the number of 
#' threads to be used.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the np slot for the object
#' np(paramDemo)
#' 
#' @export
setMethod("np", "RAIDSparam", function(x) {
  return(x@np)
})


#' Generic function for getting the listPos slot in a class
#' 
#' @description A generic function for getting the listPos slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the listPos slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(listPos="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", listPos="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # listPos(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("listPos", function(x) standardGeneric("listPos"))


#' A getter for the listPos slot in a RAIDSparam class
#' 
#' @description A function for getting the listPos slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single positive \code{integer} specifying the number of 
#' threads to be used.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the listPos slot for the object
#' listPos(paramDemo)
#' 
#' @export
setMethod("listPos", "RAIDSparam", function(x) {
  return(x@listPos)
})


#' Generic function for getting the syntheticRefDF slot in a class
#' 
#' @description A generic function for getting the syntheticRefDF slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the syntheticRefDF slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(syntheticRefDF="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", syntheticRefDF="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # syntheticRefDF(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("syntheticRefDF", function(x) standardGeneric("syntheticRefDF"))


#' A getter for the syntheticRefDF slot in a RAIDSparam class
#' 
#' @description A function for getting the syntheticRefDF slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return \code{NULL} or a \code{data.frame} containing a subset of 
#' reference profiles for each sub-population present in the Reference GDS 
#' file.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the syntheticRefDF slot for the object
#' syntheticRefDF(paramDemo)
#' 
#' @export
setMethod("syntheticRefDF", "RAIDSparam", function(x) {
  return(x@syntheticRefDF)
})


#' Generic function for getting the pruningMethod slot in a class
#' 
#' @description A generic function for getting the pruningMethod slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the pruningMethod slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pruningMethod="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pruningMethod="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pruningMethod(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pruningMethod", function(x) standardGeneric("pruningMethod"))


#' A getter for the pruningMethod slot in a RAIDSparam class
#' 
#' @description A function for getting the pruningMethod slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string representing the method that 
#' will be used to calculate the linkage disequilibrium in the 
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the pruningMethod slot for the object
#' pruningMethod(paramDemo)
#' 
#' @export
setMethod("pruningMethod", "RAIDSparam", function(x) {
  return(x@pruningMethod)
})


#' Generic function for getting the slideWindowMaxBP slot in a class
#' 
#' @description A generic function for getting the slideWindowMaxBP slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the slideWindowMaxBP slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(slideWindowMaxBP="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", slideWindowMaxBP="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # slideWindowMaxBP(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("slideWindowMaxBP", function(x) standardGeneric("slideWindowMaxBP"))

#' A getter for the slideWindowMaxBP slot in a RAIDSparam class
#' 
#' @description A function for getting the slideWindowMaxBP slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single positive \code{integer} that represents the maximum 
#' basepairs (bp) in the sliding window. This parameter is used for 
#' the LD pruning done by the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the slideWindowMaxBP slot for the object
#' slideWindowMaxBP(paramDemo)
#' 
#' @export
setMethod("slideWindowMaxBP", "RAIDSparam", function(x) {
  return(x@slideWindowMaxBP)
})



#' A getter for the pruningMethod slot in a RAIDSparam class
#' 
#' @description A function for getting the pruningMethod slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string representing the method that 
#' will be used to calculate the linkage disequilibrium in the 
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the pruningMethod slot for the object
#' pruningMethod(paramDemo)
#' 
#' @export
setMethod("pruningMethod", "RAIDSparam", function(x) {
  return(x@pruningMethod)
})


#' Generic function for getting the thresholdLD slot in a class
#' 
#' @description A generic function for getting the thresholdLD slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the thresholdLD slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(thresholdLD="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", thresholdLD="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # thresholdLD(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("thresholdLD", function(x) standardGeneric("thresholdLD"))


#' A getter for the thresholdLD slot in a RAIDSparam class
#' 
#' @description A function for getting the thresholdLD slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single positive \code{numeric} value that represents 
#' the LD threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the thresholdLD slot for the object
#' thresholdLD(paramDemo)
#' 
#' @export
setMethod("thresholdLD", "RAIDSparam", function(x) {
  return(x@thresholdLD)
})


#' Generic function for getting the specificSNV slot in a class
#' 
#' @description A generic function for getting the specificSNV slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the specificSNV slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(specificSNV="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", specificSNV="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # specificSNV(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("specificSNV", function(x) standardGeneric("specificSNV"))


#' A getter for the specificSNV slot in a RAIDSparam class
#' 
#' @description A function for getting the specificSNV slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return \code{NULL} or a \code{data.frame} containing 2 columns. The first 
#' column, called "snp.chromosome" contains the name of the chromosome where 
#' the SNV is located. The second column, called "snp.position", contains the 
#' position of the SNV on the chromosome. Optionally, the column "snvKeep" 
#' contains the SNV index in the reference or -1 if not in the reference. It 
#' is used during the pruning step. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the specificSNV slot for the object
#' specificSNV(paramDemo)
#' 
#' @export
setMethod("specificSNV", "RAIDSparam", function(x) {
  return(x@specificSNV)
})


#' Generic function for getting the genoType slot in a class
#' 
#' @description A generic function for getting the genoType slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the genoType slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(genoType="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", genoType="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # genoType(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("genoType", function(x) standardGeneric("genoType"))


#' A getter for the genoType slot in a RAIDSparam class
#' 
#' @description A function for getting the genoType slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string representing the field 
#' in the GDS profile file that contains all the genotype information. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the genoType slot for the object
#' genoType(paramDemo)
#' 
#' @export
setMethod("genoType", "RAIDSparam", function(x) {
  return(x@genoType)
})


#' Generic function for getting the phaseType slot in a class
#' 
#' @description A generic function for getting the phaseType slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the phaseType slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(phaseType="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", phaseType="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # phaseType(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("phaseType", function(x) standardGeneric("phaseType"))


#' A getter for the phaseType slot in a RAIDSparam class
#' 
#' @description A function for getting the phaseType slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{character} string TODO
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the phaseType slot for the object
#' phaseType(paramDemo)
#' 
#' @export
setMethod("phaseType", "RAIDSparam", function(x) {
  return(x@phaseType)
})


#' Generic function for getting the phase slot in a class
#' 
#' @description A generic function for getting the phase slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the phase slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(phase="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", phase="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # phase(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("phase", function(x) standardGeneric("phase"))


#' A getter for the phase slot in a RAIDSparam class
#' 
#' @description A function for getting the phase slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{logical}  TODO
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the phase slot for the object
#' phase(paramDemo)
#' 
#' @export
setMethod("phase", "RAIDSparam", function(x) {
  return(x@phase)
})


#' Generic function for getting the PCAmissingRate slot in a class
#' 
#' @description A generic function for getting the PCAmissingRate slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the PCAmissingRate slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(PCAmissingRate="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", PCAmissingRate="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # PCAmissingRate(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("PCAmissingRate", function(x) standardGeneric("PCAmissingRate"))


#' A getter for the PCAmissingRate slot in a RAIDSparam class
#' 
#' @description A function for getting the PCAmissingRate slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a positive \code{numeric} representing the maximum missing rate 
#' accepted to use SNPs in the PCA analysis done with the 
#' \link[SNPRelate]{snpgdsPCA} function. If \code{NaN}, no missing threshold.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the PCAmissingRate slot for the object
#' PCAmissingRate(paramDemo)
#' 
#' @export
setMethod("PCAmissingRate", "RAIDSparam", function(x) {
  return(x@PCAmissingRate)
})


#' Generic function for getting the PCAalgorithm slot in a class
#' 
#' @description A generic function for getting the PCAalgorithm slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the PCAalgorithm slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(PCAalgorithm="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", PCAalgorithm="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # PCAalgorithm(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("PCAalgorithm", function(x) standardGeneric("PCAalgorithm"))


#' A getter for the PCAalgorithm slot in a RAIDSparam class
#' 
#' @description A function for getting the PCAalgorithm slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{character} string representing the algorithm to use 
#' with the \link[SNPRelate]{snpgdsPCA} function. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the PCAalgorithm slot for the object
#' PCAalgorithm(paramDemo)
#' 
#' @export
setMethod("PCAalgorithm", "RAIDSparam", function(x) {
  return(x@PCAalgorithm)
})


#' Generic function for getting the eigenCount slot in a class
#' 
#' @description A generic function for getting the eigenCount slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the eigenCount slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(eigenCount="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", eigenCount="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # eigenCount(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("eigenCount", function(x) standardGeneric("eigenCount"))


#' A getter for the eigenCount slot in a RAIDSparam class
#' 
#' @description A function for getting the eigenCount slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{integer} indicating the number of 
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA} 
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the eigenCount slot for the object
#' eigenCount(paramDemo)
#' 
#' @export
setMethod("eigenCount", "RAIDSparam", function(x) {
  return(x@eigenCount)
})


#' Generic function for getting the eigenCountSyn slot in a class
#' 
#' @description A generic function for getting the eigenCountSyn slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the eigenCountSyn slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(eigenCountSyn="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", eigenCountSyn="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # eigenCountSyn(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("eigenCountSyn", function(x) standardGeneric("eigenCountSyn"))


#' A getter for the eigenCountSyn slot in a RAIDSparam class
#' 
#' @description A function for getting the eigenCountSyn slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{integer} indicating the number of 
#' eigenvectors TODO
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the eigenCountSyn slot for the object
#' eigenCountSyn(paramDemo)
#' 
#' @export
setMethod("eigenCountSyn", "RAIDSparam", function(x) {
  return(x@eigenCountSyn)
})


#' Generic function for getting the kList slot in a class
#' 
#' @description A generic function for getting the kList slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the kList slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(kList="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", kList="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # kList(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("kList", function(x) standardGeneric("kList"))


#' A getter for the kList slot in a RAIDSparam class
#' 
#' @description A function for getting the kList slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{vector} of positive \code{integer} representing the 
#' values tested for the  \emph{K} parameter. The \emph{K} parameter 
#' represents the number of neighbors used in the K-nearest neighbor analysis.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the kList slot for the object
#' kList(paramDemo)
#' 
#' @export
setMethod("kList", "RAIDSparam", function(x) {
  return(x@kList)
})


#' Generic function for getting the pcaList slot in a class
#' 
#' @description A generic function for getting the pcaList slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the pcaList slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pcaList="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pcaList="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pcaList(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pcaList", function(x) standardGeneric("pcaList"))


#' A getter for the pcaList slot in a RAIDSparam class
#' 
#' @description A function for getting the pcaList slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{vector} of positive \code{integer} representing the values 
#' tested for the  \emph{D} parameter. The \emph{D} parameter represents the 
#' number of dimensions used in the PCA analysis.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the pcaList slot for the object
#' pcaList(paramDemo)
#' 
#' @export
setMethod("pcaList", "RAIDSparam", function(x) {
  return(x@pcaList)
})


#' Generic function for getting the fieldPopInRef slot in a class
#' 
#' @description A generic function for getting the fieldPopInRef slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the fieldPopInRef slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fieldPopInRef="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fieldPopInRef="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fieldPopInRef(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fieldPopInRef", function(x) standardGeneric("fieldPopInRef"))


#' A getter for the fieldPopInRef slot in a RAIDSparam class
#' 
#' @description A function for getting the fieldPopInRef slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{character} string representing the name of the 
#' column that contains the known ancestry for the reference profiles in 
#' the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter).
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the fieldPopInRef slot for the object
#' fieldPopInRef(paramDemo)
#' 
#' @export
setMethod("fieldPopInRef", "RAIDSparam", function(x) {
  return(x@fieldPopInRef)
})


#' Generic function for getting the fieldPopInfAnc slot in a class
#' 
#' @description A generic function for getting the fieldPopInfAnc slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the fieldPopInfAnc slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fieldPopInfAnc="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fieldPopInfAnc="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fieldPopInfAnc(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fieldPopInfAnc", function(x) standardGeneric("fieldPopInfAnc"))


#' A getter for the fieldPopInfAnc slot in a RAIDSparam class
#' 
#' @description A function for getting the fieldPopInfAnc slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{character} string representing the name of 
#' the column in the data frame that contains the inferred 
#' super-population ancestry for the samples. The column should be 
#' present in the data frame.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the fieldPopInfAnc slot for the object
#' fieldPopInfAnc(paramDemo)
#' 
#' @export
setMethod("fieldPopInfAnc", "RAIDSparam", function(x) {
  return(x@fieldPopInfAnc)
})


#' Generic function for getting the fieldSubPop slot in a class
#' 
#' @description A generic function for getting the fieldSubPop slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the fieldSubPop slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fieldSubPop="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fieldSubPop="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fieldSubPop(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fieldSubPop", function(x) standardGeneric("fieldSubPop"))


#' A getter for the fieldSubPop slot in a RAIDSparam class
#' 
#' @description A function for getting the fieldSubPop slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{character} string representing the name of the column 
#' in the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter) that contains the sub-population 
#' information for the samples. The column should be present in the Population 
#' Reference GDS file. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the fieldSubPop slot for the object
#' fieldSubPop(paramDemo)
#' 
#' @export
setMethod("fieldSubPop", "RAIDSparam", function(x) {
  return(x@fieldSubPop)
})

#' Generic function for getting the anchorId slot in a class
#' 
#' @description A generic function for getting the anchorId slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the anchorId slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(anchorId="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", anchorId="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # anchorId(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("anchorId", function(x) standardGeneric("anchorId"))


#' A getter for the anchorId slot in a RAIDSparam class
#' 
#' @description A function for getting the anchorId slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{character} string representing the name of the
#' identifier of the quasi non-admixed profile from the reference.  
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the anchorId slot for the object
#' anchorId(paramDemo)
#' 
#' @export
setMethod("anchorId", "RAIDSparam", function(x) {
  return(x@anchorId)
})

#' Generic function for getting the pathADMIXTURE slot in a class
#' 
#' @description A generic function for getting the pathADMIXTURE slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the pathADMIXTURE slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pathADMIXTURE="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pathADMIXTURE="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pathADMIXTURE(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pathADMIXTURE", function(x) standardGeneric("pathADMIXTURE"))


#' A getter for the pathADMIXTURE slot in a RAIDSparam class
#' 
#' @description A function for getting the pathADMIXTURE slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a \code{character} string representing the path pf
#' the software ADMIXTURE and the binary.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the pathADMIXTURE slot for the object
#' pathADMIXTURE(paramDemo)
#' 
#' @export
setMethod("pathADMIXTURE", "RAIDSparam", function(x) {
  return(x@pathADMIXTURE)
})

#' Generic function for getting the seedADMIXTURE slot in a class
#' 
#' @description A generic function for getting the seedADMIXTURE slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the seedADMIXTURE slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(seedADMIXTURE="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", seedADMIXTURE="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # seedADMIXTURE(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("seedADMIXTURE", function(x) standardGeneric("seedADMIXTURE"))


#' A getter for the seedADMIXTURE slot in a RAIDSparam class
#' 
#' @description A function for getting the seedADMIXTURE slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{integer} seed for
#' ADMIXTURE program.
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the seedADMIXTURE slot for the object
#' seedADMIXTURE(paramDemo)
#' 
#' @export
setMethod("seedADMIXTURE", "RAIDSparam", function(x) {
  return(x@seedADMIXTURE)
})

#' Generic function for getting the cleanADMIXTUREFiles slot in a class
#' 
#' @description A generic function for getting the cleanADMIXTUREFiles slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the cleanADMIXTUREFiles slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(cleanADMIXTUREFiles="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", cleanADMIXTUREFiles="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # cleanADMIXTUREFiles(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("cleanADMIXTUREFiles", function(x) standardGeneric("cleanADMIXTUREFiles"))


#' A getter for the cleanADMIXTUREFiles slot in a RAIDSparam class
#' 
#' @description A function for getting the cleanADMIXTUREFiles slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{logical} indicating if message information should be 
#' printed. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the cleanADMIXTUREFiles slot for the object
#' cleanADMIXTUREFiles(paramDemo)
#' 
#' @export
setMethod("cleanADMIXTUREFiles", "RAIDSparam", function(x) {
  return(x@cleanADMIXTUREFiles)
})


#' Generic function for getting the verbose slot in a class
#' 
#' @description A generic function for getting the verbose slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @return a value from the verbose slot in the S4 object.
#' 
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(verbose="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", verbose="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # verbose(obj)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("verbose", function(x) standardGeneric("verbose"))


#' A getter for the verbose slot in a RAIDSparam class
#' 
#' @description A function for getting the verbose slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @return a single \code{logical} indicating if message information should be 
#' printed. 
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Extract the verbose slot for the object
#' verbose(paramDemo)
#' 
#' @export
setMethod("verbose", "RAIDSparam", function(x) {
  return(x@verbose)
})


###########################################################################
## All the setter functions for the RAIDSparam class
###########################################################################


#' Generic function for replacement of studyDF slot in a class
#' 
#' @description A generic function for replacement of studyDF slot in a 
#' S4 object.
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#' 
#' @examples
#'  
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(studyDF="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", studyDF="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # studyDF(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("studyDF<-", function(x, value) standardGeneric("studyDF<-"))


#' Replacement of studyDF slot in a \code{RAIDSparam} object
#' 
#' @description A function for replacement of the studyDF slot in a 
#' \code{RAIDSparam} class.
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{data.frame} containing the information 
#' about the study associated with the analyzed sample(s). The 
#' \code{data.frame} must have those 3 columns: "study.id", "study.desc", 
#' "study.platform". All columns must be character strings (no factors).
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Demo study data frame
#' studyDF <- data.frame(study.id="TEST1", study.desc="Study1",
#'     study.platform="CSHL", stringsAsFactors=FALSE)
#' 
#' ## Assign the new table to the studyDF slot in the object
#' studyDF(paramDemo) <- studyDF
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @importFrom methods validObject
#' @export
setMethod("studyDF<-", "RAIDSparam", function(x, value) {
  x@studyDF <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of studyDFSyn slot in a class
#' 
#' @description A generic function for replacement of studyDFSyn slot in a 
#' S4 object.
#' 
#' @param x an S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#' 
#' @examples
#'  
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(studyDFSyn="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", studyDFSyn="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # studyDFSyn(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("studyDFSyn<-", function(x, value) standardGeneric("studyDFSyn<-"))


#' Replacement of studyDFSyn slot in a \code{RAIDSparam} object
#' 
#' @description A function for replacement of the studyDFSyn slot in a 
#' \code{RAIDSparam} class.
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{data.frame} containing the information 
#' about the synthetic data for the analyzed sample(s). The \code{data.frame} 
#' must have those 3 columns: "study.id", "study.desc", "study.platform". All 
#' columns must be character strings (no factors).
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Demo synthetic study data frame
#' studyDFSyn <- data.frame(study.id="TEST1", study.desc="Synthetic",
#'     study.platform="CSHL", stringsAsFactors=FALSE)
#' 
#' ## Assign the new table to the studyDFSyn slot in the object
#' studyDFSyn(paramDemo) <- studyDFSyn
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @importFrom methods validObject
#' @export
setMethod("studyDFSyn<-", "RAIDSparam", function(x, value) {
  x@studyDFSyn <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of pedStudy slot in a class
#' 
#' @description A generic function for replacement of pedStudy slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pedStudy="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pedStudy="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pedStudy(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pedStudy<-", function(x, value) standardGeneric("pedStudy<-"))


#' Setter function for replacement of pedStudy slot in a RAIDSparam class
#' 
#' @description A function for replacement of pedStudy slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{data.frame} containing the information 
#' TODO with those mandatory columns: "Name.ID", "Case.ID", "Sample.Type", 
#' "Diagnosis", and "Source". All columns must be character strings 
#' (no factors). All row names should correspond to the values in the 
#' "Name.ID" column. 
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Demo PED study
#' pedStudy <- data.frame(Name.ID=c("TEST11", "Test12"),
#'     Case.ID=c("Case1", "Case2"), 
#'     Sample.Type=c("Cancer", "Cancer"), 
#'     Diagnosis=c("sarcoma", "sarcoma"),
#'     Source=c("CSHL", "CSHL"), stringsAsFactors=FALSE, 
#'     row.names=c("TEST11", "Test12"))
#' 
#' ## Assign the new PED study to the pedStudy slot in the object
#' pedStudy(paramDemo) <- pedStudy
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("pedStudy<-", "RAIDSparam", function(x, value) {
  x@pedStudy <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of studyType slot in a class
#' 
#' @description A generic function for replacement of studyType slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(studyType="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", studyType="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # studyType(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("studyType<-", function(x, value) standardGeneric("studyType<-"))


#' Setter function for replacement of studyType slot in a RAIDSparam class
#' 
#' @description A function for replacement of studyType slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string representing the type of 
#' study. The possible choices are: "LD" and "GeneAware". The study type 
#' affects how the allelic fraction is estimated.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the new study type to the studyType slot in the object
#' studyType(paramDemo) <- "GeneAware"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("studyType<-", "RAIDSparam", function(x, value) {
  x@studyType <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of genoSource slot in a class
#' 
#' @description A generic function for replacement of genoSource slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(genoSource="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", genoSource="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # genoSource(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("genoSource<-", function(x, value) standardGeneric("genoSource<-"))


#' Setter function for replacement of genoSource slot in a RAIDSparam class
#' 
#' @description A function for replacement of genoSource slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string corresponding to the type 
#' of file with the genotype and the allele information of the profile that 
#' will be provided in the 'profileFile' slot. The valid options are: "VCF", 
#' "generic", "snp-pileup", and "bam".
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the new study type to the genoSource slot in the object
#' genoSource(paramDemo) <- "bam"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("genoSource<-", "RAIDSparam", function(x, value) {
  x@genoSource <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of blockTypeId slot in a class
#' 
#' @description A generic function for replacement of blockTypeId slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(blockTypeId="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", blockTypeId="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # blockTypeId(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("blockTypeId<-", function(x, value) standardGeneric("blockTypeId<-"))


#' Setter function for replacement of blockTypeId slot in a RAIDSparam class
#' 
#' @description A function for replacement of blockTypeId slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string corresponding to 
#' the block type used to extract the block identifiers. The block type must 
#' be present in the GDS Reference Annotation file.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the new blockTypeId to the blockTypeId slot in the object
#' blockTypeId(paramDemo) <- "TEST"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("blockTypeId<-", "RAIDSparam", function(x, value) {
  x@blockTypeId <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of reference slot in a class
#' 
#' @description A generic function for replacement of reference slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(reference="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", reference="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # reference(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("reference<-", function(x, value) standardGeneric("reference<-"))


#' Setter function for replacement of reference slot in a RAIDSparam class
#' 
#' @description A function for replacement of reference slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string with two possible values: 
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the new reference to the reference slot in the object
#' reference(paramDemo) <- "1k_hgdpV0.1"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("reference<-", "RAIDSparam", function(x, value) {
  x@reference <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of chrInfo slot in a class
#' 
#' @description A generic function for replacement of chrInfo slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(chrInfo="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", chrInfo="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # chrInfo(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("chrInfo<-", function(x, value) standardGeneric("chrInfo<-"))


#' Setter function for replacement of chrInfo slot in a RAIDSparam class
#' 
#' @description A function for replacement of chrInfo slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{vector} of positive \code{integer} values 
#' representing the length of the chromosomes.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the new chrInfo to the chrInfo slot in the object
#' chrInfo(paramDemo) <- c(122222L, 333333L)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("chrInfo<-", "RAIDSparam", function(x, value) {
  x@chrInfo <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of paramAncestry slot in a class
#' 
#' @description A generic function for replacement of paramAncestry slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(paramAncestry="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", paramAncestry="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # paramAncestry(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("paramAncestry<-", function(x, value) 
    standardGeneric("paramAncestry<-"))


#' Setter function for replacement of paramAncestry slot in a RAIDSparam class
#' 
#' @description A function for replacement of paramAncestry slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{vector} of positive \code{integer} values 
#' representing the length of the chromosomes.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the new chrInfo to the paramAncestry slot in the object
#' paramAncestry(paramDemo) <- list(ScanBamParam=NULL, PileupParam=NULL, 
#'     yieldSize=1000023)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("paramAncestry<-", "RAIDSparam", function(x, value) {
  x@paramAncestry <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of profileFile slot in a class
#' 
#' @description A generic function for replacement of profileFile slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(profileFile="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", profileFile="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # profileFile(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("profileFile<-", function(x, value) 
    standardGeneric("profileFile<-"))


#' Setter function for replacement of profileFile slot in a RAIDSparam class
#' 
#' @description A function for replacement of profileFile slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value \code{NULL} or a \code{character} string representing the 
#' path to the file 
#' with genotype and the allele information of the profile. The format of 
#' the accepted file is determined by the value in the 'genoSource' slot. A 
#' profile would have a file with extension "vcf.gz" when 'genoSource' slot 
#' is "VCF". If 'genoSource' is "generic" or "snp-pileup", then ".txt.gz". If 
#' 'genoSource' is "bam", then ".bam" (the file needs to be indexed with 
#' an existing corresponding ".bai" file).
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Adjust the genoSource slot to the updated format file
#' genoSource(paramDemo) <- "generic"
#' 
#' ## Assign the existing profile file to the profileFile slot in the object
#' profileFile(paramDemo) <- file.path(dataDir, "tests/ex1.txt.gz")
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("profileFile<-", "RAIDSparam", function(x, value) {
  x@profileFile <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of pathProfileGDS slot in a class
#' 
#' @description A generic function for replacement of pathProfileGDS slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pathProfileGDS="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pathProfileGDS="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pathProfileGDS(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pathProfileGDS<-", function(x, value) 
    standardGeneric("pathProfileGDS<-"))


#' Setter function for replacement of pathProfileGDS slot in a RAIDSparam class
#' 
#' @description A function for replacement of pathProfileGDS slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value \code{NULL} or a single \code{character} string 
#' representing the path to the directory where the GDS Profile files will 
#' be created. If specified, the directory must exist.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the existing path to the pathProfileGDS slot in the object
#' pathProfileGDS(paramDemo) <- file.path(dataDir)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("pathProfileGDS<-", "RAIDSparam", function(x, value) {
  x@pathProfileGDS <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of fileReferenceGDS slot in a class
#' 
#' @description A generic function for replacement of fileReferenceGDS slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fileReferenceGDS="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fileReferenceGDS="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fileReferenceGDS(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fileReferenceGDS<-", function(x, value) 
    standardGeneric("fileReferenceGDS<-"))


#' Setter function for replacement of fileReferenceGDS slot in a RAIDSparam class
#' 
#' @description A function for replacement of fileReferenceGDS slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value \code{NULL} or a single \code{character} string representing 
#' the file name of the Reference GDS file. If specified, the file must exist.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the existing file to the fileReferenceGDS slot in the object
#' fileReferenceGDS(paramDemo) <- file.path(dataDir, 
#'     "tests/ex1_good_small_1KG.gds")
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("fileReferenceGDS<-", "RAIDSparam", function(x, value) {
  x@fileReferenceGDS <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of fileReferenceAnnotGDS slot in a class
#' 
#' @description A generic function for replacement of fileReferenceAnnotGDS 
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fileReferenceAnnotGDS="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fileReferenceAnnotGDS="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fileReferenceAnnotGDS(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fileReferenceAnnotGDS<-", function(x, value) 
    standardGeneric("fileReferenceAnnotGDS<-"))


#' Setter function for replacement of fileReferenceAnnotGDS slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of fileReferenceAnnotGDS slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value \code{NULL} or a \code{character} string representing the 
#' file name of the Population Reference GDS Annotation file.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the existing file to the fileReferenceAnnotGDS slot in the object
#' fileReferenceAnnotGDS(paramDemo) <- file.path(dataDir, 
#'     "tests/ex1_good_small_1KG_Annot.gds")
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("fileReferenceAnnotGDS<-", "RAIDSparam", function(x, value) {
  x@fileReferenceAnnotGDS <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of inferenceType slot in a class
#' 
#' @description A generic function for replacement of inferenceType 
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(inferenceType="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", inferenceType="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # inferenceType(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("inferenceType<-", function(x, value) 
    standardGeneric("inferenceType<-"))


#' Setter function for replacement of inferenceType slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of inferenceType slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string representing the genotype 
#' ancestry inference method. The two possible values are 'PCAknn' and 
#' 'haploAdmixture'.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the inferenceType slot in the object
#' inferenceType(paramDemo) <- "haploAdmixture"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("inferenceType<-", "RAIDSparam", function(x, value) {
  x@inferenceType <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of sampleRef slot in a class
#' 
#' @description A generic function for replacement of sampleRef 
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(sampleRef="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", sampleRef="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # sampleRef(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("sampleRef<-", function(x, value) 
    standardGeneric("sampleRef<-"))


#' Setter function for replacement of sampleRef slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of sampleRef slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value \code{NULL} or a \code{vector} of \code{character} strings 
#' that represent a subset of sample identifiers from the GDS Reference File. 
#' Those samples will be used for the ancestry inference. If \code{NULL}, all 
#' the samples are retained.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new values to the sampleRef slot in the object
#' sampleRef(paramDemo) <- c("Sample1", "Sample2")
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("sampleRef<-", "RAIDSparam", function(x, value) {
  x@sampleRef <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of batch slot in a class
#' 
#' @description A generic function for replacement of batch 
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(batch="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", batch="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # batch(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("batch<-", function(x, value) 
    standardGeneric("batch<-"))


#' Setter function for replacement of batch slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of batch slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single positive \code{integer} representing the current 
#' identifier for the batch. 
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the batch slot in the object
#' batch(paramDemo) <- 34L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("batch<-", "RAIDSparam", function(x, value) {
  x@batch <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of prefix slot in a class
#' 
#' @description A generic function for replacement of prefix 
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(prefix="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", prefix="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # prefix(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("prefix<-", function(x, value) 
    standardGeneric("prefix<-"))


#' Setter function for replacement of prefix slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of prefix slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string representing TODO.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the prefix slot in the object
#' prefix(paramDemo) <- "test"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("prefix<-", "RAIDSparam", function(x, value) {
  x@prefix <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of nbSim slot in a class
#' 
#' @description A generic function for replacement of nbSim 
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(nbSim="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", nbSim="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # nbSim(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("nbSim<-", function(x, value) standardGeneric("nbSim<-"))


#' Setter function for replacement of nbSim slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of nbSim slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a positive \code{integer} representing TODO.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the nbSim slot in the object
#' nbSim(paramDemo) <- 1L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("nbSim<-", "RAIDSparam", function(x, value) {
  x@nbSim <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of posOffset slot in a class
#' 
#' @description A generic function for replacement of posOffset 
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(posOffset="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", posOffset="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # posOffset(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("posOffset<-", function(x, value) standardGeneric("posOffset<-"))


#' Setter function for replacement of posOffset slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of posOffset slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{integer} that is added to the SNP position to
#' switch from 0-based to 1-based coordinate when needed (or reverse).
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the posOffset slot in the object
#' posOffset(paramDemo) <- 0L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("posOffset<-", "RAIDSparam", function(x, value) {
  x@posOffset <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of minCov slot in a class
#' 
#' @description A generic function for replacement of minCov 
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(minCov="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", minCov="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # minCov(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("minCov<-", function(x, value) standardGeneric("minCov<-"))


#' Setter function for replacement of minCov slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of minCov slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single positive \code{integer} representing the minimum 
#' required coverage.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the minCov slot in the object
#' minCov(paramDemo) <- 3L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("minCov<-", "RAIDSparam", function(x, value) {
  x@minCov <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of minProb slot in a class
#' 
#' @description A generic function for replacement of minProb  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(minProb="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", minProb="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # minProb(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("minProb<-", function(x, value) standardGeneric("minProb<-"))


#' Setter function for replacement of minProb slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of minProb slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{numeric} between \code{0} and \code{1} 
#' representing the probability that the calculated genotype call is correct.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the minProb slot in the object
#' minProb(paramDemo) <- 0.03
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("minProb<-", "RAIDSparam", function(x, value) {
  x@minProb <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})



#' Setter function for replacement of minCov slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of minCov slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single positive \code{integer} representing the minimum 
#' required coverage.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the minCov slot in the object
#' minCov(paramDemo) <- 3L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("minCov<-", "RAIDSparam", function(x, value) {
  x@minCov <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of seqError slot in a class
#' 
#' @description A generic function for replacement of seqError  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(seqError="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", seqError="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # seqError(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("seqError<-", function(x, value) standardGeneric("seqError<-"))


#' Setter function for replacement of seqError slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of seqError slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{numeric} between \code{0} and \code{1} 
#' representing the probability of sequencing error.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the seqError slot in the object
#' seqError(paramDemo) <- 0.03
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("seqError<-", "RAIDSparam", function(x, value) {
  x@seqError <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of seqErrorSyn slot in a class
#' 
#' @description A generic function for replacement of seqErrorSyn  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(seqErrorSyn="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", seqErrorSyn="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # seqErrorSyn(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("seqErrorSyn<-", function(x, value) standardGeneric("seqErrorSyn<-"))


#' Setter function for replacement of seqErrorSyn slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of seqErrorSyn slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{numeric} between \code{0} and \code{1} 
#' representing the probability of sequencing error for synthetic.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the seqErrorSyn slot in the object
#' seqErrorSyn(paramDemo) <- 0.03
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("seqErrorSyn<-", "RAIDSparam", function(x, value) {
  x@seqErrorSyn <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of pRecomb slot in a class
#' 
#' @description A generic function for replacement of pRecomb  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pRecomb="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pRecomb="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pRecomb(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pRecomb<-", function(x, value) standardGeneric("pRecomb<-"))


#' Setter function for replacement of pRecomb slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of pRecomb slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single positive \code{numeric} between \code{0} and 
#' \code{1} that represents the frequency of phase switching in the 
#' synthetic profiles.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the pRecomb slot in the object
#' pRecomb(paramDemo) <- 0.03
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("pRecomb<-", "RAIDSparam", function(x, value) {
  x@pRecomb <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of np slot in a class
#' 
#' @description A generic function for replacement of np  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(np="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", np="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # np(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("np<-", function(x, value) standardGeneric("np<-"))


#' Setter function for replacement of np slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of np slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single positive \code{integer} specifying the number of 
#' threads to be used.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign a new value to the np slot in the object
#' np(paramDemo) <-3L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("np<-", "RAIDSparam", function(x, value) {
  x@np <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of syntheticRefDF slot in a class
#' 
#' @description A generic function for replacement of syntheticRefDF  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(syntheticRefDF="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", syntheticRefDF="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # syntheticRefDF(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("syntheticRefDF<-", 
    function(x, value) standardGeneric("syntheticRefDF<-"))


#' Setter function for replacement of syntheticRefDF slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of syntheticRefDF slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value \code{NULL} or a \code{data.frame} containing a subset of 
#' reference profiles for each sub-population present in the Reference GDS 
#' file. The \code{data.frame} must have those columns: 
#' \describe{ 
#' \item{sample.id}{ a \code{character} string representing the sample 
#' identifier. } 
#' \item{pop.group}{ a \code{character} string representing the 
#' subcontinental population assigned to the sample. } 
#' \item{superPop}{ a \code{character} string representing the 
#' super-population assigned to the sample. }
#' }
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new values to the syntheticRefDF slot in the object
#' syntheticRefDF(paramDemo) <- data.frame(sample.id=c("A", "B"), 
#'     pop.group=c("CEU", "CEU"), superPop=c("EUR", "EUR"))
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("syntheticRefDF<-", "RAIDSparam", function(x, value) {
  x@syntheticRefDF <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of listPos slot in a class
#' 
#' @description A generic function for replacement of listPos  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(listPos="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", listPos="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # listPos(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("listPos<-", function(x, value) standardGeneric("listPos<-"))


#' Setter function for replacement of listPos slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of listPos slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value \code{NULL} or a \code{data.frame} containing 2 columns named: 
#' "snp.chromosome" and "snp.position". The first column, called 
#' "snp.chromosome", contains the name of the chromosome where the SNV 
#' is located. The second column, called "snp.position", contains the 
#' position of the SNV on the chromosome.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new values to the listPos slot in the object
#' listPos(paramDemo) <-data.frame(snp.chromosome=c("1", "1"), 
#'     snp.position=c(12232, 433333))
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("listPos<-", "RAIDSparam", function(x, value) {
  x@listPos <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of pruningMethod slot in a class
#' 
#' @description A generic function for replacement of pruningMethod  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pruningMethod="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pruningMethod="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pruningMethod(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pruningMethod<-", 
  function(x, value) standardGeneric("pruningMethod<-"))


#' Setter function for replacement of pruningMethod slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of pruningMethod slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string representing the method that 
#' will be used to calculate the linkage disequilibrium in the 
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values 
#' are: "corr", "r", "dprime", and "composite".
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the pruningMethod slot in the object
#' pruningMethod(paramDemo) <- "r"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("pruningMethod<-", "RAIDSparam", function(x, value) {
  x@pruningMethod <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of slideWindowMaxBP slot in a class
#' 
#' @description A generic function for replacement of slideWindowMaxBP  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(slideWindowMaxBP="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", slideWindowMaxBP="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # slideWindowMaxBP(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("slideWindowMaxBP<-", 
  function(x, value) standardGeneric("slideWindowMaxBP<-"))


#' Setter function for replacement of slideWindowMaxBP slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of slideWindowMaxBP slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single positive \code{integer} that represents the maximum 
#' basepairs (bp) in the sliding window. This parameter is used 
#' for the LD pruning done by the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the slideWindowMaxBP slot in the object
#' slideWindowMaxBP(paramDemo) <- 1000L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("slideWindowMaxBP<-", "RAIDSparam", function(x, value) {
  x@slideWindowMaxBP <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of thresholdLD slot in a class
#' 
#' @description A generic function for replacement of thresholdLD  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(thresholdLD="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", thresholdLD="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # thresholdLD(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("thresholdLD<-", 
  function(x, value) standardGeneric("thresholdLD<-"))


#' Setter function for replacement of thresholdLD slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of thresholdLD slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single positive \code{numeric} value that represents 
#' the LD threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. 
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the thresholdLD slot in the object
#' thresholdLD(paramDemo) <- 1000L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("thresholdLD<-", "RAIDSparam", function(x, value) {
  x@thresholdLD <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of specificSNV slot in a class
#' 
#' @description A generic function for replacement of specificSNV  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(specificSNV="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", specificSNV="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # specificSNV(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("specificSNV<-", 
  function(x, value) standardGeneric("specificSNV<-"))


#' Setter function for replacement of specificSNV slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of specificSNV slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value \code{NULL} or a \code{data.frame} containing 2 columns. The 
#' first column, called "snp.chromosome" contains the name of the chromosome 
#' wherethe SNV is located. The second column, called "snp.position", contains 
#' the position of the SNV on the chromosome. Optionally, the column "snvKeep" 
#' contains the SNV index in the reference or -1 if not in the reference.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the specificSNV slot in the object
#' specificSNV(paramDemo) <- data.frame(snp.chromosome=c("1", "1"), 
#'     snp.position=c(13333, 2222213))
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("specificSNV<-", "RAIDSparam", function(x, value) {
  x@specificSNV <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of genoType slot in a class
#' 
#' @description A generic function for replacement of genoType  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(genoType="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", genoType="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # genoType(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("genoType<-", 
  function(x, value) standardGeneric("genoType<-"))


#' Setter function for replacement of genoType slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of genoType slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string representing the field 
#' in the GDS profile file that contains all the genotype information. 
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the genoType slot in the object
#' genoType(paramDemo) <- "geno.REF"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("genoType<-", "RAIDSparam", function(x, value) {
  x@genoType <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of phaseType slot in a class
#' 
#' @description A generic function for replacement of phaseType  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(phaseType="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", phaseType="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # phaseType(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("phaseType<-", function(x, value) standardGeneric("phaseType<-"))


#' Setter function for replacement of phaseType slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of phaseType slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{character} string TODO
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the phaseType slot in the object
#' phaseType(paramDemo) <- "test"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("phaseType<-", "RAIDSparam", function(x, value) {
  x@phaseType <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})

#' Generic function for replacement of phase slot in a class
#' 
#' @description A generic function for replacement of phase  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(phase="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", phase="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # phase(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("phase<-", function(x, value) standardGeneric("phase<-"))


#' Setter function for replacement of phase slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of phase slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{logical} TODO
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the phase slot in the object
#' phase(paramDemo) <-  TRUE
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("phase<-", "RAIDSparam", function(x, value) {
  x@phase <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of PCAmissingRate slot in a class
#' 
#' @description A generic function for replacement of PCAmissingRate  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(PCAmissingRate="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", PCAmissingRate="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # PCAmissingRate(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("PCAmissingRate<-", 
  function(x, value) standardGeneric("PCAmissingRate<-"))


#' Setter function for replacement of PCAmissingRate slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of PCAmissingRate slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a positive \code{numeric} representing the maximum missing 
#' rate accepted to use SNPs in the PCA analysis done with 
#' the \link[SNPRelate]{snpgdsPCA} function. If \code{NaN}, no 
#' missing threshold.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the PCAmissingRate slot in the object
#' PCAmissingRate(paramDemo) <-  0.03322
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("PCAmissingRate<-", "RAIDSparam", function(x, value) {
  x@PCAmissingRate <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of PCAalgorithm slot in a class
#' 
#' @description A generic function for replacement of PCAalgorithm  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(PCAalgorithm="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", PCAalgorithm="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # PCAalgorithm(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("PCAalgorithm<-", 
  function(x, value) standardGeneric("PCAalgorithm<-"))


#' Setter function for replacement of PCAalgorithm slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of PCAalgorithm slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{character} string representing the algorithm to 
#' use with the \link[SNPRelate]{snpgdsPCA} function. The algorithm must be 
#' implemented and available to the \link[SNPRelate]{snpgdsPCA} function.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the PCAalgorithm slot in the object
#' PCAalgorithm(paramDemo) <-  "randomized"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("PCAalgorithm<-", "RAIDSparam", function(x, value) {
  x@PCAalgorithm <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of eigenCount slot in a class
#' 
#' @description A generic function for replacement of eigenCount  
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(eigenCount="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", eigenCount="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # eigenCount(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("eigenCount<-", 
  function(x, value) standardGeneric("eigenCount<-"))


#' Setter function for replacement of eigenCount slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of eigenCount slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{integer} indicating the number of 
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA} 
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the eigenCount slot in the object
#' eigenCount(paramDemo) <-  33L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("eigenCount<-", "RAIDSparam", function(x, value) {
  x@eigenCount <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of eigenCountSyn slot in a class
#' 
#' @description A generic function for replacement of eigenCountSyn   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(eigenCountSyn="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", eigenCountSyn="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # eigenCountSyn(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("eigenCountSyn<-", 
  function(x, value) standardGeneric("eigenCountSyn<-"))


#' Setter function for replacement of eigenCountSyn slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of eigenCountSyn slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{integer} indicating the number of 
#' eigenvectors TODO
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the eigenCountSyn slot in the object
#' eigenCountSyn(paramDemo) <-  33L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("eigenCountSyn<-", "RAIDSparam", function(x, value) {
  x@eigenCountSyn <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of kList slot in a class
#' 
#' @description A generic function for replacement of kList   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(kList="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", kList="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # kList(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("kList<-", 
  function(x, value) standardGeneric("kList<-"))


#' Setter function for replacement of kList slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of kList slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{vector} of positive \code{integer} representing the 
#' values tested for the  \emph{K} parameter. The \emph{K} parameter 
#' represents the number of neighbors used in the K-nearest neighbor analysis.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the kList slot in the object
#' kList(paramDemo) <-  c(2L, 3L, 4L, 5L)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("kList<-", "RAIDSparam", function(x, value) {
  x@kList <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of pcaList slot in a class
#' 
#' @description A generic function for replacement of pcaList   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pcaList="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pcaList="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pcaList(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pcaList<-", function(x, value) standardGeneric("pcaList<-"))


#' Setter function for replacement of pcaList slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of pcaList slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{vector} of positive \code{integer} representing the 
#' values tested for the  \emph{D} parameter. The \emph{D} parameter 
#' represents the number of dimensions used in the PCA analysis.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the pcaList slot in the object
#' pcaList(paramDemo) <-  c(2L, 3L, 4L, 5L)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("pcaList<-", "RAIDSparam", function(x, value) {
  x@pcaList <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of fieldPopInRef slot in a class
#' 
#' @description A generic function for replacement of fieldPopInRef   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fieldPopInRef="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fieldPopInRef="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fieldPopInRef(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fieldPopInRef<-", 
  function(x, value) standardGeneric("fieldPopInRef<-"))


#' Setter function for replacement of fieldPopInRef slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of fieldPopInRef slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{character} string representing the name of the 
#' column that contains the known ancestry for the reference profiles in 
#' the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter).
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the fieldPopInRef slot in the object
#' fieldPopInRef(paramDemo) <-  "newPop"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("fieldPopInRef<-", "RAIDSparam", function(x, value) {
  x@fieldPopInRef <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of fieldPopInfAnc slot in a class
#' 
#' @description A generic function for replacement of fieldPopInfAnc   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fieldPopInfAnc="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fieldPopInfAnc="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fieldPopInfAnc(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fieldPopInfAnc<-", 
  function(x, value) standardGeneric("fieldPopInfAnc<-"))


#' Setter function for replacement of fieldPopInfAnc slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of fieldPopInfAnc slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{character} string representing the name of 
#' the column in the data frame that contains the inferred 
#' super-population ancestry for the samples. The column should be 
#' present in the data frame. 
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the fieldPopInfAnc slot in the object
#' fieldPopInfAnc(paramDemo) <-  "superPOPULATION"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("fieldPopInfAnc<-", "RAIDSparam", function(x, value) {
  x@fieldPopInfAnc <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of fieldSubPop slot in a class
#' 
#' @description A generic function for replacement of fieldSubPop   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(fieldSubPop="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", fieldSubPop="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # fieldSubPop(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("fieldSubPop<-", 
  function(x, value) standardGeneric("fieldSubPop<-"))


#' Setter function for replacement of fieldSubPop slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of fieldSubPop slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{character} string representing the name of the 
#' column in the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter) that contains the sub-population 
#' information for the samples. The column should be present in the Population 
#' Reference GDS file.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the fieldSubPop slot in the object
#' fieldSubPop(paramDemo) <-  "super.POPULATION.group"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("fieldSubPop<-", "RAIDSparam", function(x, value) {
  x@fieldSubPop <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})

#' Generic function for replacement of anchorId slot in a class
#' 
#' @description A generic function for replacement of anchorId   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(anchorId="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", anchorId="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # anchorId(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("anchorId<-", 
  function(x, value) standardGeneric("anchorId<-"))


#' Setter function for replacement of anchorId slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of anchorId slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{character} string representing the name of the
#' identifier of the quasi non-admixed profile from the reference.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the anchorId slot in the object
#' anchorId(paramDemo) <-  "sp0.95vT"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("anchorId<-", "RAIDSparam", function(x, value) {
  x@anchorId <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})

#' Generic function for replacement of pathADMIXTURE slot in a class
#' 
#' @description A generic function for replacement of pathADMIXTURE slot in a 
#' S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(pathADMIXTURE="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", pathADMIXTURE="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # pathADMIXTURE(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("pathADMIXTURE<-", function(x, value) 
    standardGeneric("pathADMIXTURE<-"))


#' Setter function for replacement of pathADMIXTURE slot in a RAIDSparam class
#' 
#' @description A function for replacement of pathADMIXTURE slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a \code{character} string representing the path pf
#' the software ADMIXTURE and the binary. If specified, the 
#' file must exist.
#
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Directory where demo GDS files are located
#' dataDir <- system.file("extdata", package="RAIDS")
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign the existing path to the pathADMIXTURE slot in the object
#' pathADMIXTURE(paramDemo) <- file.path(dataDir)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("pathADMIXTURE<-", "RAIDSparam", function(x, value) {
  x@pathADMIXTURE <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})

#' Generic function for replacement of seedADMIXTURE slot in a class
#' 
#' @description A generic function for replacement of seedADMIXTURE   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(seedADMIXTURE="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", seedADMIXTURE="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # seedADMIXTURE(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("seedADMIXTURE<-", 
  function(x, value) standardGeneric("seedADMIXTURE<-"))


#' Setter function for replacement of seedADMIXTURE slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of seedADMIXTURE slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{integer} seed for
#' the ADMIXTURE program.
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the seedADMIXTURE slot in the object
#' seedADMIXTURE(paramDemo) <-  33L
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("seedADMIXTURE<-", "RAIDSparam", function(x, value) {
  x@seedADMIXTURE <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})

#' Generic function for replacement of cleanADMIXTUREFiles slot in a class
#' 
#' @description A generic function for replacement of cleanADMIXTUREFiles   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(cleanADMIXTUREFiles="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", cleanADMIXTUREFiles="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # cleanADMIXTUREFiles(obj) <- TRUE
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("cleanADMIXTUREFiles<-", 
  function(x, value) standardGeneric("cleanADMIXTUREFiles<-"))


#' Setter function for replacement of cleanADMIXTUREFiles slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of cleanADMIXTUREFiles slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{logical} indicating if message information 
#' should be printed. 
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the cleanADMIXTUREFiles slot in the object
#' cleanADMIXTUREFiles(paramDemo) <-  TRUE
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("cleanADMIXTUREFiles<-", "RAIDSparam", function(x, value) {
  x@cleanADMIXTUREFiles <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


#' Generic function for replacement of verbose slot in a class
#' 
#' @description A generic function for replacement of verbose   
#' slot in a S4 object. 
#' 
#' @param x a S4 object.
#' 
#' @param value the new value to assign or update.
#' 
#' @return the modified S4 object when the new value is valid.
#'  
#' @examples
#' 
#' # Define a dummy class to show usage
#' setClass("MyClass", slots = list(verbose="character"))
#' 
#' # Create an instance
#' obj <- new("MyClass", verbose="123")
#' 
#' # Call the generic (assuming a method is implemented)
#' # verbose(obj) <- "333"
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @export
setGeneric("verbose<-", 
  function(x, value) standardGeneric("verbose<-"))


#' Setter function for replacement of verbose slot in 
#' a RAIDSparam class
#' 
#' @description A function for replacement of verbose slot in a 
#' \code{RAIDSparam} class. 
#' 
#' @param x a \code{RAIDSparam} object.
#' 
#' @param value a single \code{logical} indicating if message information 
#' should be printed. 
#' 
#' @return the modified \code{RAIDSparam} object when the new value is valid.
#' 
#' @examples
#' 
#' ## Create a RAIDSparam object
#' paramDemo <- RAIDSparam()
#' 
#' ## Assign new value to the verbose slot in the object
#' verbose(paramDemo) <-  TRUE
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' 
#' @importFrom methods validObject
#' @export
setMethod("verbose<-", "RAIDSparam", function(x, value) {
  x@verbose <- value

  # Validate and return the modified object
  validObject(x) 
  return(x)
})


###########################################################################
## RAIDSparam function
###########################################################################

#' @title Create a RAIDSparam object 
#'
#' @description The function uses all parameters to create a 
#' \code{RAIDSparam} class object that contains all the parameters used 
#' in the RAIDS workflow.
#' 
#' @param studyDF a \code{data.frame} containing the information about the
#' study associated to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be \code{character} strings (no factors). If \code{NULL}, the 
#' following will be assigned: 
#' \code{data.frame(study.id="NotDef", study.desc="NotDef", 
#' study.platform="NotDef", stringsAsFactors=FALSE)}.
#' Default: \code{NULL}.
#' 
#' @param studyDFSyn a \code{data.frame} containing the information about the
#' synthetic data to the analysed sample(s). The \code{data.frame} must have
#' those 3 columns: "study.id", "study.desc", "study.platform". All columns
#' must be \code{character} strings (no factors). If \code{NULL}, the 
#' following will be assigned: \code{data.frame(study.id="NotDef.Synthetic", 
#' study.desc="NotDef synthetic data", study.platform="Synthetic", 
#' stringsAsFactors=FALSE)}.
#' Default: \code{NULL}.
#' 
#' @param pedStudy a \code{data.frame} containing the information TODO 
#' with those mandatory columns: "Name.ID",
#' "Case.ID", "Sample.Type", "Diagnosis", and "Source". All columns must be in
#' \code{character} strings (no factors). All row names should correspond to 
#' the values in the "Name.ID" column. If \code{NULL}, the 
#' following will be assigned:  \code{data.frame(Name.ID=c("ProfileId"), 
#' Case.ID=c("ProfileId"), Sample.Type=c("type"), Diagnosis="NotDef", 
#' Source=c("NotDef"), stringsAsFactors=FALSE, row.names = c("ProfileId"))}.
#' Default: \code{NULL}.
#' 
#' @param studyType a single \code{character} string representing the type 
#' of study. The possible choices are: "LD" and "GeneAware". The study type 
#' affects how the allelic fraction is estimated. 
#' Default: \code{"LD"}.
#' 
#' @param genoSource a single \code{character} string corresponding to the type 
#' of file with the genotype and the allele information of the profile that 
#' will be provided in the 'profileFile' slot. The valid options are: "VCF", 
#' "generic", "snp-pileup", and "bam". The "generic" format CSV file
#' must have at least these columns: 'Chromosome', 'Position', 'Ref', 'Alt', 
#' 'Count', 'File1R', and 'File1A'. The 'Count' is the depth at the 
#' specified position; 'FileR' is the depth of the reference allele, and
#' 'File1A' is the depth of the specific alternative allele.
#' Finally, in the case of a "VCF" file, the file must have at least those 
#' genotype fields: GT, AD, and DP. Default: \code{NULL}.
#' 
#' @param blockTypeId a single \code{character} string corresponding to 
#' the block type used to extract the block identifiers. The block type must 
#' be present in the GDS Reference Annotation file. 
#' Default: \code{"GeneS.Ensembl.Hsapiens.v86"}.
#' 
#' @param reference a \code{character} string with two possible values:
#' '1KGv1.0', '1k_hgdpV0.1'. It specifies the type of inference. 
#' Default: \code{"1KGv1.0"}.
#' 
#' @param genome a \code{character} string with one possible value:
#' 'HG38'. It specifies the genome used. Default: \code{"HG38"}.
#' 
#' @param chrInfo a \code{vector} of positive \code{integer} values
#' representing the length of the chromosomes. See 'details' section. 
#' If \code{NULL}, the following value is 
#' assigned: \code{seqlengths(Hsapiens)[seq_len(25)]}. Default: \code{NULL}.
#' 
#' @param paramAncestry a \code{list} of parameters related to ancestry. 
#' The \code{list} should contain those three entries: \code{ScanBamParam}, 
#' \code{PileupParam}, and \code{yieldSize}. If \code{NULL}, the default value 
#' \code{list(ScanBamParam=NULL, PileupParam=NULL, yieldSize=10000000)} will 
#' be assigned to the parameter. Default: \code{NULL}.
#' 
#' @param profileFile a \code{character} string representing the path to the
#' file with genotype and the allele information of the profile. The format 
#' of the accepted file is determined by the value in the 'genoSource' slot.
#' A profile would have a file with extension "*vcf.gz" when 'genoSource' slot 
#' is "VCF". If 'genoSource' is "generic" or "snp-pileup", then "*.txt.gz". 
#' If 'genoSource' is "bam", then "*.bam" (the file need to be indexed with a
#' existing corresponding "*.bai" file). Default: \code{NULL}.
#' 
#' @param profileFileGeno TODO
#' 
#' @param pathProfileGDS  a single \code{character} string 
#' representing the path to the directory where the GDS Profile files will 
#' be created. If specified, the directory must exist. 
#' Default: \code{NULL}.
#' 
#' @param fileReferenceGDS a single \code{character} string representing 
#' the file name of the Reference GDS file. If specified, the file  
#' must exist. Default: \code{NULL}.
#' 
#' @param fileReferenceAnnotGDS a \code{character} string representing the
#' file name of the Population Reference GDS Annotation file. If specified,
#' the file must exist. Default: \code{NULL}.
#' 
#' @param inferenceType a single \code{character} string representing the 
#' genotype ancestry inference method. The two possible values:
#' 'PCAknn' and 'haploAdmixture'. Default: \code{"PCAknn"}.
#' 
#' @param sampleRef a \code{vector} of \code{character} strings that represent 
#' a subset of sample identifiers from the GDS Reference File. Those samples 
#' will be used for the ancestry inference. If \code{NULL}, 
#' all the samples are retained. This is usefull if you want to removed some 
#' specific samples (ex: samples with admixed ancestry).
#' Default: \code{NULL}.
#' 
#' @param batch a single positive \code{integer} representing the current
#' identifier for the batch. Beware, this field is not stored anymore. 
#' Default: \code{1L}.
#' 
#' @param prefix a single \code{character} string representing TODO.
#' Default: \code{"1"}.
#' 
#' @param nbSim a positive \code{integer} representing the number  
#' of replicated simulation (Only 1 is supported for the 
#' moment). Default: \code{1L}.
#' 
#' @param posOffset a single \code{integer} that is added to the SNP position 
#' to switch from 0-based to 1-based coordinate when needed (or reverse).
#' Default: \code{-1L}.
#' 
#' @param minCov a single positive \code{integer} representing the minimum
#' required coverage. Default: \code{10L}.
#' 
#' @param minProb a single \code{numeric} between \code{0} and \code{1}
#' representing the probability that the calculated genotype call is correct.
#' Default: \code{0.999}.
#' 
#' @param seqError a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error. 
#' Default: \code{0.001}.
#' 
#' @param seqErrorSyn a single \code{numeric} between \code{0} and \code{1}
#' representing the probability of sequencing error for synthetic.
#' Default: \code{0.001}.
#' 
#' @param pRecomb a single positive \code{numeric} between \code{0} and 
#' \code{1} that represents the frequency of phase switching in the 
#' synthetic profiles. Default: \code{0.01}.
#' 
#' @param np a single positive \code{integer} specifying the number of
#' threads to be used. Default: \code{1L}.
#' 
#' @param listPos a \code{data.frame} containing 2 columns named: 
#' "snp.chromosome" and "snp.position". The first column,
#' called "snp.chromosome", contains the name of the chromosome where the
#' SNV is located. The second column, called "snp.position", contains the
#' position of the SNV on the chromosome. Default: \code{NULL}.
#' 
#' @param syntheticRefDF a \code{data.frame} containing a subset of
#' reference profiles for each sub-population present in the Reference GDS
#' file. The \code{data.frame} must have those columns:
#' \describe{
#' \item{sample.id}{ a \code{character} string representing the sample
#' identifier. }
#' \item{pop.group}{ a \code{character} string representing the
#' subcontinental population assigned to the sample. }
#' \item{superPop}{ a \code{character} string representing the
#' super-population assigned to the sample. }
#' }
#' Default: \code{NULL}.
#' 
#' @param pruningMethod a \code{character} string representing the method that 
#' will be used to calculate the linkage disequilibrium in the
#' \code{\link[SNPRelate]{snpgdsLDpruning}}() function. The 4 possible values
#' are: "corr", "r", "dprime", and "composite". Default: \code{"corr"}.
#' 
#' @param slideWindowMaxBP a single positive \code{integer} that represents
#' the maximum basepairs (bp) in the sliding window. This parameter is used
#' for the LD pruning done by the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. Default: \code{500000L}.
#' 
#' @param thresholdLD a single positive \code{numeric} value that represents 
#' the LD threshold used in the \code{\link[SNPRelate]{snpgdsLDpruning}} 
#' function. Default: \code{sqrt(0.1)}.
#' 
#' @param specificSNV a \code{data.frame} containing 2 columns. The first 
#' column, called "snp.chromosome" contains the name of the chromosome where 
#' the SNV is located. The second column, called "snp.position", contains the
#' position of the SNV on the chromosome. Optionally, the column "snvKeep" 
#' contains the SNV index in the reference or -1 if not in the reference. It 
#' is used during the pruning step. Default: \code{NULL}.
#' 
#' @param genoType a single \code{character} string representing the field 
#' in the GDS profile file that contains all the genotype information. 
#' Default: \code{"geno.ref"}.
#' 
#' @param phaseType a TODO. Default: \code{"phase.ref"}.
#' 
#' @param phase a \code{logical} indicating TODO. Default: \code{FALSE}.
#' 
#' @param PCAmissingRate a positive \code{numeric} representing the maximum 
#' missing rate retained accepted to use SNPs in the PCA analysis done 
#' with the \link[SNPRelate]{snpgdsPCA} function. If \code{NaN}, no 
#' missing threshold. Default: \code{0.025}.
#' 
#' @param PCAalgorithm a \code{character} string representing the algorithm 
#' to use with the \link[SNPRelate]{snpgdsPCA} function. The algorithm must 
#' be implemented and available to the \link[SNPRelate]{snpgdsPCA} function.
#' Default: \code{"exact"}.
#' 
#' @param eigenCount a single \code{integer} indicating the number of
#' eigenvectors that will be in the output of the \link[SNPRelate]{snpgdsPCA}
#' function; if 'eigenCount' <= 0, then all eigenvectors are returned. 
#' Default: \code{32L}.
#' 
#' @param eigenCountSyn a TODO. Default: \code{15L}.
#' 
#' @param kList a \code{vector} of positive \code{integer} representing the 
#' values tested for the  _K_ parameter. The _K_ parameter represents the
#' number of neighbors used in the K-nearest neighbor analysis. 
#' Default: \code{seq(2L, 15L, 1L)}.
#' 
#' @param pcaList a \code{vector} of positive \code{integer} representing 
#' the values tested for the  _D_ parameter. The _D_ parameter represents the
#' number of dimensions used in the PCA analysis. 
#' Default: \code{seq(2L, 15L, 1L)}.
#' 
#' @param fieldPopInRef a \code{character} string representing the name of the
#' column that contains the known ancestry for the reference profiles in
#' the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter). Default: \code{"superPop"}.
#' 
#' @param fieldPopInfAnc a \code{character} string representing the name of 
#' the column in the data frame that contains the inferred 
#' super-population ancestry for the samples. The column should be 
#' present in the data frame. Default: \code{"superPop"}.
#' 
#' @param fieldSubPop a \code{character} string representing the name of 
#' the column in the Population Reference GDS file (corresponding to the 
#' \code{fileReferenceGDS} parameter) that contains the sub-population 
#' information for the samples. The column should be present in the Population 
#' Reference GDS file. Default: \code{"pop.group"}.
#' 
#' @param anchorId a \code{character} string representing the name of the
#' identifier of the quasi non-admixed profile from the reference.
#' Default: \code{"sp0.95v0.2"}
#' 
#' @param pathADMIXTURE a \code{character} string representing the path pf
#' the software ADMIXTURE and the binary. If specified, the file must exist. 
#' Default: \code{NULL}.
#' 
#' @param seedADMIXTURE a single \code{integer} seed for 
#' the ADMIXTURE program. Default: \code{432L}.
#' 
#' @param cleanADMIXTUREFiles a single \code{logical} indicating if RAIDS
#' removes ADMIXTURE files after processing. Default: \code{TRUE}. 
#' 
#' @param verbose a \code{logical} indicating if message information should be
#' printed. Default: \code{FALSE}
#' 
#' @return an object of class \code{RAIDSparam} that contains all the 
#' required parameters needed by the RAIDS workflow.
#' 
#' @details
#'
#' The `chrInfo` parameter contains the length of the chromosomes. The
#' length of the chromosomes can be obtain through the
#' \code{\link[Seqinfo]{seqlengths}} library.
#'
#' As example, for hg38 genome:
#'
#' ```
#'
#' if (requireNamespace("Seqinfo", quietly=TRUE) &&
#'      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly=TRUE)) {
#'      chrInfo <- Seqinfo::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:25]
#' }
#'
#' ```
#' 
#' @examples
#'
#' ## New object of class "RAIDSparam" with default parameters
#' newParam1 <- RAIDSparam()
#' 
#' ## New object of class "RAIDSparam" with non-default parameters
#' newParam2 <- RAIDSparam(studyDF=data.frame(study.id="Test1",
#'     study.desc="Test Desc", study.platform="Test", stringsAsFactors=FALSE),
#'     nbSim=1L, minCov=8L)
#' 
#' @author Pascal Belleau, Astrid Deschênes and Alexander Krasnitz
#' @encoding UTF-8
#' @importFrom rlang arg_match
#' @importFrom methods new
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 Hsapiens
#' @importFrom Seqinfo seqlengths
#' @export
RAIDSparam <- function(studyDF=NULL, studyDFSyn=NULL, pedStudy=NULL, 
    studyType=c("LD", "GeneAware"), genoSource=NULL, 
    blockTypeId="GeneS.Ensembl.Hsapiens.v86", reference="1KGv1.0", 
    genome="HG38", chrInfo=NULL, paramAncestry=NULL, profileFile=NULL,
    profileFileGeno=NULL, pathProfileGDS=NULL, fileReferenceGDS=NULL,
    fileReferenceAnnotGDS=NULL, inferenceType="PCAknn", sampleRef=NULL,
    batch=1L, prefix="1", nbSim=1L, posOffset=-1L, minCov=10L, minProb=0.999,
    seqError=0.001, seqErrorSyn=0.001, pRecomb=0.01, np=1L, listPos=NULL, 
    syntheticRefDF=NULL, pruningMethod=c("corr", "r", "dprime", "composite"), 
    slideWindowMaxBP=500000L, thresholdLD=sqrt(0.1), specificSNV=NULL, 
    genoType="geno.ref", phaseType="phase.ref", phase=FALSE, 
    PCAmissingRate=0.025, PCAalgorithm="exact", eigenCount=32L, 
    eigenCountSyn=15L, kList=seq(2L, 15L, 1L), pcaList=seq(2L, 15L, 1L), 
    fieldPopInRef="superPop", fieldPopInfAnc="superPop", 
    fieldSubPop="pop.group", anchorId="sp0.95v0.2",
    pathADMIXTURE=NULL, seedADMIXTURE=432L, cleanADMIXTUREFiles=TRUE, verbose=FALSE) {
    
    batch <- as.integer(batch)
    nbSim <- as.integer(nbSim)
    posOffset <- as.integer(posOffset)
    minCov <- as.integer(minCov)
    np <- as.integer(np)
    slideWindowMaxBP <- as.integer(slideWindowMaxBP)
    eigenCount <- as.integer(eigenCount)
    eigenCountSyn <- as.integer(eigenCountSyn)
    kList <- unique(as.integer(kList))
    pcaList <- unique(as.integer(pcaList))
    seedADMIXTURE <- as.integer(seedADMIXTURE)
    cleanADMIXTUREFiles <- as.logical(cleanADMIXTUREFiles)
    verbose <- as.logical(verbose)

    studyType <- arg_match(studyType)
    pruningMethod <- arg_match(pruningMethod)
    
    if (is.null(studyDF)) {
        studyDF <- data.frame(study.id="NotDef",
                                study.desc="NotDef",
                                study.platform="NotDef",
                                stringsAsFactors=FALSE)
    }

    if (is.null(studyDFSyn)) {
        studyDFSyn <- data.frame(study.id=paste0(studyDF$study.id, 
            ".Synthetic"),
            study.desc=paste0(studyDF$study.id, " synthetic data"),
            study.platform="Synthetic", stringsAsFactors=FALSE)
    }

    if (is.null(pedStudy)) {
        pedStudy <- data.frame(Name.ID=c("ProfileId"),
                                Case.ID=c("ProfileId"),
                                Sample.Type=c("type"),
                                Diagnosis="NotDef",
                                Source=c("NotDef"),
                                stringsAsFactors=FALSE,
                                row.names=c("ProfileId"))
    }
  
    ## Assign default ancestry parameters when not assigned by user
    if(is.null(paramAncestry)) {
      paramAncestry <- list(ScanBamParam=NULL, PileupParam=NULL, 
                                yieldSize=10000000)
    }

    ## Assigne default chromosome information whe not assigned by user
    if (is.null(chrInfo)) {
        chrInfo <- seqlengths(Hsapiens)[seq_len(25)]
    }

    new("RAIDSparam", studyDF=studyDF, studyDFSyn=studyDFSyn, pedStudy=pedStudy, 
    studyType=studyType, genoSource=genoSource, blockTypeId=blockTypeId,
    reference=reference, genome=genome, chrInfo=chrInfo, 
    paramAncestry=paramAncestry, profileFile=profileFile, 
    profileFileGeno=profileFileGeno,
    pathProfileGDS=pathProfileGDS, fileReferenceGDS=fileReferenceGDS,
    fileReferenceAnnotGDS=fileReferenceAnnotGDS,
    inferenceType=inferenceType, sampleRef=sampleRef, batch=batch, 
    prefix=prefix, nbSim=nbSim, posOffset=posOffset, minCov=minCov, 
    minProb=minProb, seqError=seqError, seqErrorSyn=seqErrorSyn, 
    pRecomb=pRecomb, np=np, listPos=listPos, syntheticRefDF=syntheticRefDF, 
    pruningMethod=pruningMethod, slideWindowMaxBP=slideWindowMaxBP, 
    thresholdLD=thresholdLD,
    specificSNV=specificSNV, genoType=genoType, phaseType=phaseType, 
    phase=phase, PCAmissingRate=PCAmissingRate, PCAalgorithm=PCAalgorithm,
    eigenCount=eigenCount, eigenCountSyn=eigenCountSyn, kList=kList,
    pcaList=pcaList, fieldPopInRef=fieldPopInRef, fieldPopInfAnc=fieldPopInfAnc,
    fieldSubPop=fieldSubPop, anchorId=anchorId,
    pathADMIXTURE=pathADMIXTURE, seedADMIXTURE=seedADMIXTURE, 
    cleanADMIXTUREFiles=cleanADMIXTUREFiles, 
    verbose=verbose)
}
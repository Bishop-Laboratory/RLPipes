#' Create an RSeq Data Set
#'
#' Creates an RSeq data object that can be processed with the RSeq pipeline
#'
#' \bold{Note}: all arguments can be specified directly or passed via the \code{samples} object. Additionally, "RNH" is a
#' special phrase which should not be used in the naming of experimental samples. Finally, please run a
#' \code{test_run} to ensure correct parameters prior to analysis. See details below.
#'
#' @param mode Required. The \code{mode} parameter controls the workflow used to analyze data based
#' on best-practices for the specified R-loop mapping method. This can be overwritten
#' by manually passing parameters to \code{analysis_params}. Options include
#' \code{DRIP}, \code{DRIPc}, \code{qDRIP}, \code{ssDRIP}, \code{RDIP}, \code{R-ChIP},
#' \code{RR-ChIP}, and \code{MapR}. Please read the vignette for more details.
#'
#' @param experiment Character cotaining the GEO or SRA accession(s) of experimental samples
#' to download and analyze. Can also be the file prefix or file name of local fastq or bam file(s).
#' Can also be a directory path containing local fastq or bam files if not using public data accessions.
#' If no \code{experiment} names are provided, file prefixes will be treated as \code{experiment}
#' names for subsequent analyses. Acceptable inputs are any combination of GEO or SRA accessions
#' (GSM, GSE, SRX, SRP, SRR, SRS) with any combination of file prefixes, file names, and directories.
#' Must be specified in \code{samples} if not provided here. See \code{details} for
#' more information.
#'
#' @param control Same as \code{experiment} but for control/input samples.
#' Provide in the same order as \code{experiment} entries so that there is a 1:1 match
#' between each \code{experiment} sample and the corresponding control.
#'
#' @param condition Character specifying the biological condition of each sample (e.g., 'HeLa_WT_siCTR').
#' Should be provided in the same order as \code{experiment} entries so that there is a 1:1 match.
#' In the case of folders/studies, condition will be appended to the file prefix if provided. \code{condition} will be
#' used to evaluate the congruency of biological replicates. \bold{NOTE}: RNaseH experimental samples should not be included in
#' the same condition as experimental samples unless they are explicitely defined via \code{replicate} or \code{outname}.
#' In this case, they append the condition with the phase "RNH" (e.g., 'HeLa_WT_siCTR_RNH'). If this is not provided, it will be automatically
#' generated numerically.
#'
#' @param replicate Numeric specifying the biological replicate number or sample type. '-1' identifies RNaseH1 samples which
#' are passed as \code{experiment} rather than control samples. This will be used in combination with \code{condition} during R-loop mapping QC analysis to
#' help evaluate the congruency of replicates and the specifity of mapping. If not provided, \code{RSeq} will attempt to
#' determine RNaseH status from the \code{outname} parameter and sample metadata (for public samples). This can
#' be disabled using the \code{RNH_check} parameter. A warning will
#' be issued if the algorithm was unable to identify an RNaseH sample.
#'
#' @param RNH_check Logical. Should \code{RSeq} attempt to identify and analyze RNaseH1-treated samples for
#'  specificity analysis? Default: TRUE.
#'
#' @param genome Character specifying which genome assembly to use (e.g., "hg19", "mm9"). By default \code{RSeq} uses
#' the latest available pre-built genome annotations corresponding to the detected species (if using public data). If using
#' local data, this information must be supplied. \code{RSeq} will download and build the appropriate annotation set in
#' all cases where that is possible. If not possible or not desired, use \code{custom_genome}.
#'
#' @param custom_genome Character specifying a path to the folder containing \code{bwa} indexes for the custom
#' genome or idices for another aligner if specified in 'analysis_params'.
#' If a length-two character vector, second item should specify the path to a custom gene annotation set in bed format
#' in which the genes are the 4th column (e.g., c("data/myo_myotis_bwa_indices/", "data/myo_myotis_genes.bed"))
#'
#' @param outname Character specifying name (in the case of individual samples) or prefix (in the case of
#' folders/studies with multiple samples) to be used for plotting, reports, and output files. If not supplied,
#' \code{RSeq} will attempt to use the condition and replicate arguments (if unique) or the GEO or SRA sample
#' name if available and the file prefixes in the case of local files. If \code{replicate} does not indicate RNaseH
#'
#' @param outdir Character specifying folder locations to save data to. Default is \code{'RSeq_out'}. Can
#' be specified in \code{samples}.
#'
#' @param operations Character specifying which analysis steps should be included. Default is \code{'all'} which
#' includes \code{c('fastq_QC', 'alignment', 'post_alignment', 'call_peaks', 'annotate', 'rloop_QC')}. Can be specified
#' in \code{analysis_params}. See \code{details} below.
#'
#' @param analysis_params Nested named list used to override default parameters. See \code{details} below.
#'
#' @param samples Data frame containing sample info and analysis parameters. See \code{details} below.
#'
#' @param test_run Logical. If TRUE, will perform pre-flight checks and attempt all operations with a limited subset of data. Default: TRUE.
#'
#' @return An RSeq object containing \code{params}, \code{QC}, \code{report}, and \code{peaks}. See \code{details} below.
#'
#' @details
#'
#' \subsection{Sample input parameters} {
#'   \code{RSeq} can be used with or without specifying a \code{samples} frame. Users will find that this file
#'   provides convenience for complicated analysis. Briefly, \code{samples} should be a data frame of the form
#'   \code{sample X parameter}. The purpose is to allow for sample to be run with its own set of unique parameters if desired.
#'   Parameters specified in the \code{samples} frame supersede those specified in the typical \code{RSeq} arguments.
#'   If columns are not provided in the \code{samples} frame but are provided in arguments to \code{RSeq}, those arguments
#'   will be applied added as a new column.
#'
#'   Possible parameters:
#'   \describe{
#'     \item{experiment}{\strong{Required} column giving the accession names (GSM, GSE, SRX, SRS, SRR, SRP),
#'     file prefixes, bam files, fastq files, or file folders signifying samples or samples collections to analyze. When study accessions (GSE or SRP) are
#'     specified, then all the samples in that study will be downloaded and anlyzed. \code{RSeq} will not attempt to
#'     distinguish between R-loop-mapping samples and other sample types in that accession. It will assign sample names to
#'     the available GSM or SRX accession. Likewise, if a folder directory is specified, \code{RSeq} will blindly process
#'     all the fastq and/or bam files in that folder regardless of type. It will also assign sample names as the file prefix
#'     of the files in that folder. When individual bam or fastq files are specified, their sample names will be identified as the
#'     file prefix. Any combination of accessions, folders, prefixes, and individual files may be specified. It is worthwhile to perform
#'     a \code{test_run} to verify appropriate naming prior to analysis.
#'     }
#'     \item{control}{Column giving the sample accession(s) (GSM,  SRX, SRS, SRR), file prefixes, bam files, or fastq files
#'     which signify the control sample (e.g., genomic input) for the corresponding experimental sample. \bold{Note}: it is NOT
#'     possible to specify a file folder or study accession as a control sample. This is because control samples need to be unique
#'     samples which specifically correspond to experimental samples. \bold{However}: it IS possible to specify one control sample
#'     for a folder or study accession (e.g., SRS199908 is the control sample for SRP208881). In that case, the control sample
#'     will be individually applied to every sample in the folder/study and will be removed from the list of 'experimental' samples.
#'     This is particularly convenient in cases where a single genomic input sample was produced for an entire study.
#'     }
#'     \item{condition}{specifies the biological condition of each sample. This is not required but will significantly improve the quality
#'     of QC analysis. \code{RSeq} will use condition information to assess replicate congruity and mapping specifity}
#'     \item{replicate}{specifies the biological replicate within each condition. '-1' specifies an RNaseH-treated sample which
#'     can be used for specificity analysis. Assigning replicates is only possible for individual samples. For studies/folders,
#'     \code{RSeq} will attempt to determine RNaseH status from the \code{outname} argument and sample metadata (if available).}
#'     \item{outname} {Specifies the name used for each sample in output files, reports, and plotting. In the case of single samples, provide
#'     an output name directly. For folders or study accessions, a prefix and be supplied. If condition and replicate are supplied and are unique in combination
#'     then they will constitue the \code{outname} unless explicitely specified here. In attempting to determine RNaseH status, \code{RSeq} will
#'     check the outname parameter for the phrase 'RNH' if these samples are not specified otherwise.}
#'     \item{RNH_check} {Logical specifying whether \code{RSeq} should check for RNaseH samples.}
#'     \item{outdir} { Character specifying folder locations to save data to. Default is \code{'RSeq_out'}. Can
#'     be specified in \code{samples}}
#'   }
#'
#' }
#'
#' \subsection{Customizing the RSeq pipeline}{
#'   To offer flexibility, the \code{RSeq} pipeline allows users to fully control the workflow parameters.
#'   This can be done simply by specifying \code{operations} and \code{analysis_params}
#'
#'   \strong{operations}:
#'
#'   The default pipeline for each \code{mode} follows the same pattern:
#'   \describe{
#'     \item{fastq_QC}{Fastq reads are trimmed and filtered. Performed using \code{fastp} by default, which also provides an HTML report of read quality.}
#'     \item{alignment}{Alignment of reads to genome. Uses \code{bwa mem} with default parameters for alignment. The output is piped to \code{samtools}
#'     for conversion to bam format and duplicate removal, sorting, and indexing.}
#'     \item{post_alignment}{Calculation of alignment metrics, alignment specifity,
#'     and sample-sample correlation. Uses \code{samtools} for alignment metrics. Then,
#'     \code{deepTools} is used to construct a bin matrix and calculate alignment coverage/specificity,
#'      sample-sample correlation, principal component analysis, fragment size distribution. \code{deepTools}
#'      will also be used to calculate peak pileup around the gene body and R-loop forming sequences and
#'      to generate a bigWig coverage track using RPKM as the signal value in 10bp bins.}
#'     \item{call_peaks}{Calculates R-loop peaks using \code{macs2} with the --broad flag by default. If single-end,
#'     fragment size is estimated before peak calling, using macs2 default params}
#'     \item{annotate}{Annotate peaks with respect to genomic features (e.g., promoters) and genes. Uses \code{ChIPpeakAnno} for
#'     calculating gene annotations and \code{ChIPseeker} for calculating genomic feature overlap}
#'     \item{rloop_QC}{Based on best practices, this step involves calculating a range of quality metrics and assigning
#'     a score based on the outcome of those assessments. For additional details, please see the vignette.}
#'   }
#'
#'
#'   \strong{analysis_params}:
#'
#'   Nested list for describing alterations to the default pipeline.
#'   List
#'
#'
#' }
#'
#'
#' @examples
#' genesOfInterest <- c("ATM", "SLC7A11")
#' correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = genesOfInterest,
#'                               returnDataOnly = TRUE,
#'                               GSEA_Type = "simple",
#'                               Sample_Type = c("normal", "cancer"),
#'                               Tissue = c("respiratory", "pancreas"))
#'
#' genesOfInterest <- c("BRCA1")
#' correlationAnalyzeR::analyzeSingleGenes(genesOfInterest = genesOfInterest,
#'                               GSEA_Type = "simple", returnDataOnly = TRUE,
#'                               crossCompareMode = TRUE,
#'                               whichCompareGroups = "normal")
#'
#' @export
makeRSeqDataSet <- function(mode = NULL,
                            experiment = NULL,
                            control = NULL,
                            condition = NULL,
                            replicate = NULL,
                            RNH_check = TRUE,
                            convert_bams = FALSE,
                            cores = NULL,
                            genome = NULL,
                            outname = NULL,
                            genome_dir = "~/.RSeq_genomes/",
                            outdir = '../RSeq_out',
                            analysis_params = NULL,
                            samples = NULL) {

  # Bug testing
  setwd("modules")
  mode = NULL
  experiment = NULL
  control = NULL
  outdir = NULL
  analysis_params = NULL
  condition = NULL
  replicate = NULL
  RNH_check = TRUE
  genome = NULL
  genome_home_dir = file.path(path.expand("~"), ".RSeq_genomes")
  convert_bams = TRUE
  out_name = NULL
  cores = NULL
  source("R/utils.R")
  ## For bug testing ##
  samples <- data.frame(
    "experiment" = c("SRX2481503", "GSM2326832", "SRX2918366", "GSM3937232", "GSM3936514", "GSM3936515", "GSM1720615",
                     "GSM3936516", "SRX5129665", "GSM2104456", "GSM1720613", "GSM2550993",
                     "~/Bishop.lab/Preprocessing/DRIP_Seq/GSE145964_Developmental_Context_ChIP_DRIP/Data/Raw_Reads/SRR11185284_1.fastq+~/Bishop.lab/Preprocessing/DRIP_Seq/GSE145964_Developmental_Context_ChIP_DRIP/Data/Raw_Reads/SRR11185284_2.fastq",
                     "../EUFA_BRCA2/"),
    "control" = c("SRX2481504", "GSM2326824", "SRX2918367", NA, "GSM3936517", "GSM3936517", NA,
                  "GSM3936517", "SRX5129664", NA, NA, "GSM2550995",
                  NA, NA),
    "genome" = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                 "hg38", "hg38"),
    "strand_specific" = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE),
    stringsAsFactors = F
  )
  write.csv(samples, file = "../instance/sample_sheet_example.csv")
  #####################

  # Additional data
  read_align_pattern <- "\\.bam$|\\.f[ast]*q$"
  public_pattern <- paste0("^SR[RAXSP][0-9]+$|^GS[EM][0-9]+$|^PRJNA[0-9]+$|^SAMN[0-9]+$")
  R1_pattern <- "[\\._]+[rR]*1\\.f[ast]*q$"
  R2_pattern <- "[\\._]+[rR]*2\\.f[ast]*q$"
  fq_end_pattern <- paste0(R1_pattern, "|", R2_pattern, "|[\\._]+[rR]*\\*.f[ast]*q$|\\.f[ast]*q$")

  # Parse arguments
  tryCatch(stopifnot(! is.null(samples) | ! is.null(experiment)),
           error = function(x) stop("Must provide 'experiment' or 'samples'"))
  if (! "cores" %in% colnames(samples) &
      is.null(cores)) samples$cores <- parallel::detectCores()/2
  if (! "RNH_check" %in% colnames(samples)) samples$RNH_check <- RNH_check
  if (! "condition" %in% colnames(samples)) {samples$condition <- NA}
  if (! "paired_end" %in% colnames(samples)) {samples$paired_end <- NA}
  if (! "replicate" %in% colnames(samples)) {samples$replicate <- NA}
  if (! "effective_genome_size" %in% colnames(samples)) {samples$effective_genome_size <- NA}
  if (! "full_genome_length" %in% colnames(samples)) {samples$full_genome_length <- NA}
  if (! "genome_home_dir" %in% colnames(samples)) {samples$genome_home_dir <- genome_home_dir}

  # Mark files, dirs, and 'other'
  dirIndExp <- which(dir.exists(samples$experiment))
  fileIndExp <- grep(samples$experiment, pattern = read_align_pattern)
  publicIndExp <- grep(samples$experiment, pattern = public_pattern)

  # Must provide assembly for local files
  assembly <- samples$genome[c(fileIndExp, dirIndExp)]
  if (any(is.null(assembly))) {
    stop("If providing local files/directories, specify the genome assembly for those samples.")
  }
  badGenomeSpec <- which(is.na(samples$genome_home_dir) & is.na(samples$genome))
  if (length(badGenomeSpec)) {
    stop("User must supply a genome directory location under 'genome' or a genomes home folder as 'genome_home_dir'.",
         " For sample(s) ", samples$experiment[badGenomeSpec], ", the genome supplied is ", samples$genome[badGenomeSpec],
         " and the genome home directory is ", samples$genome_home_dir[badGenomeSpec], "../..")
  }

  # Set file type and attempt to set sample_names
  samples$file_type <- NA
  samples$file_type[publicIndExp] <- "public"
  samples$file_type[fileIndExp][grep(samples$experiment[fileIndExp],
                                     pattern = "\\.f[ast]*q$")] <- "fastq"
  samples$file_type[fileIndExp][grep(samples$experiment[fileIndExp],
                                     pattern = "\\.bam$")] <- "bam"
  samples$sample_name <- paste0(samples$condition, "_", samples$replicate)
  samples$sample_name[samples$sample_name %in%
                        samples$sample_name[duplicated(samples$sample_name)]] <- NA
  samples$read_length <- NA

  # Convert directories into file lists
  local_samples_list <- lapply(dirIndExp, FUN = function(ind_now) {
    dir_files <- list.files(samples$experiment[ind_now], recursive = TRUE,
                            pattern = read_align_pattern, full.names = TRUE)
    if (length(dir_files)) {
      samples_names_dir <- gsub(basename(dir_files),
                                pattern = read_align_pattern, replacement = "")

      # Handle fastqs
      fastqs <- grep(dir_files, pattern = "\\.f[ast]*q$")
      fq_files <- dir_files[fastqs]
      fq_sample_names <- gsub(fq_files, pattern = fq_end_pattern, replacement = "")
      table_df <- as.data.frame(table(fq_sample_names), stringsAsFactors = FALSE)
      potential_mates <- table_df$fq_sample_names[table_df$Freq == 2]
      single_ends <- table_df$fq_sample_names[table_df$Freq != 2]
      paired_ends <- c()
      for (pair in potential_mates) {
        possible_pair <- fq_files[which(fq_sample_names %in% pair)]
        fastq_1 <- possible_pair[grep(possible_pair, pattern = R1_pattern)]
        fastq_2 <- possible_pair[grep(possible_pair, pattern = R2_pattern)]
        if (! length(fastq_1) | ! length(fastq_2)) {
          single_ends <- c(single_ends, pair)
        } else {
          paired_ends <- c(paired_ends, paste0(fastq_1, "+", fastq_2))
          warning("RSeq identified ", fastq_1, " & ", fastq_2," as a mate pair based on their file names.",
                  " If this was incorrect, please disable 'smart_pairing' or specify files directly")
        }
      }


      single_end_inds <- dir_files[which(dir_files %in% fq_files[which(fq_sample_names %in% single_ends)])]
      dir_files_fq <- c(single_end_inds, paired_ends)

      # Handle bams
      bams <- grep(dir_files, pattern = "\\.bam$")
      dir_files_bam <- dir_files[bams]

      # Final naming
      dir_files_out <- c(dir_files_fq, dir_files_bam)
      samples_names_dir_out <- gsub(basename(dir_files_out),
                                    pattern = paste0(fq_end_pattern, "|", read_align_pattern),
                                    replacement = "")
      file_types <-  character(length = length(dir_files_out))
      file_types[bams] <- "bam"
      file_types[! file_types == "bam"] <- "fastq"

      # Add to samples frame
      local_samples <- samples[rep(ind_now, length(dir_files_out)),, drop = FALSE]
      local_samples$experiment <- dir_files_out
      local_samples$file_type <- file_types
      local_samples$sample_name <- samples_names_dir_out
      if (all(is.na(local_samples$condition))) {
        warning("Sample conditions cannot be determined for directory inputs.",
                " File conditions will be the name of the folder",
                " they are located in.")
        local_samples$condition <- basename(dirname(gsub(local_samples$experiment,
                                                         pattern = "\\.fastq\\+.*", replacement = "")))

      }
      return(local_samples)
    } else {
      stop("Directory ", samples$experiment[ind_now],
           " does not contain and bam or fastq files!")
    }
  })

  local_samples_frame <- data.table::rbindlist(local_samples_list)
  if (dim(local_samples_frame)[1] != 0) {
    samples <- rbind(local_samples_frame, samples[c(-dirIndExp),])
  }


  ## Handle local files ##
  # bams
  bamInd <- which(samples$file_type == "bam")
  if (length(bamInd)) {

    # -- Check bam files against assembly entered
    genome_check <- unlist(lapply(bamInd, function(ind) {
      check_bam_assembly(samples$experiment[ind], samples$genome[ind])
    }))
    if (any(! genome_check) & ! convert_bams) {
      msg <- paste0("Specified assembly (", samples$genome[bamInd], ") does not match bam file header (",
                    samples$experiment[bamInd][! genome_check], ").",
                    " Select 'convert bams = TRUE' to re-align reads to ",
                    samples$genome[bamInd][! genome_check], ".\n")
      stop(msg)
    } else if (any(! genome_check)) {
      msg <- paste0("Specified assembly (", samples$genome[bamInd][! genome_check],
                    ") does not match bam file header (",
                    samples$experiment[bamInd][! genome_check], ").",
                    " Converting and re-aligning reads to ",
                    samples$genome[bamInd][! genome_check], ".\n")
      warning(msg, immediate. = TRUE)
    }

    # -- Get the read lengths
    samples$read_length[bamInd] <- unlist(lapply(bamInd, function(ind) {
      get_bam_read_length(samples$experiment[ind])
    }))


    # -- create sample_name and condition arguments
    samples$condition[bamInd][is.na(samples$condition[bamInd])] <-
      basename(dirname(samples$experiment[bamInd][is.na(samples$condition[bamInd])]))
    samples$sample_name[bamInd][is.na(samples$sample_name[bamInd])] <-
      gsub(basename(samples$experiment[bamInd][is.na(samples$sample_name[bamInd])]),
           pattern = "\\.bam$", replacement = "")
  }

  # -- fastqs
  fqInd <- which(samples$file_type == "fastq")
  for (ind in fqInd) {
    fqNow <- samples$experiment[ind]

    plus_bool <- grepl(fqNow[1], pattern = "\\.f[ast]*q(\\+)")
    wild_bool <- grepl(fqNow[1], pattern = "\\*\\.f[ast]*q$")
    if (plus_bool) {
      fastq_1 <- gsub(fqNow, pattern = "(.+\\.f[ast]*q)\\+.+", replacement = "\\1")
      fastq_2 <- gsub(fqNow, pattern = ".+\\.f[ast]*q\\+(.+)", replacement = "\\1")
      samples$sample_name[ind] <- gsub(basename(fastq_1),
                                       pattern = fq_end_pattern, replacement = "")
      if (! is.na(samples$condition[ind]) & ! is.na(samples$replicate[ind])) {
        samples$sample_name <- paste0(samples$condition, "_", samples$replicate)
      }
      samples$paired_end[ind] <- TRUE
      samples$read_length[ind] <- get_fastq_read_length(fastq_1)
    } else if (wild_bool) {
      fqFiles <- list.files(path = dirname(fqNow), full.names = TRUE,
                            pattern = basename(gsub(fqNow, pattern = "\\*", replacement = "../..")))
      if (length(fqFiles) == 2) {
        warning("RSeq identified using the wildcard '*' as a mate pair: ", fqFiles[1], " & ", fqFiles[2], "../..",
                " If this was incorrect, please remove the wildcard and specify mates directly using",
                " 'sample_name_1.fastq+sample_name_2.fastq' nomenclature")
        fastq_1 <- fqFiles[grep(fqFiles, pattern = R1_pattern)]
        fastq_2 <- fqFiles[grep(fqFiles, pattern = R2_pattern)]
        if (! length(fastq_1) | ! length(fastq_2)) {
          stop("Could not assign files ", paste0(fqFiles, collapse = " and ")," to read 1 and 2. Please ",
               "ensure that file names follow the form 'sample_name{_1|_R1|.1|.R1}.{fastq|fq}'",
               " for read 1 and 'sample_name{_2|_R2|.2|.R2}.{fastq|fq}' for read 2. Or specify mates explicitly",
               " using 'sample_name_1.fastq+sample_name_2.fastq' nomenclature. Contact package ",
               "maintainer if more naming conventions would be useful.")
        }
        samples$sample_name[ind] <- gsub(basename(fastq_1),
                                         pattern = fq_end_pattern, replacement = "")
        if (! is.na(samples$condition[ind]) & ! is.na(samples$replicate[ind])) {
          samples$sample_name <- paste0(samples$condition, "_", samples$replicate)
        }
        samples$experiment[ind] <- paste0(fastq_1, "+", fastq_2)
        samples$paired_end[ind] <- TRUE
        samples$read_length[ind] <- get_fastq_read_length(fastq_1)
      } else {
        stop("Could not assign files ", paste0(fqFiles, collapse = " and ")," to read 1 and 2. Please ",
             "ensure that file names follow the form 'sample_name{_1|_R1|.1|.R1}.{fastq|fq}'",
             " for read 1 and 'sample_name{_2|_R2|.2|.R2}.{fastq|fq}' for read 2. Or specify mates explicitly",
             " using 'sample_name_1.fastq+sample_name_2.fastq' nomenclature. Contact package ",
             "maintainer if more naming conventions would be useful.")
      }
    } else {
      samples$sample_name[ind] <- gsub(basename(samples$experiment[ind]),
                                       pattern = fq_end_pattern, replacement = "")
      samples$read_length[ind] <- get_fastq_read_length(samples$experiment[ind])
      samples$paired_end[ind] <- FALSE

    }
  }

  # Check for out_names
  if (! 'out_name' %in% colnames(samples)) {
    samples$out_name <- gsub(samples$sample_name, pattern = "_", replacement = " ")
  } else {
    samples$out_name[is.na(samples$out_name) | is.null(samples$out_name)] <-
      gsub(samples$sample_name[is.na(samples$out_name) | is.null(samples$out_name)],
           pattern = "_", replacement = " ")
  }

  ## Process public samples ##
  if (any(samples$file_type == "public")) {
    samples_public <- samples[samples$file_type == "public",]

    # Get info for controls first
    public_ctr_accessions <- samples_public$control[! is.na(samples_public$control)]
    ctr_studies <- grep(public_ctr_accessions, pattern = "^GSE[0-9]+$|^SRP[0-9]+$|^PRJNA[0-9]+$")
    if (length(ctr_studies)) {
      stop("User input incorrect (", unique(public_ctr_accessions[ctr_studies]),"). ",
           "Control samples must be samples (e.g., SRX, SRR, GSM, SAMN) and ",
           "not studies (e.g., GSE, SRP, PRJNA).")
    }

    if (length(public_ctr_accessions)) {
      sra_info_ctr <- get_public_run_info(public_ctr_accessions)
    }

    # -- add back to public samples frame
    new_ctr_vec <- c()
    for (orig_ctr in samples_public$control) {
      if (is.na(orig_ctr)) {
        new_ctr_vec <- c(new_ctr_vec, NA)
      } else {
        new_ctr_vec <- c(new_ctr_vec, as.character(
          sra_info_ctr$experiment[sra_info_ctr$accessions_original == orig_ctr]))
      }
    }
    samples_public$control <- new_ctr_vec

    # Get info for experimental samples
    sra_info_exp <- get_public_run_info(samples_public$experiment)
    # -- collate
    public_sample_list <- list()
    for (i in 1:length(samples_public$experiment)) {
      #sample <- samples_public$experiment[3]
      sample <- samples_public$experiment[i]
      samples_public_now <- samples_public[samples_public$experiment == sample,]
      sra_info_now <- sra_info_exp[sra_info_exp$accessions_original == sample,]
      colnames(sra_info_now)[4] <- "experiment_final"
      sra_info_new <- merge(x = as.data.frame(samples_public_now),
                                y = as.data.frame(sra_info_now),
                            by.x = "experiment", by.y = "accessions_original", all = TRUE)
      # -- Fix EXP info
      sra_info_new <- sra_info_new[,c(-1)]
      colnames(sra_info_new)[which(colnames(sra_info_new) == "experiment_final")] <- "experiment"

      # TODO: make sure that user-specified genome is correct
      if (all(is.na(sra_info_new$genome.x))) {
        sra_info_new$genome <- sra_info_new$genome.y
      } else {
        if (! all(sra_info_new$genome.x == sra_info_new$genome.y)) {
          warning("RSeq has detected the newest genome assembly for ",
                  unique(sra_info_new$experiment), " as ", unique(sra_info_new$genome.y), "
                  but user has specified ", unique(sra_info_new$genome.x), "../.. ",
                  unique(sra_info_new$genome.x), " will be used instead.")
        }
        sra_info_new$genome <- sra_info_new$genome.x
      }

      if (all(is.na(sra_info_new$sample_name.x))) {
        sra_info_new$sample_name <- sra_info_new$sample_name.y
      }

      # -- sample name
      sra_info_new$sample_name <- ifelse(is.na(sra_info_new$sample_name.x),
                                         as.character(sra_info_new$sample_name.y),
                                         as.character(sra_info_new$sample_name.x))
      # -- condition
      sra_info_new$condition <- ifelse(is.na(sra_info_new$condition.x),
                                       as.character(sra_info_new$condition.y),
                                       as.character(sra_info_new$condition.x))
      # -- outname
      sra_info_new$out_name <- ifelse(is.na(sra_info_new$out_name.x),
                                      as.character(sra_info_new$out_name.y),
                                      as.character(sra_info_new$out_name.x))
      # -- paired end
      sra_info_new$paired_end <- sra_info_new$paired_end.y

      # -- read length
      sra_info_new$read_length <- sra_info_new$read_length.y

      # -- final ordering
      keepInd <- which(colnames(sra_info_new) %in% colnames(samples_public))
      sra_info_new <- as.data.frame(sra_info_new)[,keepInd]
      if(! all(colnames(samples_public) %in% colnames(sra_info_new))) {
        stop("Bug at public sample collation -- please notify package author. Should be unreachable.")
      }
      sra_info_new <- sra_info_new[,order(match(colnames(sra_info_new), colnames(samples_public)))]
      if(! all(colnames(samples_public) == colnames(sra_info_new))) {
        stop("Bug at public sample collation -- please notify package author. Should be unreachable.")
      }
      public_sample_list[[sample]] <- sra_info_new
    }
    samples_public <- data.table::rbindlist(public_sample_list)
    samples_private <- samples[samples$file_type != "public",]
    if(! all(colnames(samples_public) == colnames(samples_private))) {
      stop("Bug at public sample collation -- please notify package author. Should be unreachable.")
    }
    samples <- rbind(samples_public, samples_private)
  }

  if (any(duplicated(samples$sample_name))) {
    stop("Sample names must be unique. Please remove duplicates: ",
         paste0(samples$sample_name[duplicated(samples$sample_name)], collapse = ", "))
  }

  # Set directories and file output names
  if (! "outdir" %in% colnames(samples)) {samples$outdir <- "../RSeq_out" } else {
    samples$outdir[is.na(samples$outdir) || is.null(samples$outdir) ||
                     samples$outdir == ""] <- "../RSeq_out"
  }
  # -- fastqs
  fastq_dir <- file.path(samples$outdir, "fastqs")
  fastq_dir[samples$file_type == "fastq"] <-
    dirname(as.character(samples$experiment)[samples$file_type == "fastq"])
  samples$fastq_dir <- fastq_dir
  fastq_dir <- file.path(samples$outdir, "fastqs")
  fastq_dir[samples$file_type == "fastq"] <-
    dirname(as.character(samples$experiment)[samples$file_type == "fastq"])
  fastq_file <- file.path(fastq_dir, paste0(samples$sample_name, ".fastq"))
  fastq_file_1 <- file.path(fastq_dir, paste0(samples$sample_name, "_1.fastq"))
  fastq_file_2 <- file.path(fastq_dir, paste0(samples$sample_name, "_2.fastq"))
  fastq_file_1[samples$file_type == "fastq"] <-
    gsub(samples$experiment[samples$file_type == "fastq"],
         pattern = "(.+)\\+(.+)", replacement = "\\1")
  fastq_file_2[samples$file_type == "fastq"] <-
    gsub(samples$experiment[samples$file_type == "fastq"],
         pattern = "(.+)\\+(.+)", replacement = "\\2")
  fastq_file_2[! samples$paired_end] <- NA
  samples$fastq_dir <- fastq_dir
  samples$fastq_file <- fastq_file
  samples$fastq_file_1 <- fastq_file_1
  samples$fastq_file_2 <- fastq_file_2

  # -- bams
  bam_dir <- file.path(samples$outdir, "bams")
  bam_dir[samples$file_type == "bam"] <-
    dirname(as.character(samples$experiment)[samples$file_type == "bam"])
  bam_file <- file.path(bam_dir, paste0(samples$sample_name, ".bam"))
  bam_file[samples$file_type == "bam"] <- samples$experiment[samples$file_type == "bam"]
  samples$bam_dir <- bam_dir
  samples$bam_file <- bam_file

  # Other
  samples$peak_dir <- file.path(samples$outdir, "peaks")
  samples$coverage_dir <- file.path(samples$outdir, "coverage")
  samples$QC_dir <- file.path(samples$outdir, "QC")
  samples$bam_QC_dir <- file.path(samples$QC_dir, "bams")
  samples$fastq_QC_dir <- file.path(samples$QC_dir, "fastqs")
  samples$downstream_dir <- file.path(samples$outdir, "downstream")
  samples$report_dir <- file.path(samples$outdir, "report")

  # Misc
  RSeq_version <- paste0(unlist(utils::packageVersion("RSeq")), collapse = "../..")

  # Get effective genome size
  available_genomes <- RSeq::available_genomes
  rownames(available_genomes) <- available_genomes$UCSC_orgID
  available_genome_sizes <- available_genomes[,grep(colnames(available_genomes), pattern = "eff_genome_size")]
  sizes <- as.numeric(gsub(colnames(available_genome_sizes),
                           pattern = "eff_genome_size_([0-9]+)bp",  replacement = "\\1"))
  samples$effective_genome_size[is.na(samples$effective_genome_size)] <-
    unlist(lapply(samples$sample_name[is.na(samples$effective_genome_size)], function(sample_now) {
    len_now <- samples$read_length[samples$sample_name == sample_now]
    colInd <- which.min(abs(sizes-len_now))
    genome_now <- samples$genome[samples$sample_name == sample_now]
    eff_genome_size <- available_genome_sizes[genome_now, colInd]
    if (is.na(eff_genome_size)) {
      stop("Genome size is not available for ", genome_now, ". Please manually supply effective genome size as a",
           "column in 'samples' or variable in 'makeRSeqDataSet'.",
           " Genome sizes may be calculated with the 'unique-kmers.py script from the khmer package'")
    }
    eff_genome_size
  }))
  samples$full_genome_length <-
    unlist(lapply(samples$sample_name[is.na(samples$full_genome_length)], function(sample_now) {
    genome_now <- samples$genome[samples$sample_name == sample_now]
    genome_length <- RSeq::available_genomes$genome_length[RSeq::available_genomes$UCSC_orgID ==
                                                               genome_now]
    if (is.na(genome_length)) {
      stop("Genome length is not available for ", genome_now, ". Please manually supply genome length")
    }
    genome_length
  }))

  ## Generate RSeq_variable object for each sample
  vars_list <- lapply(samples$sample_name, FUN = function(exp_now) {
    # # Bug testing
    #exp_now <- samples$sample_name[1]

    vars_now <- new("RSeq_variables")
    samples_now <- samples[samples$sample_name == exp_now, , drop = TRUE]

    # file and directory arguments
    vars_now@fastq_1_available <- file.exists(samples_now$fastq_file_1)
    vars_now@fastq_2_available <- file.exists(samples_now$fastq_file_2)
    vars_now@fastq_available <- file.exists(samples_now$fastq_file)
    vars_now@fastq_dir <- samples_now$fastq_dir
    vars_now@fastq_QC_dir <- samples_now$fastq_QC_dir
    vars_now@bam_available <- file.exists(samples_now$bam_file)
    vars_now@bam_file <- samples_now$bam_file
    vars_now@bam_dir <- samples_now$bam_dir
    vars_now@bam_QC_dir <- samples_now$bam_QC_dir

    # Genome directory
    available_genomes <- RSeq::available_genomes
    vars_now@genome <- as.character(samples_now$genome)
    if (dir.exists(vars_now@genome)) {
      vars_now@genome_home_dir <- vars_now@genome
      if (is.na(samples_now$effective_genome_size)) {
        # TODO: Run unique-kmers.py and save as a txt file in that dir if missing.
        stop("No effective genome size provided for user-supplied genome: ", vars_now@genome, "../..",
             " Please calculate effective genome sizes using unique-kmers.py from the khmer package",
             " and supply as a column in 'samples' or as a vector to 'makeRSeqDataSet'.")
      }
    } else {
      vars_now@genome_home_dir <- samples_now$genome_home_dir
      vars_now@genome_fasta_file <- file.path(vars_now@genome_home_dir, paste0(vars_now@genome, ".fa.gz"))
      vars_now@genome_annotations_file <- file.path(vars_now@genome_home_dir, paste0(vars_now@genome, ".refGene.gtf.gz"))
      vars_now@genome_fasta_file_exists <- file.exists(vars_now@genome_fasta_file)
      vars_now@genome_annotations_file_exists <- file.exists(vars_now@genome_annotations_file)
      vars_now@genome_fasta_url <- paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/",
                                          vars_now@genome, "/bigZips/", vars_now@genome, ".fa.gz")
      vars_now@available_genome_indices <- list.dirs(vars_now@genome_home_dir, recursive = FALSE)
      vars_now@genome_annotations_available <- available_genomes$genes_available[available_genomes$UCSC_orgID ==
                                                                                   vars_now@genome]
      vars_now@genome_annotations_url <- paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/",
                                                vars_now@genome, "/bigZips/genes/", vars_now@genome, ".refGene.gtf.gz")

    }

    # Effective genome size

    # Others
    vars_now@cores <- as.integer(samples_now$cores)
    vars_now@experiments <- stringr::str_split(gsub(samples_now$experiment, pattern = ".+\\+(.+)",
                            replacement = "\\1"), ",")[[1]]
    vars_now@controls <- stringr::str_split(gsub(samples_now$control, pattern = ".+\\+(.+)",
                            replacement = "\\1"), ",")[[1]]
    vars_now@sample_name <- samples_now$sample_name
    vars_now@coverage_dir <- samples_now$coverage_dir
    vars_now@RSeq_version <- RSeq_version
    vars_now@file_format <- samples_now$file_type
    vars_now@downstream_dir <- samples_now$downstream_dir
    vars_now@out_dir <- samples_now$outdir
    vars_now@report_dir <- samples_now$report_dir
    vars_now@paired_end <- samples_now$paired_end
    vars_now@effective_genome_size <- samples_now$effective_genome_size
    vars_now@full_genome_length <- samples_now$full_genome_length
    vars_now@strand_specific <- samples_now$strand_specific
    vars_now@sample_name <- samples_now$sample_name
    vars_now@out_name <- samples_now$out_name
    slots <- methods::slotNames(vars_now)
    vars_now_list <- purrr::map(slots, ~ methods::slot(object = vars_now, name = .x))
    names(vars_now_list) <- slots
    return(vars_now_list)
  })

  names(vars_list) <- samples$sample_name

  jsonlite::write_json(vars_list, path = "../run_vars.json")
  #analysis_params <- validate_params(analysis_params)
  #rsds <- new("RSeqDataSet", samples, analysis_params)
  #return(analyze_RSeq(rsds))
}








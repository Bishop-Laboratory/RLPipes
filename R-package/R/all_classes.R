#' @rdname RSeq_variables
#' @export
setClass("RSeq_variables",
         slots = list(
           sample_name = "character",
           out_name = "character",
           out_dir = "character",
           peak_available = "logical",
           bam_available = "logical",
           bam_file = "character",
           file_format = "character",
           fastq_available = "logical",
           fastq_1_available = "logical",
           fastq_2_available = "logical",
           fastq_dir = "character",
           fastq_QC_dir = "character",
           fastq_file = "character",
           fastq_file_1 = "character",
           fastq_file_2 = "character",
           genome_home_dir = "character",
           genome = "character",
           genome_fasta_file = "character",
           genome_annotations_file = "character",
           genome_fasta_file_exists = "logical",
           genome_annotations_file_exists = "logical",
           genome_fasta_url = "character",
           available_genome_indices = "character",
           genome_annotations_available = "logical",
           genome_annotations_url = "character",
           genome_available = "logical",
           cores = "integer",
           bam_dir = "character",
           bam_QC_dir = "character",
           peaks_dir = "character",
           full_genome_length = "numeric",
           effective_genome_size = "numeric",
           experiments = "character",
           controls = "character",
           moeity = "character",
           ip_type = "character",
           coverage_dir = "character",
           RSeq_version = "character",
           all_experiment_bams = "character",
           downstream_dir = "character",
           report_dir = "character",
           read_length = "numeric",
           paired_end = "logical",
           strand_specific = "logical",
           all_experiment_labels = "character"
         ))

#' @rdname RSeq_module
#' @export
setClass("RSeq_module",
         slots = list(name = "character",
                      description = "character",
                      software = "list",
                      command = "character",
                      options = "list",
                      input_type = "list", # Ensure module compatibility
                      output_type = "list",  # Ensure module compatibility
                      status = "character",
                      variables = "S4"))

setValidity("RSeq_module", function(object) {

  ## Validate parameter R-package ##
  if (! length(object@name)) {
    stop("RSeq_modules must all contain 'name'.")
  }
  # -- must contain a command
  if (! length(object@command)) {
    stop(paste0("RSeq_module ", slot(object, name = "name"), " is missing "))
  }

  # Get required variables
  module_name <- object@name
  cmd <- object@command
  options <- object@options
  if (grep(pattern = "\\[options\\]", x = cmd)) {
    if (! length(options)) {stop("Command with '[options]' provided in ", module_name,
                                 " but options list is empty!")}
    lapply(option, FUN = function(option_name) {
      option <- options[[4]]
      if (is.list(option)) {
        slot(object, names(option)[1])
        object[["name"]]
        lapply(option, function(x) {stringr::str_match_all(x, pattern = "\\{([a-zA-Z0-9_]+)\\}")[[1]][,2]})
      }
      stringr::str_match_all(option, pattern = "\\{([a-zA-Z0-9_]+)\\}")[[1]][,2]
    })
  }


})

#' @rdname RSeq_workflow
#' @export
setClass("RSeq_workflow",
         slots = list(name = "character",
                      author = "character",
                      description = "character",
                      modules = "list",
                      active = "logical",
                      variables = "S4"))


setValidity("RSeq_workflow", function(object) {

  ## Validate parameter R-package ##
  if (! length(object@name)) {
    stop("RSeq_workflow must all contain 'name'.")
  }
  # -- must contain a command
  if (! length(object@modules)) {
    stop(paste0("RSeq_workflow ", slot(object, name = "name"), " contains no modules!"))
  }

  param_list <- object@modules
  torun <- names(param_list)
  which(! names(object@variables))


  # -- must either download or provide fastq files
  if (! "call_peaks" %in% torun & ! object@variables[["peaks_available"]]) {
    stop(paste0("sample ", object@variables[["applied_to"]], " workflow incorrect. Cannot run fastq-dependent operations -- ",
                "align_reads, bam_PCA, calculate_coverage, etc -- without providing a fastq file or selecting",
                " the 'download_sra' module in work flow."))
  }
  # -- must either align or provide a bam file
  if (! "align_reads" %in% torun & ! object@variables[["bam_available"]]) {
    stop(paste0("sample ", sample_now, " parameters incorrect. Cannot run bam-dependent operations -- ",
                "call_peaks, bam_PCA, calculate_cover, etc -- without providing a bam file or selecting",
                " the 'align_reads' module in parameter list."))
  }
  # -- must either download or provide fastq files
  if (! "download_sra" %in% torun & ! colData[sample_now, "fastq_ready"]) {
    stop(paste0("sample ", sample_now, " parameters incorrect. Cannot run fastq-dependent operations -- ",
                "align_reads, bam_PCA, calculate_coverage, etc -- without providing a fastq file or selecting",
                " the 'download_sra' module in parameter list."))
  }

})

#' @rdname RSeqDataSet
#' @export
setClass("RSeqDataSet",
         contains = "RangedSummarizedExperiment",
         slots = list(analysis_params = "list",
                      test_run = "logical",
                      cores = "integer",
                      verbose = "logical"))

setValidity("RSeqDataSet", function(object) {

  ## Validate analysis params list ##
  analysis_param <- object@analysis_params
  # -- must be one parameter set per sample
  if (! all(names(analysis_params) %in% rownames(colData(object)))) {
    stop("analysis_params must contain parameters for every sample.")
  }

  # -- requirements for each param set
  lapply(names(analysis_params), FUN = function(sample_now) {


  })


  design <- design(object)
  # 'design' is either a formula or matrix
  stopifnot(is(design, "formula") | is(design, "matrix"))

  if (is(design, "formula")) {
    designVars <- all.vars(design)
    if (!all(designVars %in% names(colData(object)))) {
      return("all variables in design formula must be columns in colData")
    }
    designVarsClass <- sapply(designVars, function(v) class(colData(object)[[v]]))
    if (any(designVarsClass == "character")) {
      return("variables in design formula are character vectors.
             convert these columns of colData(object) to factors before including in the design formula")
    }
    designFactors <- designVars[designVarsClass == "factor"]
    # levels would duplicate after make.names()
    if (any(sapply(designFactors,function(v) {
      factor.lvls <- levels(colData(object)[[v]])
      factor.nms <- make.names(factor.lvls)
      any(duplicated(factor.nms))
    }))) {
      return("levels of factors in the design have non-unique level names after make.names() is applied.
             best to only use letters and numbers for levels of factors in the design")
    }
    # levels contain characters other than letters, numbers, and underscore
    if (any(sapply(designFactors,function(v) {
      factor.lvls <- levels(colData(object)[[v]])
      any(!grepl("^[A-Za-z0-9_.]+$",factor.lvls))
    }))) {
      # just a warning for now
      message("  Note: levels of factors in the design contain characters other than
              letters, numbers, '_' and '.'. It is recommended (but not required) to use
              only letters, numbers, and delimiters '_' or '.', as these are safe characters
              for column names in R. [This is a message, not a warning or an error]")
    }
    } else if (is(design, "matrix")) {
      # TODO add some more tests for if 'design' is matrix
      stopifnot(nrow(design) == ncol(object))
    }

  TRUE
    })

#' RSeqDataSet object and constructors
#'
#' \code{RSeqDataSet} is a subclass of \code{RangedSummarizedExperiment},
#' used to store the input values, intermediate calculations and results of an
#' analysis of R-loop mapping.
#' The constructor functions create a RSeqDataSet object
#' from various types of input:
#' a RangedSummarizedExperiment and a list of R-package corresponding to
#' workflows for each sample.
#' See the vignette for examples of construction from different types.
#'
#' Note on the error message "assay colnames() must be NULL or equal colData rownames()":
#' this means that the colnames of countData are different than the rownames of colData.
#' Fix this with: \code{colnames(countData) <- NULL}
#'
#' @param se a \code{RangedSummarizedExperiment} with columns of variables
#' indicating sample information in \code{colData},
#' and the counts as the first element in the assays list, which will
#' be renamed "counts". A \code{RangedSummarizedExperiment} object can be
#' generated by the function \code{summarizeOverlaps} in the GenomicAlignments
#' package.
#' @param design a \code{formula} or \code{matrix}.
#' the \code{formula} expresses how the counts for each gene
#' depend on the variables in \code{colData}. Many R \code{formula} are valid,
#' including designs with multiple variables, e.g., \code{~ group + condition},
#' and designs with interactions, e.g., \code{~ genotype + treatment + genotype:treatment}.
#' See \code{\link{results}} for a variety of designs and how to extract results tables.
#' By default, the functions in this package will use
#' the last variable in the formula for building results tables and plotting.
#' \code{~ 1} can be used for no design, although users need to remember
#' to switch to another design for differential testing.
#' @param countData for matrix input: a matrix of non-negative integers
#' @param colData for matrix input: a \code{DataFrame} or \code{data.frame} with at least a single column.
#' Rows of colData correspond to columns of countData
#' @param tidy for matrix input: whether the first column of countData is the rownames for the count matrix
#' @param sampleTable for htseq-count: a \code{data.frame} with three or more columns. Each row
#' describes one sample. The first column is the sample name, the second column
#' the file name of the count file generated by htseq-count, and the remaining
#' columns are sample metadata which will be stored in \code{colData}
#' @param txi for tximport: the simple list output of the \code{tximport} function
#' @param directory for htseq-count: the directory relative to which the filenames are specified. defaults to current directory
#' @param ignoreRank use of this argument is reserved for DEXSeq developers only.
#' Users will immediately encounter an error upon trying to estimate dispersion
#' using a design with a model matrix which is not full rank.
#' @param ... arguments provided to \code{SummarizedExperiment} including rowRanges and metadata. Note that
#' for Bioconductor 3.1, rowRanges must be a GRanges or GRangesList, with potential metadata columns
#' as a DataFrame accessed and stored with \code{mcols}. If a user wants to store metadata columns
#' about the rows of the countData, but does not have GRanges or GRangesList information,
#' first construct the RSeqDataSet without rowRanges and then add the DataFrame with \code{mcols(dds)}.
#'
#' @return A RSeqDataSet object.
#'
#' @aliases RSeqDataSet RSeqDataSet-class RSeqDataSetFromMatrix RSeqDataSetFromHTSeqCount
#'
#' @references See \url{http://www-huber.embl.de/users/anders/HTSeq} for htseq-count
#'
#' @docType class
#'
#' @examples
#'
#' countData <- matrix(1:100,ncol=4)
#' condition <- factor(c("A","A","B","B"))
#' dds <- RSeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
#'
#' @rdname RSeqDataSet
#' @importFrom utils packageVersion
#' @export
RSeqDataSet <- function(se, design, ignoreRank=FALSE) {
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }
  if (is.null(assayNames(se)) || assayNames(se)[1] != "counts") {
    message("renaming the first element in assays to 'counts'")
    assayNames(se)[1] <- "counts"
  }
  # special tximeta processing
  if ("tximetaInfo" %in% names(metadata(se))) {
    se <- processTximeta(se)
  }
  # before validity check, try to convert assay to integer mode
  if (any(is.na(assay(se))))
    stop("NA values are not allowed in the count matrix")
  if (any(assay(se) < 0)) {
    stop("some values in assay are negative")
  }
  if (!is.integer(assay(se))) {
    if (!is.numeric(assay(se))) {
      stop(paste("counts matrix should be numeric, currently it has mode:", mode(assay(se))))
    }
    if (any(round(assay(se)) != assay(se))) {
      stop("some values in assay are not integers")
    }
    message("converting counts to integer mode")
    mode(assay(se)) <- "integer"
  }

  if (all(assay(se) == 0)) {
    stop("all samples have 0 counts for all genes. check the counting script.")
  }

  if (all(rowSums(assay(se) == assay(se)[,1]) == ncol(se))) {
    warning("all genes have equal values for all samples. will not be able to perform differential analysis")
  }

  if (any(duplicated(rownames(se)))) {
    warning(sum(duplicated(rownames(se)))," duplicate rownames were renamed by adding numbers")
    rnms <- rownames(se)
    dups <- unique(rnms[duplicated(rnms)])
    for (rn in dups) {
      idx <- which(rnms == rn)
      rnms[idx[-1]] <- paste(rnms[idx[-1]], c(seq_len(length(idx) - 1)), sep= "../..")
    }
    rownames(se) <- rnms
  }

  if (is(design, "formula")) {
    designVars <- all.vars(design)
    if (!all(designVars %in% names(colData(se)))) {
      stop("all variables in design formula must be columns in colData")
    }
    for (v in designVars) {
      if (any(is.na(colData(se)[[v]])))
        stop(paste("variables in design formula cannot contain NA:", v))
    }

    designVarsClass <- sapply(designVars, function(v) class(colData(se)[[v]])[1])

    if (any(designVarsClass == "character")) {
      warning("some variables in design formula are characters, converting to factors")
      for (v in designVars[designVarsClass == "character"]) {
        colData(se)[[v]] <- factor(colData(se)[[v]])
      }
    }

    if (length(designVars) == 1) {
      var <- colData(se)[[designVars]]
      if (all(var == var[1])) {
        stop("design has a single variable, with all samples having the same value.
             use instead a design of '~ 1'. estimateSizeFactors, rlog and the VST can then be used")
      }
      }

    designVarsNumeric <- sapply(designVars, function(v) is.numeric(colData(se)[[v]]))
    if (any(designVarsNumeric)) {
      msgIntVars <- FALSE
      msgCenterScale <- FALSE
      for (v in designVars[designVarsNumeric]) {
        if (all(colData(se)[[v]] == round(colData(se)[[v]]))) {
          msgIntVars <- TRUE
        }
        if (mean(colData(se)[[v]]) > 5 | sd(colData(se)[[v]]) > 5) {
          msgCenterScale <- TRUE
        }
      }
      if (msgIntVars) {
        message("  the design formula contains one or more numeric variables with integer values,
                specifying a model with increasing fold change for higher values.
                did you mean for this to be a factor? if so, first convert
                this variable to a factor using the factor() function")
      }
      if (msgCenterScale) {
        message("  the design formula contains one or more numeric variables that have mean or
                standard deviation larger than 5 (an arbitrary threshold to trigger this message).
                it is generally a good idea to center and scale numeric variables in the design
                to improve GLM convergence.")
      }
      }

    if (any(designVarsClass == "ordered")) {
      stop("the design formula contains an ordered factor. The internal steps
           do not work on ordered factors as a formula. Instead you should provide a matrix to
           the 'design' slot or to the 'full' argument of RSeq(), constructed using model.matrix.")
    }

    designFactors <- designVars[designVarsClass == "factor"]
    missingLevels <- sapply(designFactors, function(v) any(table(colData(se)[[v]]) == 0))
    if (any(missingLevels)) {
      message("factor levels were dropped which had no samples")
      for (v in designFactors[missingLevels]) {
        colData(se)[[v]] <- droplevels(colData(se)[[v]])
      }
    }

    singleLevel <- sapply(designFactors, function(v) all(colData(se)[[v]] == colData(se)[[v]][1]))
    if (any(singleLevel)) {
      stop("design contains one or more variables with all samples having the same value,
           remove these variables from the design")
    }

    # if the last variable in the design formula is a
    # factor, and has a level 'control', check if it is
    # the reference level and if not print a message
    lastDV <- length(designVars)
    if (length(designVars) > 0 && designVarsClass[lastDV] == "factor") {
      lastDVLvls <- levels(colData(se)[[designVars[lastDV]]])
      controlSynonyms <- c("control","Control","CONTROL")
      for (cSyn in controlSynonyms) {
        if (cSyn %in% lastDVLvls) {
          if (cSyn != lastDVLvls[1]) {
            message(paste0("  it appears that the last variable in the design formula, '",designVars[lastDV],"',
                           has a factor level, '",cSyn,"', which is not the reference level. we recommend
                           to use factor(...,levels=...) or relevel() to set this as the reference level
                           before proceeding. for more information, please see the 'Note on factor levels'
                           in vignette('RSeq2')."))
          }
          }
          }
          }

    modelMatrix <- stats::model.matrix.default(design, data=as.data.frame(colData(se)))
          } else if (is(design, "matrix")) {
            modelMatrix <- design
          } else {
            stop("'design' should be a formula or a matrix")
          }

  if (!ignoreRank) {
    checkFullRank(modelMatrix)
  }

  # Add columns on the columns
  mcolsCols <- DataFrame(type=rep("input",ncol(colData(se))),
                         description=rep("",ncol(colData(se))))
  mcols(colData(se)) <- if (is.null(mcols(colData(se)))) {
    mcolsCols
  } else if (all(names(mcols(colData(se))) == c("type","description"))) {
    mcolsCols
  } else {
    cbind(mcols(colData(se)), mcolsCols)
  }

  object <- new("RSeqDataSet", se, design = design)

  # now we know we have at least an empty GRanges or GRangesList for rowRanges
  # so we can create a metadata column 'type' for the mcols
  # and we label any incoming columns as 'input'

  # this is metadata columns on the rows
  mcolsRows <- DataFrame(type=rep("input",ncol(mcols(object))),
                         description=rep("",ncol(mcols(object))))
  mcols(mcols(object)) <- if (is.null(mcols(mcols(object)))) {
    mcolsRows
  } else if (all(names(mcols(mcols(object))) == c("type","description"))) {
    mcolsRows
  } else {
    cbind(mcols(mcols(object)), mcolsRows)
  }

  # stash the package version
  metadata(object)[["version"]] <- packageVersion("RSeq2")

  return(object)
      }

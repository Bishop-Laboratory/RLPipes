#' RSeq Pipeline
#'
#' A flexible pipeline for R-loop mapping
#'
#' @param mode Required. The \code{mode} parameter controls the workflow used to analyze data based
#' on best-practices for the specified R-loop mapping method. This can be overwritten
#' by manually passing parameters to \code{analysis_params}. Options include
#' \code{DRIP}, \code{DRIPc}, \code{qDRIP}, \code{ssDRIP}, \code{RDIP}, \code{R-ChIP},
#' \code{RR-ChIP}, and \code{MapR}. These options can also be specifid in \code{samples}.
#' Please read the vignette for more details.
#'
#' @param experiment Character cotaining the GEO or SRA accession(s) (GSM or SRX) of experimental samples
#' to download and analyze. Can also be the file prefix of local fastq or bam file(s). Must be specified
#' in \code{samples} if not provided here.
#'
#' @param control Same as \code{experiment} but for control/input samples.
#' Provide in the same order as \code{experiment} samples so that there is a 1:1 match
#' between each \code{experiment} sample and the corresponding control. Can be specified in \code{samples}.
#'
#' @param operations Character specifying which analysis steps should be included. Default is \code{'all'} which
#' includes \code{c('fastq_QC', 'alignment', 'post_alignment', 'call_peaks', 'rloop_QC', 'annotate', 'report')}. Can be specified
#' in \code{samples}. See \code{details} below.
#'
#' @param outdir Character specifying folder locations to save data to. Default is \code{'RSeq_out'}. Can
#' be specified in \code{samples}.
#'
#' @param analysis_params Named list specifying additional parameters to pass to each analysis step.
#' Can be specified in \code{samples} instead. See \code{details} below.
#'
#' @param samples Data frame containing sample info and analysis parameters. See \code{details} below.
#'
#' @return An RSeq object containing \code{params}, \code{QC}, \code{report}, and \code{peaks}. See \code{details} below.
#'
#' @details
#'
#' \subsection{The RSeq Workflow}{
#'   To offer flexibility, the \code{RSeq} pipeline allows users to fully control the workflow parameters.
#'   This can be done simply by specifying \code{operations} and \code{analysis_params}
#'
#'   \strong{operations}:
#'
#'   The default pipeline for each \code{mode} follows the same pattern:
#'   \describe{
#'     \item{fastq_QC}{Fastq reads are trimmed and filtered. Performed using \code{fastp} by default, which also provides an HTML report of read quality.}
#'     \item{alignment}{Alignment of reads to genome. Uses \code{bwa mem} with default parameters for alignment. The output is piped to \code{samtools}
#'     for conversion to bam format and sorting/indexing.}
#'     \item{post_alignment}{Duplicate removal and calculation of alignment metrics, alignment specifity,
#'     and sample-sample correlation. Uses \code{samtools} for duplicate removal and alignment metrics. Then,
#'     \code{deepTools} is used to construct a bin matrix and calculate alignment coverage/specificity,
#'      sample-sample correlation, principal component analysis, fragment size distribution. \code{deepTools}
#'      will also be used to calculate peak pileup around the gene body and R-loop forming sequences and
#'      to generate a bigWig coverage track using RPKM as the signal value in 10bp bins.}
#'     \item{call_peaks}{Calculates R-loop peaks using \code{macs2} with the --broad flag by default. If single-end,
#'     fragment size is estimated before peak calling, using macs2 default params}
#'     \item{}
#'   }
#'
#'
#' }
#'
#'
#'
#'
#' \code{RSeq} can be used with or without specifying a \code{samples} file. Users will file that this file
#' provides convenience for complicated analysis. The 'samples' file allows complete flexibility with respect
#' to analysis parameters.
#'
#' \section{Simple Analysis: without 'samples' file}{
#'   e
#' }
#'
#'
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
RSeq <- function(mode = NULL,
                 experiment = NULL,
                 control = NULL,
                 operations = 'all',
                 outdir = 'RSeq_out',
                 analysis_params = NULL,
                 samples = NULL) {

}









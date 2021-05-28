#' Prepare Report
#' Generate final RMarkdown Report
prepare_report <- function(helper_dir, configs, sample_name, output_html, output_data, rlfs_enrichment, bam_stats,
                           homer_annotations, correlation_analysis, read_qc_data, peak_compilation_data, final_peaks) {

  print(c(helper_dir, configs, sample_name, output_html, output_data, rlfs_enrichment, bam_stats,
                           homer_annotations, correlation_analysis, read_qc_data, peak_compilation_data, final_peaks))
  
  # # Bug testing
  # helper_dir <- "/home/millerh1/PycharmProjects/RSeq/bin/../src"
  # configs <- "rseq_out_input/config.json"
  # sample_name <- "SRX1025890_TC32_NT_DRIP"
  # output_html <- "rseq_out_input/RSeq_report/SRX1025890_TC32_NT_DRIP_hg38__RSeq_Report.html"
  # output_data <- "rseq_out_input/RSeq_report/SRX1025890_TC32_NT_DRIP_hg38__RSeq_Report.rda"
  # rlfs_enrichment <- "rseq_out_input/RLFS_analysis/SRX1025890_TC32_NT_DRIP/SRX1025890_TC32_NT_DRIP_hg38__rlfs_enrichment.rda"
  # bam_stats <- "rseq_out_input/bam_stats/SRX1025890_TC32_NT_DRIP/SRX1025890_TC32_NT_DRIP_hg38__bam_stats.txt"
  # homer_annotations <- "rseq_out_input/homer_annotations/SRX1025890_TC32_NT_DRIP/SRX1025890_TC32_NT_DRIP_hg38__feature_overlaps.txt"
  # correlation_analysis <- "rseq_out_input/correlation_analysis/SRX1025890_TC32_NT_DRIP/SRX1025890_TC32_NT_DRIP_hg38__correlation_analysis.rda"
  # read_qc_data <- "rseq_out_input/QC/fastq/json/SRX1025890_TC32_NT_DRIP.hg38.json"
  # peak_compilation_data <- "rseq_out_input/peaks/SRX1025890_TC32_NT_DRIP/SRX1025890_TC32_NT_DRIP_hg38__compiled_peaks.rda"
  # final_peaks <- "rseq_out_input/peaks/SRX1025890_TC32_NT_DRIP/SRX1025890_TC32_NT_DRIP_hg38__compiled_peaks.bed"

  # Pull configs for this sample
  configs <- jsonlite::read_json(configs, simplifyVector = TRUE)
  ind <- which(configs$sample_name == sample_name)
  configlist <- lapply(configs, function(fig) {
    if (length(fig) > 1) {
      fig[ind]
    } else (
      fig
    )
  })
  print(configlist)

  suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))

  # Get correlation results
  load(correlation_analysis)
  corr_data <- list(annoNow = annoNow,
                    corMat = corMat)

  # Get homer annotations
  anno_data <- suppressMessages(read_tsv(homer_annotations)) %>%
      dplyr::filter(! grepl(Annotation, pattern = "\\?"))

  print("AFTER")
  # Parse RLFS data
  load(rlfs_enrichment)
  rlfs_data <- list(pt, lz)

  # Parse peak_compile data
  load(peak_compilation_data)
  final_peaks <- suppressMessages(readr::read_tsv(final_peaks))
  total_peaks <- length(final_peaks$end)

  # Parse bam stats
  bam_stats_raw <- suppressMessages(read_lines(bam_stats))
  bam_stats <- list(
    "reads_aligned" = as.numeric(gsub(bam_stats_raw[1], pattern = "^([0-9]+) \\+ .+", replacement = "\\1")),
    "duplicate_reads" = as.numeric(gsub(bam_stats_raw[4], pattern = "^([0-9]+) \\+ .+", replacement = "\\1"))
  )

  # Parse QC data file
  print(configlist$file_type)
  if (! configlist$file_type %in% c("bam", "peak_coverage")) {
    read_qc_data <- jsonlite::read_json(read_qc_data, simplifyVector = TRUE)
    bam_stats[["total_reads"]] <- read_qc_data$filtering_result$passed_filter_reads
  }

  # TODO: get number of macs2 and epic2 peaks called. 
  # Make markdown report TODO: Needs to include more info...
  md_template <- file.path(helper_dir, "templates/report_template.Rmd")

  # Write files
  data_list <- list(corr_data = corr_data,
                    anno_data = anno_data,
                    rlfs_data = rlfs_data,
                    bam_stats = bam_stats,
                    peak_ol = peak_ol,
                    total_peaks = total_peaks,
                    read_qc_data = read_qc_data,
                    configlist = configlist,
                    helper_dir = helper_dir)
  save(data_list, file = output_data)
  print("Starting render!")

  rmarkdown::render(md_template, 
                    params = data_list, 
                    output_format = "html_document", 
                    output_dir = normalizePath(dirname(output_html)),
                    output_file = output_html)

}

# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)

# Get output JSON
print(arg)
suppressMessages(prepare_report(helper_dir = arg[1],
                                configs = arg[2],
                                sample_name = arg[3],
                                output_html = arg[4],
                                output_data = arg[5],
                                rlfs_enrichment = arg[stringr::str_which(arg, pattern = "__rlfs_enrichment\\.rda$")],
                                bam_stats = arg[stringr::str_which(arg, pattern = "__bam_stats\\.txt$")],
                                homer_annotations = arg[stringr::str_which(arg, pattern = "__feature_overlaps\\.txt$")],
                                correlation_analysis = arg[stringr::str_which(arg, pattern = "__correlation_analysis\\.rda$")],
                                read_qc_data = arg[stringr::str_which(arg, pattern = ".*/QC/fastq/*\\.json$")],
                                peak_compilation_data = arg[stringr::str_which(arg, pattern = "__compiled_peaks\\.rda$")],
                                final_peaks = arg[stringr::str_which(arg, pattern = "__compiled_peaks\\.bed$")]))

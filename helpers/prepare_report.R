#' Prepare Report
#' Generate final RMarkdown Report
prepare_report <- function(input, sample_name, configs) {
  
  # ## Bug testing ##
  # input <- c(
  #   "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/RSeq_out16/SRX113812_Ntera2_DNA/SRX113812_Ntera2_DNA.final_report.tmp.json"
  # )
  # sample_name <- "SRX113812_Ntera2_DNA"
  # configs <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/RSeq_out16/rseqVars.json"
  # ###########
  
  js <- jsonlite::read_json(input, simplifyVector = TRUE)
  
  suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
  
  
  # Parse configs
  configlist <- jsonlite::read_json(configs, simplifyVector = TRUE)
  configlist <- configlist[[sample_name]]
  helpers_dir <- configlist$helpers_dir
  output_html <- paste0(configlist$out_dir, "/", configlist$sample_name, "/", paste0(configlist$sample_name, "_", 
                                                      configlist$genome, ".QC_report.html"))
  
  # Parse cor data
  if (grepl(js$correlation_output[2], pattern = "\\.rda")) {
    load(js$correlation_output[2])
    corr_data <- list(annoNow = annoNow, 
                      corMat = corMat)
  } else {
    corr_data <- NA
  }
  
  # Parse anno data
  if (grepl(js$anno_output[1], pattern = "\\.feature_overlaps\\.txt")) {
    anno_data <- suppressMessages(read_tsv(js$anno_output[1])) %>%
      filter(! grepl(Annotation, pattern = "\\?"))
  } else {
    anno_data <- NA
  }
  
  # Parse QC data file
  read_qc_data <- jsonlite::read_json(js$fastpdata, simplifyVector = TRUE)
  
  # Parse RLFS data
  if (grepl(js$rlfs_output[1], pattern = "\\.rlfs_data\\.txt")) {
    rlfs_data <- suppressMessages(read_tsv(js$rlfs_output[1])) 
  } else {
    rlfs_data <- NA
  }
  
  # Parse cons data
  if (grepl(js$rlcons_output[1], pattern = "\\.rlcons_data\\.txt")) {
    rlcons_data <- suppressMessages(read_tsv(js$rlcons_output[1])) 
  } else {
    rlcons_data <- NA
  }
  
  # Parse bam stats
  bam_stats_raw <- suppressMessages(read_lines(js$bam_stats_output[[1]])) 
  bam_stats <- list(
    "total_reads" = read_qc_data$filtering_result$passed_filter_reads,
    "reads_aligned" = as.numeric(gsub(bam_stats_raw[1], pattern = "^([0-9]+) \\+ .+", replacement = "\\1")),
    "duplicate_reads" = as.numeric(gsub(bam_stats_raw[4], pattern = "^([0-9]+) \\+ .+", replacement = "\\1"))
  )
  
  # TODO: get number of macs2 and epic2 peaks called. 
  
  # Make markdown report TODO: Needs to include more info...
  md_template <- file.path(helpers_dir, "report_template.Rmd")
  output_html <- gsub(output_html, pattern = "//", replacement = "/")
  
  rmarkdown::render(md_template, 
                    params = list(corr_data = corr_data,
                                  anno_data = anno_data,
                                  read_qc_data = read_qc_data,
                                  rlfs_data = rlfs_data,
                                  rlcons_data = rlcons_data,
                                  bam_stats = bam_stats,
                                  configlist = configlist), 
                    output_format = "html_document", 
                    output_dir = normalizePath(dirname(output_html)),
                    output_file = output_html)
  
}

# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)
# print(arg)

# Get output JSON
suppressMessages(prepare_report(input = arg[1],
                                sample_name = arg[2],
                                configs = arg[3]))

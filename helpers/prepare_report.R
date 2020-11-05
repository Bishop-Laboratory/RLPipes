#' Prepare Report
#' Generate final RMarkdown Report
prepare_report <- function(input, sample_name, configs) {
  
  # ## Bug testing ##
  # input <- c("/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/SRX1070678_NT2_DRIP-seq_1/SRX1070678_NT2_DRIP-seq_1.final_report.tmp.json")
  # sample_name <- "SRX1070678_NT2_DRIP-seq_1"
  # configs <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/rseqVars.json"
  # 
  # input <- "data/SRX113812_Ntera2_DNA/SRX113812_Ntera2_DNA.final_report.tmp.json"
  # sample_name <- "SRX113812_Ntera2_DNA"
  # configs <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/rseqVars.json"
  
  # ###########
  input <- gsub(input, pattern = "//", replacement = "/")
  js <- jsonlite::read_json(input, simplifyVector = TRUE)
  
  suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
  
  # Parse configs
  configs <- gsub(configs, pattern = "//", replacement = "/")
  configlist <- jsonlite::read_json(configs, simplifyVector = TRUE)
  configlist <- configlist[[sample_name]]
  helpers_dir <- configlist$helpers_dir
  output_html <- paste0(configlist$out_dir, "/", configlist$sample_name, "/", paste0(configlist$sample_name, "_", 
                                                      configlist$genome, ".QC_report.html"))
  output_rda <- paste0(configlist$out_dir, "/", configlist$sample_name, "/", paste0(configlist$sample_name, "_", 
                                                                                     configlist$genome, ".QC_report.rda"))
  
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
    i <- length(js$anno_output)  # Get stranded if available
    anno_data <- suppressMessages(read_tsv(js$anno_output[i])) %>%
      dplyr::filter(! grepl(Annotation, pattern = "\\?"))
    if(! nrow(anno_data)) {
      anno_data <- NA
    }
  } else {
    anno_data <- NA
  }
  
  # Parse QC data file
  read_qc_data <- jsonlite::read_json(js$fastpdata, simplifyVector = TRUE)
  
  # Parse RLFS data
  if (grepl(js$rlfs_output[1], pattern = "\\.rlfs_enrichment\\.rda")) {
    # i <- length(js$rlfs_output)  # Get stranded if available
    load(js$rlfs_output[1])
    rlfs_data <- list(pt, lz)
  } else {
    rlfs_data <- NA
  }
  
  # Parse cons data
  if (grepl(js$rlcons_output[1], pattern = "\\.rlcons_data\\.txt")) {
    rlcons_data <- suppressMessages(read_tsv(js$rlcons_output[1])) 
  } else {
    rlcons_data <- NA
  }
  
  # Parse peak_compile data
  i <- length(js$peak_compile_output)  # Get stranded if available
  load(js$peak_compile_output[i])
  
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
  output_rda <- gsub(output_rda, pattern = "//", replacement = "/")
  
  data_list <- list(corr_data = corr_data,
                    anno_data = anno_data,
                    read_qc_data = read_qc_data,
                    rlfs_data = rlfs_data,
                    rlcons_data = rlcons_data,
                    bam_stats = bam_stats,
                    configlist = configlist,
                    peak_ol = peak_ol)
  save(data_list, file = output_rda)
  rmarkdown::render(md_template, 
                    params = data_list, 
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

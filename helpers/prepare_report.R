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
  # input <- "data/SRX1761639_UNKNOWN/SRX1761639_UNKNOWN.final_report.tmp.json"
  # sample_name <- "SRX1761639_UNKNOWN"
  
  # input <- "data/SRX8908661_Control_rep_2_Rnase-H_DRIP-seq/SRX8908661_Control_rep_2_Rnase-H_DRIP-seq.final_report.tmp.json"
  # sample_name <- "SRX8908661_Control_rep_2_Rnase-H_DRIP-seq"
  # configs <- "data/rseqVars.json"
  
  # ###########
  
  input <- gsub(input, pattern = "//", replacement = "/")
  js <- jsonlite::read_json(input, simplifyVector = TRUE)
  
  suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
  suppressWarnings(suppressPackageStartupMessages(library(xgboost)))
  
  
  # Parse configs
  configs <- gsub(configs, pattern = "//", replacement = "/")
  configlist <- jsonlite::read_json(configs, simplifyVector = TRUE)
  configlist <- configlist[[sample_name]]
  helpers_dir <- configlist$helpers_dir
  output_html <- paste0(configlist$out_dir, "/", configlist$sample_name, "/", paste0(configlist$sample_name, "_", 
                                                      configlist$genome, ".QC_report.html"))
  output_rda <- paste0(configlist$out_dir, "/", configlist$sample_name, "/", paste0(configlist$sample_name, "_", 
                                                                                     configlist$genome, ".QC_report.rda"))
  tmp_folder <- paste0(configlist$out_dir, "/", configlist$sample_name, "/tmp")
  
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
  
  
  # Add in scores
  show_anno <- ifelse(tibble::is_tibble(anno_data), TRUE, FALSE)
  show_corr <- ifelse(is.list(corr_data), TRUE, FALSE)
  correct_mode <- configlist$mode %in% c("DRIP", "DRIPc", "sDRIP", 'qDRIP')
  if (show_anno & show_corr & correct_mode) {
    rmap_samples <- file.path(helpers_dir, 'data', 'RMapDB_samples_10_22_2020.csv')
    rmap_samples <- read_csv(rmap_samples)
    rmap_samples$mode_group <- ifelse(rmap_samples$mode %in% c("RNH-CnR", "MapR"), "MapR",
                                      ifelse(rmap_samples$mode %in% c("DRIP", "DRIPc", "sDRIP", 'qDRIP'), 'DRIP',
                                             ifelse(rmap_samples$mode %in% c("RDIP"), "RDIP", 
                                                    ifelse(rmap_samples$mode %in% c("ssDRIP"), "ssDRIP", 
                                                           ifelse(rmap_samples$mode %in% c("R-ChIP"), 'R-ChIP', "misc")))))
    # Get corr_median
    my_sample <- rownames(corr_data$annoNow)[which(corr_data$annoNow$Source != "RMapDB")]
    corr_median <- data.frame("clean_name" = rownames(corr_data$annoNow), "R" = corr_data$corMat[my_sample,]) %>% 
      left_join(y = rmap_samples, by = "clean_name") %>%
      dplyr::filter((is.na(mode_group) | mode_group == "DRIP") & 
                      ! Condition %in% c("IgG", "RNaseH1", "RNASEH1", "RNH-high", "WKKD")) %>%
      pull(R) %>% median()
    
    # Get annotations
    annos <- anno_data %>%
      pivot_longer(!Annotation) %>%
      mutate(feature = paste0(Annotation, "__", name)) %>%
      dplyr::select(feature, value) 
    annos <- as.data.frame(t(annos))
    colnames(annos) <- as.character(unlist(annos[1,,drop=TRUE]))
    annos <- annos[-1,]
    
    # Pct aligned  
    total_reads <- bam_stats$total_reads 
    reads_aligned <- bam_stats$reads_aligned
    pct_aligned <- 100*(reads_aligned/total_reads)
    
    # Compile into mat
    test_x <- annos
    test_x$corr_median <- corr_median
    test_x$pct_aligned <- pct_aligned
    test_x <- as.matrix(test_x)
    test_x <- apply(test_x, 1:2, as.numeric)
    
    # XGBoost
    xgb_file <- file.path(helpers_dir, 'data', "xgb_DRIP_group_HS_binary_11_09_2020.model")
    xgb_feats <- c("corr_median", "TTS__Log2 Ratio (obs/exp)", "SINE__Log2 Ratio (obs/exp)",
                   "Exon__Log2 Ratio (obs/exp)", "Intron__Log2 Ratio (obs/exp)")
    if (any(! xgb_feats %in% colnames(test_x))) {
      summary_scores <- NA
    } else {
      test_x_xgb <- test_x[ , xgb_feats, drop = FALSE]
      bst <- xgb.load(modelfile = xgb_file)
      xgbprob <- predict(bst, newdata = test_x_xgb)
      xgbverdict <- ifelse(xgbprob > .5, 'pass', 'fail')
      
      # Linear Regression
      lm_file <- file.path(helpers_dir, 'data', "lm_DRIP_group_HS_nonbinary_11_09_2020.rda")
      test_x <- as.data.frame(test_x)
      colnames(test_x) <- gsub(colnames(test_x), pattern = " |\\(|\\)|\\/", replacement = "_")
      colnames(test_x) <- gsub(colnames(test_x), pattern = "^([0-9]+)", replacement = "d\\1")
      load(lm_file)
      lmscore <- predict(lmfit, newdata = test_x)
      
      # Compile
      summary_scores <- data.frame(xgbprob, xgbverdict, lmscore)
    }
    
  } else {
    summary_scores <- NA
  }
  
  # Write files
  data_list <- list(corr_data = corr_data,
                    anno_data = anno_data,
                    read_qc_data = read_qc_data,
                    rlfs_data = rlfs_data,
                    rlcons_data = rlcons_data,
                    bam_stats = bam_stats,
                    configlist = configlist,
                    peak_ol = peak_ol,
                    summary_scores = summary_scores)
  save(data_list, file = output_rda)
  print("Starting render!")
  rmarkdown::render(md_template, 
                    params = data_list, 
                    output_format = "html_document", 
                    output_dir = normalizePath(dirname(output_html)),
                    output_file = output_html)
  
  
  if (dir.exists(tmp_folder)) {
    print(paste0("Deleting tmp folder: ", tmp_folder))
    system(paste0("rm -rf ", tmp_folder))
  }
  
}

# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)
# print(arg)

# Get output JSON
suppressMessages(prepare_report(input = arg[1],
                                sample_name = arg[2],
                                configs = arg[3]))

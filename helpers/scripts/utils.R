#' Check bam assembly
#'
#' Checks bam file header against chrom.sizes for 'assembly' from UCSC. Returns logical.
check_bam_assembly <- function(bam_file, assembly) {

  # # Bug testing #
  # assembly <- "hg38"
  # bam_file <- "../EUFA_BRCA2//Data/Bam_Files/EUFA_BRCA2_DRIP/EUFA_BRCA2_1_S37_unique_sorted.bam"
  # convert_bams = TRUE

  # Extract chrom sizes
  sizes <- Rsamtools::scanBamHeader(bam_file)[[1]]$targets
  names(sizes) <- gsub(names(sizes), pattern = "chr", replacement = "")
  # Compare with known chr_sizes
  if (! assembly %in% names(chrom_sizes_list) ||
      ! length(chrom_sizes_list[[assembly]])) {
    warning("Was unable to find chrom.sizes for species ", assembly, " and will continue processing anyways.")
    return(TRUE)
  }
  available_chrs <- chrom_sizes_list[[assembly]]
  names(available_chrs) <- gsub(names(available_chrs), pattern = "chr", replacement = "")
  match_ind <- which(names(sizes) %in% names(available_chrs))
  if (! length(match_ind)) {
    # TODO: handle errors arising from species with some contigs in annotation and others in bam
    return(FALSE)
  } else if (! all(sizes[match_ind] %in% available_chrs)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Get bam read length
get_bam_read_length <- function(bam_file) {
  ## Bug testing #
  #bam_file <- "../EUFA_BRCA2//Data/Bam_Files/EUFA_BRCA2_DRIP/EUFA_BRCA2_1_S37_unique_sorted.bam"

  cmd <- paste0("samtools view ", bam_file,
                " | awk '{print length($10)}' | head -10000")
  lens <- as.numeric(system(cmd, intern = TRUE))
  return(mean(lens))
}

#' Get fastq read length
get_fastq_read_length <- function(fastq_file) {
  ## Bug testing #
  #fastq_file <- "~/Bishop.lab/Preprocessing/DRIP_Seq/GSE145964_Developmental_Context_ChIP_DRIP/Data/Raw_Reads/SRR11185284_1.fastq"

  cmd <- paste0("head -n 40000 ", fastq_file,
                " | awk '{if(NR%4==2) print length($1)}'")
  lens <- as.numeric(system(cmd, intern = TRUE))
  return(mean(lens))
}

# Helper function for querying public data accession, convert to SRA, and return run table
get_public_run_info <- function(accessions) {
  suppressWarnings(suppressPackageStartupMessages(require(XML)))


  httr::set_config(httr::config(http_version = 2))
  # set the HTTP version to 1.1 (none, 1.0, 1.1, 2)
  #### Bug testing ##
  #accessions <- c("SRX2481503", "SRX2481504", "GSE134101", "SRP150774", "GSE127329", "SRS1466492")
  #accessions <- c("SRX2918366", "SRX2918367", "GSM3936516", "SRX5129664")
  #accessions <- c("SRX2918366", "SRX2918367", "GSM3936517", "GSM3936517", "GSM3936517", "SRX5129664", "GSM2550995")
  #accessions <- c("SRR2019278")
  # accessions <- public_ctr_accessions
  # accessions <- samples_public$experiment
  # accessions <- public_ctr_accessions
  # accessions <- "SRR3504393"
  # accessions <- public_ctr_accessions
  # accessions <- samples_public$experiment
  ###################


  accessions <- unique(accessions)
  badMsg <- c("HTTP error: Status 429; Too Many Requests",
              "HTTP error: Status 500; Internal Server Error",
              "XML parse error: StartTag: invalid element name\n")

  # Convert GEO series to BioProject accessions
  acc_gse <- accessions[grep(accessions, pattern = "^GSE[0-9]+$")]
  if (length(acc_gse)) {
    fail <- TRUE
    while (fail) {
      esearch_gse <- reutils::esearch(acc_gse, db = "gds")
      if (length(as.character(esearch_gse$errors$error))) {
        if (as.character(esearch_gse$errors$error) %in% badMsg) {
          fail <- TRUE
          Sys.sleep(1)
        }
      } else {
        fail <- FALSE
      }
    }

    if (length(esearch_gse$errors$errmsg)) {
      stop(paste0(paste0(esearch_gse$errors$errmsg, collapse = ", ")), " not found")
    }
    fail <- TRUE
    while (fail) {
      content_bp <- reutils::efetch(esearch_gse)
      if (length(as.character(content_bp$errors$error))) {
        if (as.character(content_bp$errors$error) %in% badMsg) {
          fail <- TRUE
          Sys.sleep(1)
        }
      } else {
        fail <- FALSE
        content_bp <- content_bp$content
      }
    }
    res_bp <- stringr::str_match_all(string = content_bp,
                                     pattern = "/Traces/study/\\?acc=(PRJNA[0-9]+)")[[1]][,2]
    res_gse <- stringr::str_match_all(string = content_bp,
                                      pattern = "geo/series/GSE[0-9]+nnn/(GSE[0-9]+)/")[[1]][,2]
    convert_ord <- order(match(res_gse, acc_gse)) # Need original series order
    map_1 <- data.frame("accessions_original" = accessions,
                        "accessions_cleaned" = NA, stringsAsFactors = FALSE)
    accessions[grep(accessions, pattern = "^GSE[0-9]+$")] <- res_bp[convert_ord]
    map_1$accessions_cleaned <- accessions
  } else {
    map_1 <- data.frame("accessions_original" = accessions,
                        "accessions_cleaned" = accessions, stringsAsFactors = FALSE)
  }

  # Build chunks of 50 accessions
  acc_list <- split(accessions, ceiling(seq_along(accessions)/50))
  res_list_full <- lapply(acc_list, function (accnow) {
    # accnow <- acc_list[[1]]

    # Query all accessions in SRA
    fail <- TRUE
    fail_counter <- 0
    while (fail) {
      Sys.sleep(.5)
      if (fail_counter > 1) {
        stop("Unable to contact NCBI servers. Please use only local files for now and please contact the package maintainer!")
      }
      esearch_sra <- reutils::esearch(accnow, db = "sra")
      if (length(as.character(esearch_sra$errors$error))) {
        if (as.character(esearch_sra$errors$error) %in% badMsg) {
          fail <- TRUE
          if (as.character(esearch_sra$errors$error) == badMsg[2]) {
            fail_counter <- fail_counter + 1
          }
          Sys.sleep(1)
        }
      } else {
        fail <- FALSE
      }
    }
    if (length(esearch_sra$errors$errmsg)) {
      stop(paste0(paste0(esearch_sra$errors$errmsg, collapse = ", ")), " not found")
    }
    fail <- TRUE
    fail_counter <- 0
    while (fail) {
      Sys.sleep(.5)
      if (fail_counter > 1) {
        stop("Unable to contact NCBI servers. Please use only local files for now and please contact the package maintainer!")
      }
      # This is convenient, but might produce an error. 
      # See: https://github.com/gschofl/reutils/issues/14
      content_sra <- suppressMessages(reutils::efetch(esearch_sra, retmode = "xml"))
      if (length(as.character(content_sra$errors$error))) {
        if (as.character(content_sra$errors$error) %in% badMsg) {
          fail <- TRUE
          if (as.character(content_sra$errors$error) %in% badMsg[1:2]) {
            # CASE: Error due to HTTP issue. 
            # Could be because too many requests or temporary gateway down.
            # Will continue trying 2 more times + sys.sleep().
            Sys.sleep(1)
            fail_counter <- fail_counter + 1
          } else if (as.character(content_sra$errors$error) == badMsg[3]) {
            # CASE: XML parsing error. 
            # See: https://github.com/gschofl/reutils/issues/14
            ct <- XML::xmlToList(esearch_sra$content)
            ids <- unlist(ct$IdList, use.names = FALSE)
            url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=",
                          paste0(ids, collapse = ","))
            content_sra <- XML::xmlParse(httr::content(httr::GET(url), "text"))
            fail <- FALSE
          }
        }
      } else {
        fail <- FALSE
        content_sra <- reutils::content(content_sra)
      }
    }
    result <- XML::xmlToList(content_sra)

    # Need to map back to original entries -- very annoying...
    # TODO: Figure out a better way to do this without user error issues
    fltr <- reutils::make_flattener()
    acc_matchs <- lapply(result, function(res_now) {
      flt_res <- unlist(fltr(res_now))
      acc_match <- accnow[which(accnow %in% flt_res)]
      names(acc_match) <- res_now$EXPERIMENT$IDENTIFIERS$PRIMARY_ID
      if (length(acc_match) > 1) {
        stop("Multiple provided accessions (",paste0(acc_match, collapse = ", "),
             ") mapped to identical SRA entries: ",
             res_now$EXPERIMENT$IDENTIFIERS$PRIMARY_ID, ". This interferes with assigning",
             " parameters for experimental conditions. Please remove these duplcates.")
      }
      return(acc_match)
    })
    acc_matchs <- unlist(acc_matchs, use.names = TRUE)
    names(acc_matchs) <- gsub(names(acc_matchs), pattern = "EXPERIMENT_PACKAGE\\.", replacement = "")

    map_2 <- data.frame("accessions_cleaned" = acc_matchs,
                        "sra_experiment" = names(acc_matchs), stringsAsFactors = FALSE)
    map_merge <- merge(x = map_1, y = map_2, by = "accessions_cleaned", all = TRUE)

    # Unpack experiment list
    resList <- lapply(result, FUN = function(exp_now) {
      # exp_now <- result[[1]]

      SRX <- exp_now$EXPERIMENT$IDENTIFIERS$PRIMARY_ID
      SRRs <- unlist(exp_now$RUN_SET, use.names = FALSE)
      SRRs <- unique(SRRs[grep(SRRs, pattern = "^[SE]RR[0-9]+$")])
      new_exp_name <- paste0(SRRs, collapse = ",")
      # TODO: see if there's some way to get condition info from SRA
      condition <- exp_now$STUDY$IDENTIFIERS$PRIMARY_ID
      out_name <- exp_now$SAMPLE$TITLE
      if (is.null(out_name)) {
        out_name <- "UNKNOWN"
      }
      sample_name <- paste0(SRX, "_", clean_str(out_name))
      sample_tx <- exp_now$SAMPLE$SAMPLE_NAME$TAXON_ID
      possible_genomes <- available_genomes[available_genomes$taxId == sample_tx,]
      if (sample_tx %in% c("4932", "559292")) {
        # Special processing for yeast genome... which is incorrectly identified
        possible_genomes <- available_genomes[available_genomes$taxId %in% c("4932", "559292"),]
      }
      lib_type <- names(exp_now$EXPERIMENT$DESIGN$LIBRARY_DESCRIPTOR$LIBRARY_LAYOUT)
      if ("PAIRED" %in% lib_type) {
        paired_end <- TRUE
      } else {
        paired_end <- FALSE
      }
      if (! length(possible_genomes$UCSC_orgID)) {
        stop("No UCSC genome assembly available for sample ", SRX, " with taxonomy ID: ", sample_tx)
      }
      genomes_available <- possible_genomes[which(possible_genomes$genes_available),]
      if (! length(genomes_available$UCSC_orgID)) {
        warning("No UCSC genome annotations available for sample ", SRX, " with taxonomy ID: ", sample_tx,
                ". Will not run steps which require annotations.")
        genome <- possible_genomes$UCSC_orgID[which.max(possible_genomes$year)]
      } else {
        genome <- genomes_available$UCSC_orgID[which.max(genomes_available$year)]
      }
      read_length <- as.numeric(exp_now$RUN_SET$RUN$Statistics$Read["average"])
      res_now <- data.frame(
        sra_experiment = SRX,
        experiment = new_exp_name,
        genome = genome,
        condition = condition,
        paired_end = paired_end,
        sample_name = sample_name,
        out_name = out_name,
        read_length = read_length
      )
      return(res_now)
    })

    res_df <- dplyr::bind_rows(resList)
    map_final <- merge(x = map_merge, y = res_df, by = "sra_experiment", all = TRUE)
  })

  res_df_full <- dplyr::bind_rows(res_list_full)
  res_df_full <- res_df_full[! is.na(res_df_full$experiment),]
  return(res_df_full)
}


clean_str <- function(str) {
  # str <- "24hr input S 9.6 DRIP Seq"

  # Unsafe regex patterns to replace
  repList <- list(" " = "_",
                  "\\+|\\&" = "and",
                  '\\"|\\{|\\}|\\[|\\]|\\(|\\)|<|>|,|\\=|\\`|\\`' = "",
                  "\\||:|;|'|\\.|~|!|@|\\?|#|\\$|%|\\^|\\*|/|\\\\" = "")

  for (repNow in names(repList)) {
    toRep <- repNow
    repWith <- repList[[repNow]]
    str <- gsub(str, pattern = toRep, replacement = repWith)
  }
  return(str)
}

# Get genome sizes using unique-kmers.py script from khmer
get_genome_sizes <- function() {

  genome_size_list <- list()
  for (i in seq(available_genomes$UCSC_orgID)) {
    genome_now <- available_genomes$UCSC_orgID[i]
    print(genome_now)
    fasta_file <- paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/",
                         genome_now, "/bigZips/", genome_now, ".fa.gz")
    dir.create(file.path(system.file(package = "../../RSeq"),
                         "../extra/genomes", genome_now), recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(system.file(package = "../../RSeq"),
                          "../extra/genomes", genome_now,
                          paste0(genome_now, ".fa.gz"))
    if (! file.exists(out_file) & ! file.exists(paste0(out_file, "_300.txt"))) {

      err <- tryCatch(download.file(fasta_file, destfile = out_file), error = function(e) {
        return(e)
      })
      if ("error" %in% class(err)) {
        print(paste0(genome_now, ' link is broken...'))
        next
      }
    } else {
      print("Info already found or fasta already found!")
    }

    lengths <- c(36, 50, 75, 100, 125, 150, 200, 250, 300)
    sizes <- character(length = length(lengths))
    for (length_now in lengths) {
      cmd <- paste0("unique-kmers.py -k ", length_now, " -R ",
                    out_file, "_", length_now, ".txt ", out_file)
      if (! file.exists(paste0(out_file, "_", length_now, ".txt"))) {
        system(cmd)
      } else {
        print(paste0("Already ran length: ", length_now))
      }
      lines <- readLines(paste0(out_file, "_", length_now, ".txt"))
      lines <- lines[grep(x = lines,"number of unique k-mers: ")]
      size_now <- gsub(lines, pattern = "number of unique k-mers: \t([0-9]+)$",
                       replacement = "\\1")
      sizes <- c(sizes, size_now)
    }
    sizes <- as.numeric(sizes)
    genome_size_list_df <- data.frame(
      "UCSC_orgID" = genome_now,
      "eff_genome_size_36bp" = sizes[1],
      "eff_genome_size_50bp" = sizes[2],
      "eff_genome_size_75bp" = sizes[3],
      "eff_genome_size_100bp" = sizes[4],
      "eff_genome_size_125bp" = sizes[5],
      "eff_genome_size_150bp" = sizes[6],
      "eff_genome_size_200bp" = sizes[7],
      "eff_genome_size_250bp" = sizes[8],
      "eff_genome_size_300bp" = sizes[9]
    )
    genome_size_list[[genome_now]] <- genome_size_list_df
    if (file.exists(out_file)) {
      file.remove(out_file)
    }
  }
  genome_sizes <- dplyr::bind_rows(genome_size_list)
  save(genome_sizes, file = "genome_sizes.rda")
  available_genomes <- available_genomes
  available_genome_info <- merge(x = available_genomes, y = genome_sizes, by = "UCSC_orgID")
  # load("full_len_list.rda")
  #available_genomes <- available_genomes
  #len_data_frame <- data.frame("UCSC_orgID" = names(full_len_list),
  #                             "genome_length" = unlist(full_len_list))
  #available_genomes <- merge(x = available_genomes, y = len_data_frame, by = "UCSC_orgID")
  # available_genomes <- available_genome_info
  # usethis::use_data(available_genomes, overwrite = TRUE)
  return(available_genome_info)
}



#' Get chrome sizes
#' Helper function to return list with all available chrom.sizes.
get_chrom_sizes <- function() {

  # from http://rstudio-pubs-static.s3.amazonaws.com/562103_092b7264f392482e827440cf1521363c.html
  api_genome <- restfulr::RestUri("http://api.genome.ucsc.edu/list")
  response <- restfulr::read(api_genome$ucscGenomes)
  available_genomes <- do.call(rbind.data.frame, response$ucscGenomes)
  available_genomes <- cbind(rownames(available_genomes), available_genomes)
  colnames(available_genomes)[1] <- "UCSC_orgID"
  # Check to see if genome predictions are available
  chrom_sizes_list <- list()
  full_len_list <- list()
  for (i in seq(available_genomes$UCSC_orgID)) {
    genome_now <- available_genomes$UCSC_orgID[i]
    print(as.character(genome_now))

    chr_now <- tryCatch(read.table(paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/",
                                          genome_now, "/bigZips/", genome_now, ".chrom.sizes"),
                                   stringsAsFactors = FALSE),
                        error = function(e) e)
    if ("error" %in% class(chr_now)) {next}

    print(length(chr_now$V1))
    full_len_list[[as.character(genome_now)]] <- sum(chr_now$V2)
    chr_now <- chr_now[grep(chr_now$V1, pattern = "chrUn|random|_alt|_hap|Het$", invert = TRUE),]
    print(length(chr_now$V1))
    if (length(chr_now$V1) > 200) {n <- 200} else {n <- length(chr_now$V1)}
    chrom_sizes <- chr_now$V2[1:n]
    names(chrom_sizes) <- chr_now$V1[1:n]
    chrom_sizes_list[[as.character(genome_now)]] <- chrom_sizes
  }

  # # For internal use
  # usethis::use_data(chrom_sizes_list)
  save(full_len_list, file = "full_len_list.rda")
  return(chrom_sizes_list)
}

#' Get gold standard genes
#' Helper function to get the gold standard genes in all different species
get_gs_genes <- function() {

  hg19gs <- import("RSeq_CLI/helpers/data/correlation_genes_100kb.bed")

  require(rtracklayer)

  # hg38
  genome_now <- "Hg38"
  tmp <- tempfile()
  tmpgz <- paste0(tmp, ".gz")
  download.file(url = paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19To", genome_now, ".over.chain.gz"),
                destfile = tmpgz)
  R.utils::gunzip(tmpgz)
  chain <- import.chain(tmp)
  hg38gs <- unique(unlist(liftOver(chain = chain, x = hg19gs)))

  gs_list <- list("hg19" = hg19gs,
       "hg38" = hg38gs)
  save(gs_list, file = "RSeq_CLI/helpers/data/gs_list.rda")

}

# Check to make sure homer annotations available
check_homer_anno <- function(available_genomes, genome = NULL) {
  if (is.null(genome) || ! "homer_anno_available" %in% colnames(available_genomes)) {
    # Check all possible genomes
    res <- sapply(available_genomes$UCSC_orgID, function(genome_now) {
      urlnow <- paste0("http://homer.ucsd.edu/homer/data/genomes/", genome_now, ".v6.4.zip")
      RCurl::url.exists(urlnow)
    })
    available_genomes$homer_anno_available <- res
    if (all(! res)) {
      warning("Unable to locate any Homer Annotations. If found, please report to package maintainer!!")
    }
  } else {
    urlnow <- paste0("http://homer.ucsd.edu/homer/data/genomes/", genome, ".v6.4.zip")
    res <- RCurl::url.exists(urlnow)
    available_genomes$homer_anno_available[which(available_genomes$UCSC_orgID == genome)] <- res
  }
  
  return(available_genomes)
}


# Get RLFSs
get_rlfs <- function() {
  helpers_dir <- paste0(path.expand("~"), "/Bishop.lab/Projects/RSeq/helpers/")
  
  load(file.path(helpers_dir, "../data", "available_genomes.rda"))
  script <- file.path(helpers_dir, "../external", "QmRLFS-finder.py")
  outdir <- file.path(helpers_dir, "../data", "rlfs")
  dir.create(outdir, showWarnings = FALSE)
  
  genomes <- available_genomes$UCSC_orgID[available_genomes$homer_anno_available]
  
  for (genome_now in genomes) {
    print(genome_now)
    if (! file.exists(file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt"))) && 
        ! file.exists(file.path(outdir, paste0(genome_now, ".fa")))) {
      print("Getting genome!")
      download.file(paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/", genome_now, "/bigZips/", genome_now, ".fa.gz"),
                    destfile = file.path(outdir, paste0(genome_now, ".fa.gz")))
      R.utils::gunzip(file.path(outdir, paste0(genome_now, ".fa.gz")))
    }
    cmd <- paste0("python ", script, " -i ", file.path(outdir, paste0(genome_now, ".fa")),
                  " -o ", file.path(outdir, paste0(genome_now, ".rlfs")))
    if (! file.exists(file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt")))) {
      print("Calculating RLFS")
      system(cmd, wait = FALSE)
    }
    
    # Conver to bed6 and liftOver
    
    cmd <- paste0("awk 'FNR > 1 {print $3\"\\t\"$21}' ", file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt")),
                  "| awk '{gsub(\":|-\",\"\\t\", $1); print $1\"\\t\"\".\"\"\\t\"\".\"\"\\t\"$2}'",
                  " | bedtools sort -i stdin | mergeBed -i stdin -s -c 4,5,6 -o distinct > ", 
                  file.path(outdir, paste0(genome_now, ".rlfs.bed")))
    if (file.exists(file.path(outdir, paste0(genome_now, ".rlfs.out.table.txt"))) &&
        ! file.exists(file.path(outdir, paste0(genome_now, ".rlfs.bed")))) {
      print("Convert to bed!")
      system(cmd)
    }
  }
  
  # Clean up
  bed_files <- list.files(outdir, pattern = "\\.bed", full.names = F)
  lapply(bed_files, function(filenow) {
    old_file <- file.path(outdir, filenow)
    new_file <- file.path(paste0(outdir, '-beds'), filenow)
    file.copy(old_file, new_file, overwrite = TRUE)
  })
  
}
# # Collate data
# collate_rda <- function() {
#   helpers_dir <- paste0(path.expand("~"), "/Bishop.lab/Projects/RSeq/helpers/")
#   outdir <- file.path(helpers_dir, "export")
#   dir.create(outdir, showWarnings = FALSE)
#   outfile <- file.path(outdir, "RMapDB.h5")
#   
#   qc_data <- list.files(path = file.path(dirname(helpers_dir), "tests"), full.names = TRUE,
#                         pattern = "QC_report.rda", recursive = TRUE)
#   
#   library(rhdf5)
#   library(tidyverse)
#   library(readxl)
#   
#   if (! file.exists(outfile)) {
#     h5createFile(outfile)
#   }
#   
#   h5createGroup(file = outfile, "meta")
#   rmap_samples <- read_csv(file.path(dirname(helpers_dir), "misc/analysis/rmapsamples_10_05_2020.csv"))
#   
#   
#   for (i in 1:length(qc_data)) {
#     dnow <- qc_data[i]
#     load(dnow)
#     
#   }
#   
#   
# }



# Fix MACS2 and EPIC2 filenames -- needed for switch to new naming schema 11/5/2020
fix_peak_files <- function() {
  # Macs unstranded
  bp <- list.files("~/Bishop.lab/Projects/RMapDB/data/", recursive = TRUE, full.names = TRUE, 
                   pattern = ".+_peaks.broadPeak")
  bp <- bp[grep(bp, pattern = "peaks_macs_unstranded")]
  newbp <- gsub(bp, pattern = "_peaks(\\.broadPeak)", replacement = ".unstranded\\1")
  lapply(seq(length(bp)), function(i) {file.rename(bp[i], newbp[i])})
  
  # Macs stranded
  bp <- list.files("~/Bishop.lab/Projects/RMapDB/data/", recursive = TRUE, full.names = TRUE, 
                   pattern = ".+_[plus|minus]+_peaks.broadPeak")
  bp <- bp[grep(bp, pattern = "peaks_macs_stranded")]
  newbp <- gsub(bp, pattern = "_([plus|minus]+_peaks\\.broadPeak)", replacement = ".\\1")
  lapply(seq(length(bp)), function(i) {file.rename(bp[i], newbp[i])})
  
  # Epic unstranded
  bp <- list.files("~/Bishop.lab/Projects/RMapDB/data/", recursive = TRUE, full.names = TRUE, 
                   pattern = "_[a-zA-Z]+[0-9]+\\.bed")
  bp <- bp[grep(bp, pattern = "peaks_epic_unstranded")]
  newbp <- gsub(bp, pattern = "\\.bed", replacement = ".unstranded.bed")
  lapply(seq(length(bp)), function(i) {file.rename(bp[i], newbp[i])})
  
  # Epic stranded
  bp <- list.files("~/Bishop.lab/Projects/RMapDB/data/", recursive = TRUE, full.names = TRUE, 
                   pattern = ".+_[plus|minus]+\\.bed")
  bp <- bp[grep(bp, pattern = "peaks_epic_stranded")]
  newbp <- gsub(bp, pattern = "_([plus|minus]+\\.bed)", replacement = ".\\1")
  lapply(seq(length(bp)), function(i) {file.rename(bp[i], newbp[i])})
  
  
}

  
# Send peaks to AWS
fix_peak_files <- function() {
  # Peaks final unstranded
  bp <- list.files("~/Bishop.lab/Projects/RMapDB/data/", recursive = TRUE, full.names = TRUE, 
                   pattern = "_[a-zA-Z]+[0-9]+\\.unstranded\\.bed")
  bp <- bp[grep(bp, pattern = "peaks_final_unstranded")]
  
  # Copy files to new folder
  dir.create("helpers/export/final-peaks-unstranded")
  
  for (i in 1:length(bp)) {
    print(i)
    filenow <- bp[i]
    file.copy(filenow, to = "helpers/export/final-peaks-unstranded")
    
  }
  
}

# Send quality reports to AWS
prepare_reports_to_upload <- function() {
  # Peaks final unstranded
  bp <- list.files("~/Bishop.lab/Projects/RMapDB/data/", recursive = TRUE, full.names = TRUE, 
                   pattern = "._report\\.html")
  
  # Copy files to new folder
  dir.create("helpers/export/rseq-quality-reports")
  
  for (i in 1:length(bp)) {
    print(i)
    filenow <- bp[i]
    file.copy(filenow, to = "helpers/export/rseq-quality-reports")
    
  }
  
}


prepare_unstranded_bw_to_upload <- function() {
  # Peaks final unstranded
  bp <- list.files("~/Bishop.lab/Projects/RMapDB/data/", recursive = TRUE, full.names = TRUE, 
                   pattern = "\\.bw")
  bp <- bp[grep(bp, pattern = "coverage_unstranded")]
  
  # Copy files to new folder
  dir.create("helpers/export/rseq-coverage-unstranded")
  
  for (i in 1:length(bp)) {
    print(i)
    filenow <- bp[i]
    file.copy(filenow, to = "helpers/export/rseq-coverage-unstranded")
    
  }
  
}

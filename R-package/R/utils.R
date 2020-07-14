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
  if (! assembly %in% names(RSeq::chrom_sizes_list) ||
      ! length(RSeq::chrom_sizes_list[[assembly]])) {
    warning("Was unable to find chrom.sizes for species ", assembly, " and will continue processing anyways.")
    return(TRUE)
  }
  available_chrs <- RSeq::chrom_sizes_list[[assembly]]
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


#' Get available genomes
#' Helper function to return data frame with info on available genomes from UCSC.
get_available_genomes <- function() {

  # from http://rstudio-pubs-static.s3.amazonaws.com/562103_092b7264f392482e827440cf1521363c.html
  api_genome <- restfulr::RestUri("http://api.genome.ucsc.edu/list")
  response <- restfulr::read(api_genome$ucscGenomes)
  available_genomes <- do.call(rbind.data.frame, response$ucscGenomes)
  available_genomes <- cbind(rownames(available_genomes), available_genomes)
  colnames(available_genomes)[1] <- "UCSC_orgID"
  # Check to see if genome predictions are available
  resList <- list()
  for (i in 1:length(available_genomes$UCSC_orgID)) {
    genome_now <- available_genomes$UCSC_orgID[i]
    # print(i)
    resList[[i]] <- tryCatch(! httr::http_error(paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/",
                                                       genome_now, "/bigZips/genes/README.txt")),
                             error = function(e) {
                               if (class(e) == "simpleError") {
                                 FALSE
                               } else {
                                 stop(e)
                               }})

  }

  available_genomes$genes_available <- unlist(resList)
  available_genomes
  available_genomes <- data.frame(lapply(available_genomes, as.character), stringsAsFactors=FALSE)
  available_genomes$taxId <- as.numeric(available_genomes$taxId)
  available_genomes$genes_available <- as.logical(available_genomes$genes_available)
  years <- gsub(available_genomes$description, pattern = ".+ (20[0-9][0-9])[ /].+", replacement = "\\1")
  available_genomes$year <- as.numeric(years)
  # # For internal use
  # usethis::use_data(available_genomes, overwrite = TRUE)
  return(available_genomes)
}



# Helper function for querying public data accession, convert to SRA, and return run table
get_public_run_info <- function(accessions) {

   ### Bug testing ##
   #accessions <- c("SRX2481503", "SRX2481504", "GSE134101", "SRP150774", "GSE127329", "SRS1466492")
   #accessions <- c("SRX2918366", "SRX2918367", "GSM3936516", "SRX5129664")
   #accessions <- c("SRX2918366", "SRX2918367", "GSM3936517", "GSM3936517", "GSM3936517", "SRX5129664", "GSM2550995")
   ##################

  accessions <- unique(accessions)

  # Convert GEO series to BioProject accessions
  acc_gse <- accessions[grep(accessions, pattern = "^GSE[0-9]+$")]
  if (length(acc_gse)) {
    fail = TRUE
    while (fail) {
      esearch_gse <- reutils::esearch(acc_gse, db = "gds")
      if (length(as.character(esearch_gse$errors$error))) {
        if (as.character(esearch_gse$errors$error) == "HTTP error: Status 429; Too Many Requests") {
          fail = TRUE
          Sys.sleep(1)
        }
      } else {
        fail = FALSE
      }
    }


    if (length(esearch_gse$errors$errmsg)) {
      stop(paste0(paste0(esearch_gse$errors$errmsg, collapse = ", ")), " not found")
    }
    fail = TRUE
    while (fail) {
      content_bp <- reutils::efetch(esearch_gse)
      if (length(as.character(content_bp$errors$error))) {
        if (as.character(content_bp$errors$error) == "HTTP error: Status 429; Too Many Requests") {
          fail = TRUE
          Sys.sleep(1)
        }
      } else {
        fail = FALSE
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


  # Query all accessions in SRA
  fail = TRUE
  while (fail) {
    Sys.sleep(.5)
    esearch_sra <- reutils::esearch(accessions, db = "sra")
    if (length(as.character(esearch_sra$errors$error))) {
      if (as.character(esearch_sra$errors$error) == "HTTP error: Status 429; Too Many Requests") {
        fail = TRUE
        Sys.sleep(1)
      }
    } else {
      fail = FALSE
    }
  }
  if (length(esearch_sra$errors$errmsg)) {
    stop(paste0(paste0(esearch_sra$errors$errmsg, collapse = ", ")), " not found")
  }
  fail = TRUE
  while (fail) {
    Sys.sleep(.5)
    content_sra <- reutils::efetch(esearch_sra)
    if (length(as.character(content_sra$errors$error))) {
      if (as.character(content_sra$errors$error) == "HTTP error: Status 429; Too Many Requests") {
        fail = TRUE
        Sys.sleep(1)
      }
    } else {
      fail = FALSE
      content_sra <- content_sra$content
    }
  }

  result <- XML::xmlToList(XML::xmlParse(content_sra))

  # Need to map back to original entries -- very annoying...
  # TODO: Figure out a better way to do this without user error issues
  fltr <- reutils::make_flattener()
  acc_matchs <- lapply(result, function(res_now) {
    flt_res <- unlist(fltr(res_now))
    acc_match <- accessions[which(accessions %in% flt_res)]
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
    # exp_now <- result[[30]]
    SRX <- exp_now$EXPERIMENT$IDENTIFIERS$PRIMARY_ID
    SRRs <- unlist(exp_now$RUN_SET, use.names = FALSE)
    SRRs <- unique(SRRs[grep(SRRs, pattern = "^SRR[0-9]+$")])
    new_exp_name <- paste0(SRX, "+", paste0(SRRs, collapse = ","))
    # TODO: see if there's some way to get condition info from SRA
    condition <- exp_now$STUDY$IDENTIFIERS$PRIMARY_ID
    out_name <- exp_now$SAMPLE$TITLE
    sample_name <- paste0(SRX, "_", clean_str(out_name))
    sample_tx <- exp_now$SAMPLE$SAMPLE_NAME$TAXON_ID
    possible_genomes <- RSeq::available_genomes[RSeq::available_genomes$taxId == sample_tx,]
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
      genes_available <- FALSE
      genome <- possible_genomes$UCSC_orgID[which.max(possible_genomes$year)]
    } else {
      genes_available <- TRUE
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

  res_df <- data.table::rbindlist(resList)
  map_final <- merge(x = map_merge, y = res_df, by = "sra_experiment", all = TRUE)

  return(map_final)
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
  available_genomes <- RSeq::available_genomes

  genome_size_list <- list()
  for (i in 1:length(available_genomes$UCSC_orgID)) {
    genome_now <- available_genomes$UCSC_orgID[i]
    print(genome_now)
    fasta_file <- paste0("ftp://hgdownload.soe.ucsc.edu/goldenPath/",
                         genome_now, "/bigZips/", genome_now, ".fa.gz")
    dir.create(file.path(system.file(package = "RSeq"),
                         "../extra/genomes", genome_now), recursive = TRUE, showWarnings = FALSE)
    out_file <- file.path(system.file(package = "RSeq"),
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
    sizes <- c()
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
  genome_sizes <- data.table::rbindlist(genome_size_list)
  save(genome_sizes, file = "genome_sizes.rda")
  available_genomes <- RSeq::available_genomes
  available_genome_info <- merge(x = available_genomes, y = genome_sizes, by = "UCSC_orgID")
  # load("full_len_list.rda")
  #available_genomes <- RSeq::available_genomes
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
  for (i in 1:length(available_genomes$UCSC_orgID)) {
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





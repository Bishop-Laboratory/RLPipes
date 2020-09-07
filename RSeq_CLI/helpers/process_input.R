processInput <- function(mode = NULL,
                         outdir = "RSeq_out/",
                         genome = NULL,
                         genome_home_dir = file.path(path.expand("~"), ".RSeq_genomes"),
                         cores = NULL,
                         no_dedupe = FALSE,
                         no_fastp = FALSE,
                         samples = NULL) {


  ### For bug testing ##
  #mode = arg[-1][1]
  #outdir = arg[-1][2]
  #genome = arg[-1][3]
  #genome_home_dir = arg[-1][4]
  #cores = arg[-1][5]
  #samples <- arg[-1][6]
  ## #####################

  # Additional data
  # TODO: add support for bigWig and bedGraph
  read_align_pattern <- "\\.bam$|\\.f[ast]*q$"
  public_pattern <- paste0("^SR[RAXSP][0-9]+$|^GS[EM][0-9]+$|^PRJNA[0-9]+$|^SAMN[0-9]+$")
  R1_pattern <- "[\\._]+[rR]*1\\.f[ast]*q$"
  R2_pattern <- "[\\._]+[rR]*2\\.f[ast]*q$"
  fq_end_pattern <- paste0(R1_pattern, "|", R2_pattern, "|[\\._]+[rR]*\\*.f[ast]*q$|\\.f[ast]*q$")

  mode_df <- data.frame(
  mode = c("DRIP", "ssDRIP", "qDRIP",
           "sDRIP", "S1-DRIP", "bisDRIP",
           "SMRF", "RDIP", "DRIPc",
           "R-ChIP", "RR-ChIP", "DRIVE",
           "MapR", "RNH-CnR"),
  strand_specific = c(FALSE, rep(TRUE, 3), FALSE, rep(TRUE, 6), rep(FALSE, 3)),
  ip_type = c(rep("S9.6", 9), rep("RNaseH1", 5)),
  moeity = c(rep("DNA", 7), rep("RNA", 7)), stringsAsFactors = FALSE
  )
  modes <- c(mode_df$mode, "custom")

  # Parse arguments
  tryCatch(stopifnot(! is.null(samples) || ! is.null(experiment)),
         error = function(x) stop("Must provide 'experiment' or 'samples'"))
  if (! is.null(samples)) {
    # Check if samples is a csv file
    if (typeof(samples) == "character") {
      if (file.exists(samples)) {
        samples <- read.csv(samples, as.is = TRUE)
      } else {
        stop("Samples is provided as a character, but no such file exists: ", samples, ". Supply either a valid",
             " file name or a data frame as the sample sheet.")
      }
    }

    # Remove any empty rows and columns
    samplesbin <- as.data.frame(apply(samples, 1:2, function(x) {
      if (x == "" || is.na(x)) {
        x <- 0
      } else {
        x <- 1
      }
      return(x)
    }))
    rm_cols <- which(colSums(samplesbin) == 0 | colnames(samplesbin) == "X")
    rm_rows <- which(rowSums(samplesbin) == 0)
    if (length(rm_cols)) {
      samples <- samples[, -rm_cols, drop=FALSE]
    }
    if (length(rm_rows)) {
      samples <- samples[-rm_rows, , drop=FALSE]
    }

    # Replace any empty cells with NA
    samples <- as.data.frame(apply(samples, 1:2, function(x) {
      if (x == "" || is.na(x)) {
        x <- NA
      }
      return(x)
    }))

    # Replace factor with char
    i <- sapply(samples, is.factor)
    samples[i] <- lapply(samples[i], as.character)
    if (! "experiment" %in% colnames(samples)) {
      stop("Column 'experiment' is required but not provided in sample sheet.")
    }
  }

  # Fix genome home dir if using tilde
  genome_home_dir <- gsub(genome_home_dir, pattern = "~", replacement = path.expand("~"))

  # Check for arguments. If not set, replace with defaults
  if (! "cores" %in% colnames(samples) &
    is.null(cores)) samples$cores <- parallel::detectCores()/2
  if (! "control" %in% colnames(samples)) {samples$control <- NA}
  if (! "paired_end" %in% colnames(samples)) {samples$paired_end <- NA}
  if (! "effective_genome_size" %in% colnames(samples)) {samples$effective_genome_size <- NA}
  if (! "full_genome_length" %in% colnames(samples)) {samples$full_genome_length <- NA}
  if (! "genome_home_dir" %in% colnames(samples)) {samples$genome_home_dir <- genome_home_dir}
  if (! "genome" %in% colnames(samples)) {samples$genome <- NA}

  # Deal with strand specificity and mode type
  if (! is.null(mode) & ! "mode" %in% colnames(samples)) {
    samples$mode <- mode
  } else if (! "mode" %in% colnames(samples) & is.null(mode)) {
    stop("Mode was not specified in arguments or in sample sheet. Please see usage for more details.")
  }
  if (! all(samples$mode %in% modes)) {
    bad_mode <- unique(samples$mode[! samples$mode %in% modes])
    stop(paste0("Mode(s) specified in sample sheet are not available: ", paste0(bad_mode, collapse = ", ")))
  }
  if (! "strand_specific" %in% colnames(samples)) {samples$strand_specific <- NA}
  if (! "moeity" %in% colnames(samples)) {samples$moeity <- NA}
  if (! "ip_type" %in% colnames(samples)) {samples$ip_type <- NA}
  for (i in 1:length(samples$experiment)) {
    mode_now <- samples$mode[i]
    if(is.na(samples$strand_specific[i])) {samples$strand_specific[i] <- mode_df$strand_specific[mode_df$mode == mode_now]}
    if(is.na(samples$moeity[i])) {samples$moeity[i] <- mode_df$moeity[mode_df$mode == mode_now]}
    if(is.na(samples$ip_type[i])) {samples$ip_type[i] <- mode_df$ip_type[mode_df$mode == mode_now]}
  }
  if (! "cores" %in% colnames(samples)) {samples$cores <- cores}

  # Add in genome if provided
  if (! is.null(genome)) {
    if (! genome %in% available_genomes$UCSC_orgID) {
      stop("Specified genome, ", genome, ", is not a valid UCSC genome ID. Please see the manual for more information.")
    } else {
      samples$genome <- genome
    }
  }

  # Mark files, and public samples
  fileIndExp <- grep(samples$experiment, pattern = read_align_pattern)
  publicIndExp <- grep(samples$experiment, pattern = public_pattern)

  # Check for missing
  bad_ind <- samples$experiment[c(-fileIndExp, -publicIndExp)]
  if (length(bad_ind)) {
    stop(paste0("Could not find experiment(s): ", paste0(bad_ind, collapse = ",")))
  }

  # Must provide assembly for local files
  assembly <- samples$genome[fileIndExp]
  if (any(is.null(assembly) | is.na(assembly))) {
    stop("If providing local files/directories, specify the genome assembly for those samples.")
  }
  # Custom genome.
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
  if (! "sample_name" %in% colnames(samples)) {
    samples$sample_name <- NA
  }
  samples$read_length <- NA

  ## Handle local files ##
  # TODO: Add bigWig and bedGraph
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
    if (plus_bool) {
      fastq_1 <- gsub(fqNow, pattern = "(.+\\.f[ast]*q)\\+.+", replacement = "\\1")
      if (is.na(samples$sample_name)) {
        samples$sample_name[ind] <- gsub(basename(fastq_1),
                                       pattern = fq_end_pattern, replacement = "")
      }
      samples$paired_end[ind] <- TRUE
      samples$read_length[ind] <- get_fastq_read_length(fastq_1)
    } else {
      samples$sample_name[ind] <- gsub(basename(samples$experiment[ind]),
                                       pattern = fq_end_pattern, replacement = "")
      samples$read_length[ind] <- get_fastq_read_length(samples$experiment[ind])
      samples$paired_end[ind] <- FALSE
    }
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
    samples_public <- samples_public[! is.na(samples_public$experiment),]
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
  if (! "outdir" %in% colnames(samples)) {
    samples$outdir <- outdir
  }

  # Get effective genome size
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
    genome_length <- available_genomes$genome_length[available_genomes$UCSC_orgID ==
                                                               genome_now]
    if (is.na(genome_length)) {
      stop("Genome length is not available for ", genome_now, ". Please manually supply genome length")
    }
    genome_length
  }))

  ## Generate RSeq_variable object for each sample
  vars_list <- lapply(samples$sample_name, FUN = function(exp_now) {
    # # Bug testing
    # exp_now <- samples$sample_name[1]

    vars_now <- new("RSeq_variables")
    samples_now <- samples[samples$sample_name == exp_now, , drop = TRUE]

    # Genome directory
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
    }

    # Others
    vars_now@cores <- as.integer(samples_now$cores)
    vars_now@experiments <- samples_now$experiment
    vars_now@controls <- ifelse(is.na(samples_now$control), "None", samples_now$control)
    vars_now@sample_name <- samples_now$sample_name
    vars_now@out_dir <- samples_now$outdir
    vars_now@paired_end <- samples_now$paired_end
    vars_now@mode <- samples_now$mode
    vars_now@effective_genome_size <- samples_now$effective_genome_size
    vars_now@full_genome_length <- samples_now$full_genome_length
    vars_now@strand_specific <- samples_now$strand_specific
    vars_now@read_length <- samples_now$read_length
    vars_now@moeity <- samples_now$moeity
    vars_now@ip_type <- samples_now$ip_type
    vars_now@no_dedupe <- no_dedupe
    vars_now@no_fastp <- no_fastp
    vars_now@sample_type <- samples_now$file_type
    vars_now@sample_name <- samples_now$sample_name
    slots <- methods::slotNames(vars_now)
    vars_now_list <- purrr::map(slots, ~ methods::slot(object = vars_now, name = .x))
    names(vars_now_list) <- slots
    return(vars_now_list)
  })

  names(vars_list) <- samples$sample_name

  # Generate the JSON output
  dir.create(outdir, showWarnings = FALSE)
  output_json <- file.path(outdir, "rseqVars.json")
  jsonlite::write_json(vars_list, path = output_json)

  return(output_json)

}

# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)



# Source helpers
source(file.path(arg[1], "utils.R"))
# Load required data objects
load(file.path(arg[1], "data", "available_genomes.rda"))
load(file.path(arg[1], "data", "chrom_sizes_list.rda"))

# Get output JSON
result <- processInput(mode = arg[2],
                       outdir = arg[3],
                       genome = arg[4],
                       genome_home_dir = arg[5],
                       cores = arg[6],
                       no_dedupe = as.logical(arg[7]),
                       no_fastp = as.logical(arg[8]),
                       samples = arg[9])
cat(result)



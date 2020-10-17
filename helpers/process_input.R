options(warn=-1) # Prevent warnings

processInput <- function(helpers_dir,
                         mode = NULL,
                         outdir = "RSeq_out/",
                         genome = NULL,
                         genome_home_dir = file.path(path.expand("~"), ".RSeq_genomes"),
                         cores = NULL,
                         no_dedupe = FALSE,
                         no_fastp = FALSE,
                         samples = NULL,
                         keepTmp = FALSE,
                         dryRun = FALSE,
                         dag = FALSE,
                         force = FALSE,
                         reason = FALSE,
                         available_genomes) {

  # ### For bug testing ##
  # mode = NULL
  # outdir = "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_out"
  # genome = NULL
  # genome_home_dir = "/home/UTHSCSA/millerh1/.RSeq_genomes"
  # cores = NULL
  # no_dedupe <- FALSE
  # no_fastp <- FALSE
  # # samples <- "RSeq_CLI/tests/manifest_for_RSeq_testing_11092020.csv"
  # # samples <- "RSeq_CLI/tests/sampleSheet_test5.csv"
  # samples <- data.frame(
  #   "experiment" = "SRR2019281",
  #   "control" = "SRR2019278",
  #   "mode" = "DRIP",
  #   "cores" = 1,
  #   "outdir" = "RSeq_out/",
  #   "genome_home_dir" = "/home/UTHSCSA/millerh1/.RSeq_genomes"
  # )
  # source("~/Bishop.lab/Projects/RSeq/RSeq_CLI/helpers/utils.R")
  # load("~/Bishop.lab/Projects/RSeq/RSeq_CLI/helpers/data/available_genomes.rda")
  # ### #####################

  # Additional data
  # TODO: add support for bigWig and bedGraph
  read_align_pattern <- "\\.bam$|\\.f[ast]*q$"
  public_pattern <- paste0("^[ES]R[RAXSP][0-9]+$|^GS[EM][0-9]+$|^PRJNA[0-9]+$|^SAMN[0-9]+$")
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

  # Unique
  samples <- unique(samples)

  # Remove erroneous rows where control is in the experiment spot
  # TODO: remove for production
  if (! is.null(samples$control)) {samples <- samples[! samples$experiment %in% samples$control, , drop = FALSE]}
  
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
  for (i in seq(samples$experiment)) {
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
    # -- Get the read lengths
    samples$read_length[bamInd] <- unlist(lapply(bamInd, function(ind) {
      if (! file.exists(samples$experiment[ind])) {stop("Could not find bam file ", samples$experiment[ind])}
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
      fastq_1 <- gsub(fqNow, pattern = "(.+\\.f[ast]*q)\\+(.+\\.f[ast]*q)", replacement = "\\1")
      fastq_2 <- gsub(fqNow, pattern = "(.+\\.f[ast]*q)\\+(.+\\.f[ast]*q)", replacement = "\\2")
      if (! file.exists(fastq_1) | ! file.exists(fastq_2)) {stop("Could not find fastq files ", fastq_1, " ", fastq_2)}
      if (is.na(samples$sample_name)) {
        samples$sample_name[ind] <- gsub(basename(fastq_1),
                                       pattern = fq_end_pattern, replacement = "")
      }
      samples$paired_end[ind] <- TRUE
      samples$read_length[ind] <- get_fastq_read_length(fastq_1)
    } else {
      samples$sample_name[ind] <- gsub(basename(samples$experiment[ind]),
                                       pattern = fq_end_pattern, replacement = "")
      if (! file.exists(samples$experiment[ind])) {stop("Could not find fastq file ", samples$experiment[ind])}
      samples$read_length[ind] <- get_fastq_read_length(samples$experiment[ind])
      samples$paired_end[ind] <- FALSE
    }
  }

  ## Process public samples ##
  if (any(samples$file_type == "public")) {

    # Need to standardize with output of public sample search
    samples$experiment_orig <- NA
    samples$control_orig <- NA

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
    orig_ctr_vec <- samples_public$control
    for (i in seq(samples_public$control)) {
      orig_ctr <- samples_public$control[i]
      if (is.na(orig_ctr)) {
        new_ctr_vec <- c(new_ctr_vec, NA)
      } else {
        new_ctr_vec <- c(new_ctr_vec, as.character(
          as.character(sra_info_ctr$experiment[which(sra_info_ctr$accessions_original == orig_ctr)])))
      }
    }
    samples_public$control <- new_ctr_vec
    samples_public$control_orig <- orig_ctr_vec



    # Get info for experimental samples
    samples_public <- samples_public[! is.na(samples_public$experiment),]

    sra_info_exp <- get_public_run_info(samples_public$experiment)

    # -- collate
    public_sample_list <- list()
    for (i in seq(samples_public$experiment)) {
      #sample <- samples_public$experiment[3]
      sample <- samples_public$experiment[i]
      samples_public_now <- samples_public[samples_public$experiment == sample,]
      sra_info_now <- sra_info_exp[sra_info_exp$accessions_original == sample,]
      colnames(sra_info_now)[4] <- "experiment_final"
      sra_info_new <- merge(x = as.data.frame(samples_public_now),
                                y = as.data.frame(sra_info_now),
                            by.x = "experiment", by.y = "accessions_original", all = TRUE)
      # -- Fix EXP info
      sra_info_new <- sra_info_new[, ! colnames(sra_info_new) %in% c("experiment_orig")]
      colnames(sra_info_new)[1] <- "experiment_orig"
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
    # Determine whether Homer annotations available 
    # TODO: What if custom genome?
    available_genomes <- check_homer_anno(available_genomes, genome=vars_now@genome)
    vars_now@homer_anno_available <- available_genomes$homer_anno_available[which(available_genomes$UCSC_orgID == vars_now@genome)]
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
    vars_now@experiments <- as.character(samples_now$experiment)
    vars_now@controls <- as.character(ifelse(is.na(samples_now$control), "None", samples_now$control))
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
    vars_now@helpers_dir <- helpers_dir
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
  samples$outdir <- samples$outdir[1]
  dir.create(samples$outdir[1], showWarnings = FALSE)
  output_json <- file.path(samples$outdir[1], "rseqVars.json")
  output_json <- gsub(output_json, pattern = "//", replacement = "/")
  vars_list[["dryrun"]] <- dryRun
  vars_list[["dag"]] <- dag
  vars_list[["keepTmp"]] <- keepTmp
  vars_list[["force"]] <- force
  vars_list[["reason"]] <- reason
  jsonlite::write_json(vars_list, path = output_json)
  output_csv <- file.path(samples$outdir[1], "rseqVars.csv")
  write.csv(samples, file = output_csv)

  return(output_json)

}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)

# # print(args)
# args <- c("--experiment", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_R1.fastq",
#           "-c", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_ctr_R1.fastq",
#          "-g", "mm10", "-n", "my_experiment", "-o", "RSeq_out/", "-t", "20", "--noDedupe", "-m", "DRIP",
#          "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/helpers")
# args <- c("-e1", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_R1.fastq",
#           "-2", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_R2.fastq",
#           "-c1", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_ctr_R1.fastq",
#           "-c2", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_ctr_R1.fastq",
#           "-m", "qDRIP", "-s", "RSeq_CLI/tests/manifest_for_RSeq_testing_09092020_small.csv",
#           "-g", "hg38", "-n", "my_experiment", "-o", "RSeq_out/", "-t", "20",
#           "--dryRun", "--keepTmp",
#           "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/helpers")

argument_possibles <- data.frame(
  short = c("e", "e1", "e2", "c", "c1", "c2", "m", "g", "n", "s", "o", "G", "t", "r", "v", "h",
            NA, NA, NA, NA, NA, NA, NA),
  long = c("experiment", "experiment_R1", "experiment_R2", "control", "control_R1", "control_R2",
           "mode", "genome", "name", "sampleSheet", "outputDir",
            "genomeDir", "threads", "rseqVars", "version", "help", 
           "dryRun", "noFastp", "force", "reason", "noDedupe", "keepTmp", "dag"),
  type = c(rep("valued", 14), rep("valueless", 9)), stringsAsFactors = FALSE
)

# Loop for parsing command line shell arguments
collect_list <- list()
name <- NULL
unrecognized_arguments <- c()
errors <- c()
for (i in 1:(length(args))) {
  arg <- args[i]
  if (substr(arg, 1, 1) == "-" || grepl(arg, pattern = "RSeq_CLI/helpers")) {
    # CASE: it is a flag -- colect following

    # Collect from collector if name isn't in NULL state (collector will be full)
    if (! is.null(name) && length(type) > 0 && ! type == "valueless") {
      collect_list[[name]] <- collector
    }
    if (grepl(arg, pattern = "RSeq_CLI/helpers")) {break}
    collector <- c()  # re-init collector

    # Get the flag text (could be short-form or long-form flag)
    nameRaw <- gsub(arg, pattern = "^[-]{1,2}", replacement = "")

    # Get the long-form name if not already provided
    name <- argument_possibles$long[which(argument_possibles$short == nameRaw)]
    if (! length(name)) {
      name <- nameRaw
    }

    # Get the argument type
    type <- argument_possibles$type[which(argument_possibles$short == name)]
    if (! length(type)) {
      type <- argument_possibles$type[which(argument_possibles$long == name)]
      if (! length(type)) {
        unrecognized_arguments <- c(unrecognized_arguments, arg)
        next
      }
    }
    # Collect valueless arguments as logical
    if (type == "valueless") {
      collect_list[[name]] <- TRUE
    }
    next
  } else {
    collector <- c(collector, arg)
  }
}

if (length(unrecognized_arguments)) {
  errors <- c(errors, paste0("Unrecognized argument(s): '",paste0(unrecognized_arguments, collapse = "', '"), "'."))
}


# Exit for usage and version info
if (! is.null(collect_list$help)) {cat("usage"); quit(status = 0, save = "no")}
if (! is.null(collect_list$version)) {cat("version"); quit(status = 0, save = "no")}

# Check if rseqVars given -- bypass processing if so
rseqVars <- collect_list$rseqVars
if (is.null(rseqVars)) {
  # Compile experiment and control 
  experiment <- collect_list$experiment
  experiment_R1 <- collect_list$experiment_R1
  experiment_R2 <- collect_list$experiment_R2
  control <- collect_list$control
  control_R1 <- collect_list$control_R1
  control_R2 <- collect_list$control_R2
  
  # Set experiment values
  if (! is.null(experiment_R1)) {
    if (! is.null(experiment_R2)) {
      experiment <- paste0(experiment_R1, "+", experiment_R2)
    } else {
      errors <- c(errors, "Fastq Experiment R1 provided, but not R2!")
    }
  }
  
  # Set control values
  if (! is.null(control_R1)) {
    if (! is.null(control_R2)) {
      control <- paste0(control_R1, "+", control_R2)
    } else {
      errors <- c(errors, "Fastq control R1 provided, but not R2!")
    }
  }
  
  
  # Other variables
  mode <- collect_list$mode
  outdir <- collect_list$outputDir
  genome <- collect_list$genome
  genome_home_dir <- collect_list$genome_home_dir
  samples <- collect_list$sampleSheet
  cores <- collect_list$threads
  helpers_dir <- args[length(args)]
  
  # Read in samples if they exist
  if (! is.null(samples)) {samples <- read.csv(samples)}
  
  # Collect errors
  if (is.null(experiment) && is.null(samples$experiment)) {errors <- c(errors, "No experiment or sampleSheet specified!"); experiment <- NA}
  if (is.null(mode) && is.null(samples$mode)) {errors <- c(errors, "No mode specified!"); mode <- NA}
  
  # Set defaults
  if (is.null(outdir)) {outdir <- "RSeq_out/"}
  if (is.null(genome_home_dir)) {genome_home_dir <- file.path(path.expand("~"), ".RSeq_genomes")}
  if (is.null(cores)) {cores <- 1} else {cores <- as.numeric(cores)}
  
  # Compile sampleSheet
  if (is.null(samples)) {samples <- data.frame(experiment)}
  if (is.null(samples$mode)) {samples$mode <- mode}
  if (is.null(samples$outdir)) {samples$outdir <- outdir}
  if (is.null(samples$cores)) {samples$cores <- cores}
  if (is.null(samples$genome_home_dir)) {samples$genome_home_dir <- genome_home_dir}
  
  
  # Check for control and genome
  if (is.null(samples$control) && ! is.null(control)) {
    if (length(control) != length(samples$experiment)) {
      errors <- c(errors, "Number of experiment samples is different from the number of controls!")
    } else {
      samples$control <- control
    }
  }
  if (! is.null(genome)) {samples$genome <- genome}
  
  # Handle any errors and return to shell
  if (length(errors) > 0) {
    cat(paste0("ERROR(s):\n - ", paste0(errors, collapse = " \n - ")))
    quit(status = 1, save = "no")
  }
  
  
  # Source helpers
  source(file.path(helpers_dir, "utils.R"))
  # Load required data objects
  load(file.path(helpers_dir, "data", "available_genomes.rda"))
  
  if (! "keepTmp" %in% names(collect_list)) {keepTmp <- FALSE} else {keepTmp <- collect_list$keepTmp}
  if (! "dryRun" %in% names(collect_list)) {dryRun <- FALSE} else {dryRun <- collect_list$dryRun}
  if (! "dag" %in% names(collect_list)) {dag <- FALSE} else {dag <- collect_list$dag}
  if (! "force" %in% names(collect_list)) {force <- FALSE} else {force <- collect_list$force}
  if (! "reason" %in% names(collect_list)) {reason <- FALSE} else {reason <- collect_list$reason}
  
  # # Get output JSON
  result <- processInput(samples = samples, keepTmp = keepTmp, dryRun = dryRun,
                         helpers_dir = helpers_dir, force = force, dag = dag, reason=reason,
                         available_genomes = available_genomes)
} else {
  if (! file.exists(rseqVars)) {
    stop("Cannot locate the user-supplied rseqVars file: ", rseqVars)
  }
  rseqVarList <- jsonlite::read_json(rseqVars, simplifyVector = TRUE)
  if (! "keepTmp" %in% names(collect_list)) {keepTmp <- FALSE} else {keepTmp <- collect_list$keepTmp}
  if (! "dryRun" %in% names(collect_list)) {dryRun <- FALSE} else {dryRun <- collect_list$dryRun}
  if (! "dag" %in% names(collect_list)) {dag <- FALSE} else {dag <- collect_list$dag}
  if (! "force" %in% names(collect_list)) {force <- FALSE} else {force <- collect_list$force}
  if (! "reason" %in% names(collect_list)) {reason <- FALSE} else {reason <- collect_list$reason}
  if (! "threads" %in% names(collect_list)) {cores <- 1} else {cores <- as.numeric(collect_list$threads)}
  rseqVarList[["dryrun"]] <- dryRun
  rseqVarList[["dag"]] <- dag
  rseqVarList[["keepTmp"]] <- keepTmp
  rseqVarList[["force"]] <- force
  rseqVarList[["reason"]] <- reason
  # Change the number of cores to match current arguments
  for (xnow in names(rseqVarList)) {
    if (! xnow %in% c("dryrun", "dag", "keepTmp", "force", "reason")) {
      rseqVarList[[xnow]]$cores <- cores
    }
  }
  rseqVarList[["cores"]] <- cores
  od <- rseqVarList[[1]]$out_dir
  rseqVarList[["dryrun"]]
  result <- file.path(od, "rseqVars.json")
  jsonlite::write_json(rseqVarList, result)
}


cat(result)



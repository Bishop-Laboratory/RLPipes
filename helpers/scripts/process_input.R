#globals
#options(warn=-1) # Prevent warnings

#libraries
library(tidyverse)

processInput <- function(samples,
                         available_genomes) {

  # ### For bug testing ##
  # samples <- read.csv("tests/manifest_for_RSeq_testing_09092020.csv", fileEncoding="UTF-8-BOM")
  # source("helpers/utils.R")
  # load("helpers/data/available_genomes.rda")
  # ### #####################

  suppressPackageStartupMessages(require(dplyr))

  # Additional data
  # TODO: add support for bigWig and bedGraph
  read_align_pattern <- "\\.bam$|\\.f[ast]*q$"
  public_pattern <- paste0("^[ES]R[RAXSP][0-9]+$|^GS[EM][0-9]+$|^PRJNA[0-9]+$|^SAMN[0-9]+$")
  R1_pattern <- "[\\._]+[rR]*1\\.f[ast]*q$"
  R2_pattern <- "[\\._]+[rR]*2\\.f[ast]*q$"
  fq_end_pattern <- paste0(R1_pattern, "|", R2_pattern, "|[\\._]+[rR]*\\*.f[ast]*q$|\\.f[ast]*q$")

  mode_df <- data.frame(
    mode = c("DRIP", "ssDRIP", "qDRIP",
             "sDRIP", "S1-DRIP", # "bisDRIP", "SMRF",
             "RDIP", "DRIPc",
             "R-ChIP", "RR-ChIP", "DRIVE",
             "MapR", "RNH-CnR"),
    strand_specific = c(FALSE, rep(TRUE, 3), FALSE, rep(TRUE, 4), rep(FALSE, 3)),
    ip_type = c(rep("S9.6", 7), rep("RNaseH1", 5)),
    moeity = c(rep("DNA", 5), rep("RNA", 7)), stringsAsFactors = FALSE
  )
  modes <- c(mode_df$mode, "custom")

  # Unique
  samples <- unique(samples)

  # Collate controls as separate experiments
  if (! "control" %in% colnames(samples)) {samples$control <- NA}
  samples_ctr <- samples %>%
    filter(! is.na(control)) %>%
    mutate(experiment = control) %>%
    mutate(control = NA)
  if (length(samples_ctr$experiment)) {
    samples <- samples %>%
      dplyr::bind_rows(samples_ctr)
  }


  # Check for arguments. If not set, replace with defaults
  if (! "paired_end" %in% colnames(samples)) {samples$paired_end <- NA}
  if (! "effective_genome_size" %in% colnames(samples)) {samples$effective_genome_size <- NA}
  if (! "full_genome_length" %in% colnames(samples)) {samples$full_genome_length <- NA}
  if (! "genome" %in% colnames(samples)) {samples$genome <- NA}


  # Deal with strand specificity and mode type
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
  # Catch failure to provide genome for local files.
  badGenomeSpec <- which(is.na(samples$genome) & samples$file_type != "public")
  if (length(badGenomeSpec)) {
    stop("User must supply a UCSC genome id or genome directory location under 'genome'.",
         " For sample(s) ", samples$experiment[badGenomeSpec], ", the genome supplied is ",
         samples$genome[badGenomeSpec])
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
          as.character(sra_info_ctr$sample_name[which(sra_info_ctr$accessions_original == orig_ctr)])))
      }
    }
    samples_public$control <- new_ctr_vec
    samples_public$control_orig <- orig_ctr_vec

    # Get info for experimental samples
    samples_public <- samples_public[! is.na(samples_public$experiment),]

    # TODO: More efficient way without rerunning all the control samples
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
    samples_public <- dplyr::bind_rows(public_sample_list)
    samples_private <- samples[samples$file_type != "public",]
    if(! all(colnames(samples_public) == colnames(samples_private))) {
      stop("Bug at public sample collation -- please notify package author. Should be unreachable.")
    }
    samples <- rbind(samples_public, samples_private)
  }

  # dups <- duplicated(samples$sample_name) & samples$file_type != "public"
  # if (any(dups)) {
  #   stop("Sample names must be unique. Please remove duplicates: ",
  #        paste0(samples$sample_name[dups], collapse = ", "))
  # }

  # Mark control samples
  samples <- samples %>%
    mutate(ip_type = case_when(
      sample_name %in% (samples %>% filter(! is.na(control)) %>%
                          pull(control)) ~ "Input",
      TRUE ~ ip_type
    ))

  # Get effective genome size
  rownames(available_genomes) <- available_genomes$UCSC_orgID
  available_genome_sizes <- available_genomes[,grep(colnames(available_genomes),
                                                    pattern = "eff_genome_size")]
  sizes <- as.numeric(gsub(colnames(available_genome_sizes),
                         pattern = "eff_genome_size_([0-9]+)bp",
                         replacement = "\\1"))
  samples$effective_genome_size[is.na(samples$effective_genome_size)] <- unlist(
      lapply(samples$sample_name[is.na(samples$effective_genome_size)],
             function(sample_now) {
      len_now <- samples$read_length[samples$sample_name == sample_now]
      colInd <- which.min(abs(sizes-len_now))
      genome_now <- samples$genome[samples$sample_name == sample_now]
      eff_genome_size <- available_genome_sizes[genome_now, colInd]
      if (is.na(eff_genome_size)) {
        stop("Genome size is not available for ", genome_now,
             ". Please manually supply effective genome size as a",
             "column in 'samples' or variable in 'makeRSeqDataSet'.",
             " Genome sizes may be calculated with the 'unique-kmers.py",
             " script from the khmer package'")
      }
      eff_genome_size
    })
  )

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

  # Finish compiling data
  configs <- samples %>%
    mutate(homer_anno_available = genome %in% (available_genomes %>%
                                                 filter(homer_anno_available) %>%
                                                 pull(UCSC_orgID))) %>%
    mutate(control = as.character(ifelse(is.na(control), "None", control))) %>%
    as.list()

  # TODO: Handle custom genomes
  # if (dir.exists(vars_now@genome)) {
  #   vars_now@genome_home_dir <- vars_now@genome
  #   if (is.na(samples_now$effective_genome_size)) {
  #     # TODO: Run unique-kmers.py and save as a txt file in that dir if missing.
  #     stop("No effective genome size provided for user-supplied genome: ", vars_now@genome, "../..",
  #          " Please calculate effective genome sizes using unique-kmers.py from the khmer package",
  #          " and supply as a column in 'samples' or as a vector to 'makeRSeqDataSet'.")
  #   }
  # } else {
  #   vars_now@genome_home_dir <- samples_now$genome_home_dir
  # }

  # Generate the JSON output
  return(configs)
}

# Parse shell args
args <- commandArgs(trailingOnly=TRUE)
# print(paste0(args, collapse = "', '"))
# print(args)
# args <- c("-r", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/configs.json",
#           '-o', "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/",
#           '-t', '80', '--dag', "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/helpers")

# args <- c("-r", "RSeq_out/config.json",
#           "-o", "outdir2",
#           "-t", "10")

# args <- c("--experiment", "tests/qDRIP_R1.fastq",
#           "-c", "tests/qDRIP_ctr_R1.fastq",
#          "-g", "mm10", "-n", "my_experiment", "-o", "RSeq_out/", "-t", "20", "-m", "DRIP",
#          "helpers/")
# args <- c("-e1", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_R1.fastq",
#           "-2", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_R2.fastq",
#           "-c1", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_ctr_R1.fastq",
#           "-c2", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/tests/qDRIP_ctr_R1.fastq",
#           "-m", "qDRIP", "-s", "RSeq_CLI/tests/manifest_for_RSeq_testing_09092020_small.csv",
#           "-g", "hg38", "-n", "my_experiment", "-o", "RSeq_out/", "-t", "20",
#           "--dryrun", "--keepTmp",
#           "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/helpers")
# args <- c("-e", "whatevs" , "-c",
#            "--dag","-t","/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/helpers")
# args <- c('-e', 'GSM4276887', 'GSM4276888', '-c', 'GSM4276889', 'GSM4276890', '-o', '../RMapDB/data/', '-g', 'hg38',
#           '-t', '80', "-d", "10000", '--dryrun', '-m', 'qDRIP', '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/helpers')

# args <- c('-s', '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/rmapsamples_10_05_2020.csv',
#           '--dryrun', '-o', '../RMapDB/data/', '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/helpers')

# args <- c("-r", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/tests/RSeq_out3/configs.json",
#           "-o", "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/tests/RSeq_out9/",
#           '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/helpers')

# args <- c("-s", "tests/sampleSheet_test9.csv", "-m", "DRIP",  "-S", "dryrun=False", "targets=[asd,", "asdd]",
#           "-t", "30", "helpers/")

# args <- c('-e', 'GSM4276887', 'GSM4276888', '-c', 'GSM4276889', 'GSM4276890', '-o', 'data/', '-g', 'hg38',
#           '-t', '20', '--dryrun', '-m', 'qDRIP', '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/helpers')

# args <- c("-s", "tests/manifest_for_RSeq_testing_09092020.csv",
#           "-m", "DRIP",
#           "-t", "30",
#           "-S", "dryrun=True", "notemp=True", "printreason=True",
#           "-o", "test4", "helpers")

# args <- c("-s", "tests/manifest_for_RSeq_testing_09092020_small.csv",
#           "-o", "test5", "helpers")

# args <- c("-s", "tests/manifest_for_RSeq_testing_11092020.csv",
#           "-o", "test7", "helpers")

# args <- c("-e", "SRX1025890", "-m", "DRIP", "-o", "outfodler", "-b", "/mnt/c/Users/mille/RSeq/helpers")

# Dataframe of mappings between possible arguments and whether they have a value or not
argument_possibles <- data.frame(
  short = c("e", "e1", "e2", "c", "c1", "c2", "d", "m", "g", "n", "s", "o", "G", "t", "r", "S", "b", "v", "h"),
  long = c("experiment", "experiment_R1", "experiment_R2", "control", "control_R1", "control_R2",
           "downsample",  "mode", "genome", "name", "sampleSheet", "outdir",
           "genome_home_dir", "threads", "configs", "snake_args", "bwa_mem2", "version", "help"),
  type = c(rep("valued", 16), rep("valueless", 3)), stringsAsFactors = FALSE
)

# Loop for parsing command line shell arguments
collect_list <- list()
name <- NULL
unrecognized_arguments <- c()
errors <- c()
collector <- c()  # init collector
for (i in 1:(length(args))) {
  arg <- args[i]
  if ((substr(arg, 1, 1) == "-" && ! grepl(substr(arg, 2, 2), pattern = "[0-9]+|[Ii]+")) ||
      grepl(arg, pattern = "RSeq/helpers")) {
    # CASE: it is a flag -- collect following

    # Collect from collector if name isn't in NULL state (collector will be full)
    if (! is.null(name) && length(type) > 0 && type != "valueless") {
      if (is.null(collector)) {
        # if collector is null, this means a valued flag was given but no value supplied
        errors <- c(errors,  paste0("'", args[i-1], "' was invoked but no value was supplied!"))
        break
      }
      collect_list[[name]] <- collector
    }
    if (grepl(arg, pattern = "RSeq/helpers")) {break}
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

# TODO: Need better solution to this which can handle arbitrarily-nested dicts/lists
# Check for list/dict
if (! is.null(collect_list$snake_args)) {
  snkarg <- collect_list$snake_args
  snkarg_final <- c()
  collect <- FALSE
  item_list <- NA
  for (i in seq(snkarg)) {
    arg <- snkarg[i]
    if (grepl(arg, pattern = "\\=\\[|\\=\\{")) {
      collect <- TRUE
      item_list <- arg
    } else if (collect) {
      item_list <- c(item_list, arg)
    } else {
      snkarg_final <- c(snkarg_final, arg)
    }

    if (grepl(arg, pattern = "\\}$|\\]$")) {
      compiled_item <- paste0(item_list, collapse = " ")
      compiled_item <- gsub(compiled_item, pattern = "\\[", replacement = "['")
      compiled_item <- gsub(compiled_item, pattern = "\\]", replacement = "']")
      compiled_item <- gsub(compiled_item, pattern = "\\{", replacement = "{'")
      compiled_item <- gsub(compiled_item, pattern = "\\}", replacement = "'}")
      compiled_item <- gsub(compiled_item, pattern = ", ", replacement = "', '")
      snkarg_final <- c(snkarg_final, compiled_item)
      collect <- FALSE
    }
  }
  collect_list$snake_args <- snkarg_final
}

if (length(unrecognized_arguments)) {
  errors <- c(errors, paste0("Unrecognized argument(s): '",paste0(unrecognized_arguments, collapse = "', '"), "'."))
}

# Exit for usage and version info
if (! is.null(collect_list$help)) {cat("usage"); quit(status = 0, save = "no")}
if (! is.null(collect_list$version)) {cat("version"); quit(status = 0, save = "no")}

# Set global defaults
outdir <- collect_list$outdir
threads <- collect_list$threads
genome_home_dir <- collect_list$genome_home_dir
snake_args <- collect_list$snake_args
bwa_mem2 <- collect_list$bwa_mem2
if (is.null(snake_args)) {snake_args <- "use_conda=True"}
helpers_dir <- args[length(args)]
# Set defaults
if (is.null(outdir)) {outdir <- "rseq_out/"}
if (is.null(genome_home_dir)) {genome_home_dir <- file.path(path.expand("~"), ".rseq_genomes")}
if (is.null(threads)) {threads <- 1} else {threads <- as.numeric(threads)}
if (is.null(snake_args)) {snake_args <- ""}
if (is.null(bwa_mem2)) {bwa_mem2 <- FALSE}
#create genome home directory if it does not exist
dir.create(genome_home_dir, showWarnings = FALSE)
#create output directory and specify path to json output
dir.create(outdir, showWarnings = FALSE)
output_json <- file.path(outdir, "config.json")
output_json <- gsub(output_json, pattern = "//", replacement = "/")

# Check if configs given -- bypass processing if so
configs <- collect_list$configs
if (is.null(configs)) {
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
  genome <- collect_list$genome
  samples <- collect_list$sampleSheet

  # Read in samples if they exist
  if (! is.null(samples)) {
    # UTF-8-BOM option removes the byte order mark from Windows-generated files
    # See: https://stackoverflow.com/questions/21624796/read-a-utf-8-text-file-with-bom
    samples <- read.csv(samples, fileEncoding = "UTF-8-BOM")
    if ("X" %in% colnames(samples)) {
      samples <- samples[, ! colnames(samples) == "X", drop=FALSE]
    }
  }


  # Collect errors
  if (is.null(experiment) && is.null(samples$experiment)) {errors <- c(errors, "No experiment or sampleSheet specified!"); experiment <- NA}
  if (is.null(mode) && is.null(samples$mode)) {errors <- c(errors, "No mode specified!"); mode <- NA}

  # Compile sampleSheet
  if (is.null(samples)) {samples <- data.frame(experiment)}
  if (is.null(samples$mode)) {samples$mode <- mode}

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
  source(file.path(helpers_dir, "scripts", "utils.R"))
  # Load required data objects
  load(file.path(helpers_dir, "data", "available_genomes.rda"))


  # # # Debug
  # readr::write_csv(samples, "samples.csv")
  # samples <- readr::read_csv("samples.csv")

  # # Get output JSON
  configs <- processInput(samples = samples,
                          available_genomes = available_genomes)
  configs <- c(configs, list(
    outdir=outdir,
    genome_home_dir=genome_home_dir,
    threads=threads,
    helpers_dir=helpers_dir,
    snake_args=snake_args,
    bwa_mem2=bwa_mem2
  ))

  jsonlite::write_json(configs, path = output_json)
} else {
  if (! file.exists(configs)) {
    errors <- c(errors, paste0("Cannot locate the user-supplied configs file: ", configs))
  }
  if (length(errors) > 0) {
    cat(paste0("ERROR(s):\n - ", paste0(errors, collapse = " \n - ")))
    quit(status = 1, save = "no")
  }
  configs <- jsonlite::read_json(configs, simplifyVector = TRUE)
  configs[['outdir']] <- outdir
  configs[['genome_home_dir']] <- genome_home_dir
  configs[['helpers_dir']] <- helpers_dir
  configs[['threads']] <- threads
  configs[['snake_args']] <- snake_args
  configs[['bwa_mem2']] <- bwa_mem2
  jsonlite::write_json(configs, path = output_json)
}


cat(output_json)

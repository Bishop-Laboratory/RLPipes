# Compile peaks from MACS2 and EPIC2 based on mode of sequencing
compile_peaks <- function(mode, control, macs2, epic2, output_prefix) {

  # mode <- "DRIP"
  # control <- "None"
  # macs2 <- "rseq_out_input/peaks/SRX1025893_TC32_Input/macs2/SRX1025893_TC32_Input_hg38__peaks.broadPeak"
  # epic2 <- "rseq_out_input/peaks/SRX1025893_TC32_Input/epic2/SRX1025893_TC32_Input_hg38__peaks_epic2.bed"
  # output_prefix <- "rseq_out_input/peaks/SRX1025890_TC32_NT_DRIP/SRX1025890_TC32_NT_DRIP_hg38__compiled_peaks"

  # TODO: Need to revisit this approach. Simpler would be to just ouput all macs2 that ol with epic2
  # TODO: We could also test to see the # of peaks called by each approach on input samples and estimate error rate
  # TODO: Should output pval as same format for macs2 and epic2
  control <- control != "None"
  no_macs <- FALSE
  no_epic <- FALSE
  
  # Initialize empty GR
  empty_gr <- as.data.frame(t(matrix(c("chrM", "0", "10", "Empty", "..", "*"))))
  colnames(empty_gr) <- c("seqnames", "start", "end", "name", "score", "strand")
  empty_gr <- ChIPpeakAnno::toGRanges(empty_gr)
  
  macs2 <- readr::read_tsv(macs2, col_names = c("seqnames", "start", "end", "name", "score", "strand",
                                                "signalValue", "pValue", "qValue"))
  if (! length(colnames(macs2))) {
    no_macs <- TRUE
    macs2 <- empty_gr
    warning("No MACS2 peaks found!")
  } else {
    macs2 <- dplyr::filter(macs2, start < end & ! is.na(score))
    macs2 <- ChIPpeakAnno::toGRanges(as.data.frame(macs2))
    macs2$score <- as.numeric(macs2$score)
  }
  epic2 <- readr::read_tsv(epic2)

  if (control) {
    colnames(epic2) <- c("seqnames", "start", "end", "pValue", "score", "strand",
                                                "ChIPCount", "InputCount", "FDR", "log2FoldChange")
  } else {
    colnames(epic2) <- c("seqnames", "start", "end", "ChIPCount", "score", "strand")
  }

  if (! length(colnames(epic2))) {
    no_epic <- TRUE
    epic2 <- empty_gr
    warning("No EPIC2 peaks found!")
  } else {
    epic2$name <- seq(length(epic2$start))
    epic2 <- dplyr::filter(epic2, start < end & ! is.na(score))
    epic2 <- ChIPpeakAnno::toGRanges(as.data.frame(epic2))
    epic2$score <- as.numeric(epic2$score)
  }
  
  # Overlap and save in RData file
  peak_ol <- ChIPpeakAnno::findOverlapsOfPeaks(macs2, epic2, ignore.strand = FALSE)
  save(peak_ol, file = paste0(output_prefix, ".rda"))
  outbed <- peak_ol$overlappingPeaks$`macs2///epic2`

  if (no_macs & no_epic) {
    warning("No peaks found! Outputting empty range set.")
    output <- "epic"
  } else if (no_macs) {
    warning("Outputting epic2 ranges only.")
    output <- "epic"
  } else if (no_epic) {
    warning("Outputting macs2 ranges only.")
    output <- "macs"
  } else {
    if (mode %in% c("DRIP", "DRIPc", "sDRIP", 'qDRIP', 'RDIP', 'ssDRIP', 'S1-DRIP')) {
      # CASE: DRIP family / S9.6 based
      if (control) {
        # CASE: Control sample is available -- using the overlap
        if (! length(outbed$peaks1)) {
          # CASE No overlapping peaks
          output <- "macs"
        } else {
          # CASE: Both peaksets available, taking macs2 that overlaps with epic2 and saving
          output <- "filtered_macs"
        }
      } else {
        # CASE: No contol, using EPIC2 only
        output <- "epic"
      }
    } else if (mode %in% c("R-ChIP", "RR-ChIP", "RNH-CnR", "MapR", 'DRIVE')) {
      # CASE: RNH1 family of modes
      if (! length(outbed$peaks1)) {
        output <- "macs"
      } else {
        # CASE: epic2 and macs2 available, select macs2 peaks that overlap with epic2
        output <- "filtered_macs"
      }
    }
  }

  if (output == "filtered_macs") {
    outbedgr <- as.data.frame(macs2[names(macs2) %in% outbed$peaks1,])
    colnames(outbedgr)[1] <- paste0("#", colnames(outbedgr)[1])
    outbedgr$name <- paste0("peak_", seq(outbedgr$width))
    outbedgr <- outbedgr[,c("#seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")]
    readr::write_tsv(outbedgr, file = paste0(output_prefix, ".bed"))
  } else if (output == "macs") {
    macs2 <- as.data.frame(macs2)
    colnames(macs2)[1] <- paste0("#", colnames(macs2)[1])
    macs2$name <- paste0("peak_", seq(macs2$width))
    macs2 <- macs2[,c("#seqnames", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue")]
    readr::write_tsv(macs2, file = paste0(output_prefix, ".bed"))
  } else if (output == "epic") {
    epic2 <- as.data.frame(epic2)
    colnames(epic2)[1] <- paste0("#", colnames(epic2)[1])
    epic2$name <- paste0("peak_", seq(epic2$width))
    epic2 <- epic2[,c("#seqnames", "start", "end", "name", "score", "strand", "ChIPCount")]
    readr::write_tsv(epic2, file = paste0(output_prefix, ".bed"))
  }
}

# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)
# print(arg)

suppressMessages(compile_peaks(mode = arg[1],
                               control = arg[2],
                               macs2 = arg[3],
                               epic2 = arg[4],
                               output_prefix = arg[5]))

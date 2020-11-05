# Compile peaks from MACS2 and EPIC2 based on mode of sequencing
compile_peaks <- function(configs, sample_name, macs2, epic2, output_rda, output_bed) {
  # configs <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/rseqVars.json"
  # sample_name <- 'SRX2675009_HKE293-WKKD-V5ChIP'
  # macs2 <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/SRX113813_Ntera2_DNA/peaks_macs_unstranded/SRX113813_Ntera2_DNA_hg38.unstranded.broadPeak"
  # epic2 <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/SRX113813_Ntera2_DNA/peaks_epic_unstranded/SRX113813_Ntera2_DNA_hg38.unstranded.bed"
  # macs2 <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/SRX6427717_DMSO_qDRIP-seq_1/peaks_macs_stranded/SRX6427717_DMSO_qDRIP-seq_1_hg38.stranded.broadPeak"
  # epic2 <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/SRX6427717_DMSO_qDRIP-seq_1/peaks_epic_stranded/SRX6427717_DMSO_qDRIP-seq_1_hg38.stranded.bed"  # output <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/tests/RSeq_out3/qDRIP/peaks_final_unstranded/qDRIP_hg38.bed"
  # macs2 <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/SRX2675009_HKE293-WKKD-V5ChIP/peaks_macs_unstranded/SRX2675009_HKE293-WKKD-V5ChIP_hg38.unstranded.broadPeak"
  # epic2 <- "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/SRX2675009_HKE293-WKKD-V5ChIP/peaks_epic_unstranded/SRX2675009_HKE293-WKKD-V5ChIP_hg38.unstranded.bed"
  configs <- gsub(configs, pattern = "//", replacement = "/")
  configlist <- jsonlite::read_json(configs, simplifyVector = TRUE)
  configlist <- configlist[[sample_name]]
  
  mode <- configlist$mode
  control <- ifelse(configlist$controls == 'None', FALSE, TRUE)
  
  macs2 <- readr::read_tsv(macs2, col_names = c("seqnames", "start", "end", "name", "score", "strand"))
  macs2 <- dplyr::filter(macs2, start < end & ! is.na(score))
  MACS2 <- ChIPpeakAnno::toGRanges(as.data.frame(macs2))
  MACS2$score <- as.numeric(MACS2$score)
  
  epic2 <- readr::read_tsv(epic2, col_names = c("seqnames", "start", "end", "name", "score", "strand"))
  epic2$name <- seq(length(epic2$name))
  epic2 <- dplyr::filter(epic2, start < end & ! is.na(score))
  EPIC2 <- ChIPpeakAnno::toGRanges(as.data.frame(epic2))
  EPIC2$score <- as.numeric(EPIC2$score)
  
  peak_ol <- ChIPpeakAnno::findOverlapsOfPeaks(MACS2, EPIC2, ignore.strand = FALSE)
  save(peak_ol, file = output_rda)
  
  if (mode %in% c("DRIP", "DRIPc", "sDRIP", 'qDRIP', 'RDIP', 'ssDRIP')) {
    if (control) {
      outbed <- peak_ol$overlappingPeaks$`MACS2///EPIC2`
      if (! length(outbed$peaks1)) {
        rtracklayer::export(MACS2, con = output_bed)
      } else {
        outbedgr <- ChIPpeakAnno::toGRanges(outbed[,c(2:4, 7, 6)])
        rtracklayer::export(outbedgr, con = output_bed)
      }
    } else {
      rtracklayer::export(EPIC2, con = output_bed)
    }
  } else if (mode %in% c("R-ChIP", "RR-ChIP", "RNH-CnR", "MapR", 'DRIVE')) {
    outbed <- peak_ol$overlappingPeaks$`MACS2///EPIC2`
    if (! length(outbed$peaks1)) {
      rtracklayer::export(MACS2, con = output_bed)
    } else {
      outbedgr <- ChIPpeakAnno::toGRanges(outbed[,c(2:4, 7, 6)])
      rtracklayer::export(outbedgr, con = output_bed)
    }
  }
  
  # bed_awk <- "awk ' NR>1 {print $1\"\t\"$2\"\t\"$3\"\t\"\".\"\"\t\"\".\"\"\t\"$6} '"  # Gives the file in bed6 format
  # if (mode %in% c("DRIP", "DRIPc", "sDRIP", 'qDRIP', 'RDIP', 'ssDRIP')) {
  #   if (control) {
  #     cmd <- paste0("bedtools intersect -a ", macs2, " -b ", epic2, " | ", bed_awk, " > ", output)
  #   } else {
  #     cmd <- paste0(bed_awk, " ", epic2, " > ", output)
  #   }
  # } else if (mode %in% c("R-ChIP", "RR-ChIP", "RNH-CnR", "MapR", 'DRIVE')) {
  #   cmd <- paste0("bedtools intersect -a ", macs2, " -b ", epic2, " | ", bed_awk, " > ", output)
  # }
  # system(cmd)
}


# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)
# print(arg)

suppressMessages(compile_peaks(configs = arg[1],
                               sample_name = arg[2],
                               macs2 = arg[3],
                               epic2 = arg[4],
                               output_rda = arg[5],
                               output_bed = arg[6]))

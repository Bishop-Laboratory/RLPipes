#' Calculate enrichment of peaks at RLFS
perm_test <- function(cores, genome, peaks, rlfs, outfile) {
  
  # cores = 80
  # genome = "mm10"
  # peaks = "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RMapDB/data/SRX1674681_3T3_DRIPc-seq/peaks_final_stranded/SRX1674681_3T3_DRIPc-seq_mm10.stranded.bed"
  # rlfs = "~/.RSeq_genomes/mm10/rloop_predictions/mm10.rlfs.bed"
  
  peaks <- regioneR::toGRanges(peaks)
  rlfs <- regioneR::toGRanges(rlfs)
  # Prevent stranded assignment
  GenomicRanges::strand(rlfs) <- "*"
  rlfs <- GenomicRanges::reduce(rlfs)
  # Mask makes no difference. Turns out R-ChIP is finding genuine R-loops.
  chrom_sizes <- readr::read_tsv(paste0('http://hgdownload.soe.ucsc.edu/goldenPath/',
                                        genome, '/bigZips/', genome, '.chrom.sizes'),
                                 col_names = FALSE)
  pt <- regioneR::permTest(A=peaks, B=rlfs, ntimes=500, force.parallel = TRUE, 
                           genome=as.data.frame(chrom_sizes), mc.cores=cores, allow.overlaps = FALSE,
                           randomize.function=regioneR::circularRandomizeRegions, 
                           evaluate.function=regioneR::numOverlaps, alternative = "greater")
  lz <- regioneR::localZScore(pt=pt, A=peaks, B=rlfs, window = 5000, step = 50)
  save(pt, lz, file = outfile)
}

# Parse shell args
arg <- commandArgs(trailingOnly=TRUE)
# print(arg)

# Get output JSON
suppressMessages(perm_test(cores = arg[1],
                           genome = arg[2],
                           peaks = arg[3],
                           rlfs = arg[4],
                           outfile = arg[5]))

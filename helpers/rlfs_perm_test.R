#' Calculate enrichment of peaks at RLFS
perm_test <- function(cores, genome, peaks, rlfs, outfile) {
  peaks <- regioneR::toGRanges(peaks)
  rlfs <- regioneR::toGRanges(rlfs)
  # Mask makes no difference. Turns out R-ChIP is finding genuine R-loops.
  genome <- regioneR::getGenome(genome)
  pt <- regioneR::permTest(A=peaks, B=rlfs, ntimes=500, force.parallel = TRUE, 
                           genome=genome, mc.cores=cores, allow.overlaps = FALSE,
                           randomize.function=regioneR::circularRandomizeRegions, 
                           evaluate.function=regioneR::numOverlaps, alternative = "greater")
  lz <- regioneR::localZScore(pt=pt, A=peaks, B=rlfs, window = 5000)
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

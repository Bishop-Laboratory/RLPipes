##### Correlation analysis #####
library(tidyverse)
library(ChIPpeakAnno)

# Get the metadata
sample_metadata <- read_csv("misc/rmap_full_11_25.csv")

# Get the human genome
download.file('http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes',
              destfile = "misc/conservation_analysis/data/hg38.chrom.sizes")

# Divide the genome into 10kb windows with a 2.5kb step size 
if (! file.exists("misc/conservation_analysis/data/hg38.10kb.tiles.bed")) {
  cmd <- paste0('bedtools makewindows -g misc/conservation_analysis/data/hg38.chrom.sizes ',
                '-w 10000 -s 2500 > misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.bed')
  system(cmd)
}

# Count number of overlaps with human experimental samples (removes control & non-human samples)
samples_to_intersect <- sample_metadata$sample_name[
  which(sample_metadata$genome == "hg38" &
          ! sample_metadata$Condition %in% c("WKKD", "IgG", "RNH", "ACTD"))
]

# Compile file paths for these peaks
peaks_to_intersect <- paste0('misc/conservation_analysis/data/final-peaks-unstranded/',
                             samples_to_intersect, '_hg38.unstranded.bed', collapse = " ")

# Count the overlaps between peaks and genome widows
cmd <- paste0('bedtools intersect -a misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.bed -b ',
              peaks_to_intersect,
              " -c > misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed")
if (! file.exists('misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed')) {
  system(cmd)
}

# Read in the results as a GRanges object
rl_cons <- read_tsv("misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed",
                    col_names = c("seqnames", "start", "end", "score"))
rl_cons <- toGRanges(as.data.frame(rl_cons))




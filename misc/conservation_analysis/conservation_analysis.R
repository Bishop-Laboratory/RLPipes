##### Correlation analysis #####
library(tidyverse)
library(ChIPpeakAnno)

# Get the metadata
sample_metadata <- read_csv("misc/rmap_full_11_25.csv")

# Steps to get the data/hg38_10kb-widow_2.5kb-step_tiles.bed file ##

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
              " -C > misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed")
if (! file.exists('misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed')) {
  system(cmd)
}

# Pivot the resulting table into a matrix
rl_cons_raw <- read_tsv("misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.intersected.bed",
                        col_names = c("seqnames", "start", "end", "intersected_file", "number_of_intersects"))
number_of_peak_files <- 252
rl_cons_mat <- rl_cons_raw %>%
  mutate(id = rep(paste0("window_", seq(length(rl_cons_raw$seqnames) / number_of_peak_files)), each = number_of_peak_files)) %>%
  mutate(number_of_intersects = ifelse(number_of_intersects > 0, 1, 0)) %>%
  mutate(intersected_file = paste0("peakfile_", intersected_file)) %>%
  select(id, intersected_file, number_of_intersects) %>%
  pivot_wider(names_from = intersected_file, values_from = number_of_intersects)
rl_cons_mat2 <- column_to_rownames(rl_cons_mat, var = "id")
rs <- rowSums(rl_cons_mat2)

window_key <- rl_cons_raw %>%
  mutate(id = rep(paste0("window_", seq(length(rl_cons_raw$seqnames) / number_of_peak_files)), each = number_of_peak_files)) %>%
  select(id, seqnames, start, end) %>%
  distinct(id, .keep_all = TRUE)
conservation_simple <- window_key %>%
  mutate(score = rs)
write_tsv(conservation_simple, path = "misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.conservation_levels.bed")

## Steps to do the conservation analysis -- start here ##

# Read in the results as a GRanges object
rl_cons <- read_tsv("misc/conservation_analysis/data/hg38_10kb-widow_2.5kb-step_tiles.conservation_levels.bed", skip = 1,
                    col_names = c("name", "seqnames", "start", "end", "score"))
rl_cons <- toGRanges(as.data.frame(rl_cons))

# The score column of rl_cons represents the number of samples overlapping with that 10kb range
# See the histogram of scores
hist(rl_cons$score)
# Total number of samples is 252
total_number_of_samples <- 252

# Get the pct conservation
rl_cons$pct_cons <- 100 * (rl_cons$score / 252)

# See histogram of pct conservation
hist(rl_cons$pct_cons)

# Some windows have many R-loops overlapping with them, some have very few
top_ranges <- as.data.frame(rl_cons) %>%
  top_n(10, score) %>%
  rownames_to_column() %>%
  pull(rowname)
rl_cons[top_ranges,]  # Top 10 Ranges with the maximum score

# Continue with the rest of the analysis here...






library(tidyverse)
library(readxl)

# Wrangle samples
info_sheet <- read_xlsx("RSeq_CLI/analysis/R_loop_map_accessions_11092020.xlsx")
samples_available <- list.dirs("RSeq_CLI/tests/RSeq_out16/", recursive = FALSE, full.names = F)
rmap_samples <- read_csv("RSeq_CLI/tests/RSeq_out16/rseqVars.csv") %>%
  select(-X1) %>%
  filter(sample_name %in% samples_available) %>%
  rename(GSM = experiment_orig) %>%
  left_join(y = info_sheet, by = "GSM") %>%
  filter(! duplicated(GSM)) %>%
  select(-Type, -Issues, -Reference, -AddInfo, -SpikeIn, -ControlSample, -cores, -file_type,
         -Run, -outdir, -control_orig, -genome_home_dir, -effective_genome_size, -full_genome_length) %>%
  mutate_at(.vars = "Group", str_to_title) %>%
  mutate(ControlType = case_when(ControlType == "NA" ~ "None",
                                 TRUE ~ ControlType))
rmap_samples$Condition[rmap_samples$Condition == "Input"] <- c("S9.6", "S9.6", "RNaseH1", "RNaseH1", "RNaseH1")
rmap_samples$Condition[rmap_samples$Condition %in% c("RNH")] <- "RNaseH1"
rmap_samples$Condition[rmap_samples$Condition %in% c("S96")] <- "S9.6"
write_csv(rmap_samples, "RSeq_CLI/analysis/rmapsamples_10_05_2020.csv")  
rmap_samples <- read_csv("RSeq_CLI/analysis/rmapsamples_10_05_2020.csv")

## Summary stats ##
# Type
rmap_samples %>%
  group_by(mode) %>%
  tally() %>%
  mutate(pct = round(n/sum(n)*100)) %>%
  mutate(label = case_when(pct > 1 ~ as.character(n),
                           TRUE ~ "")) %>%
  arrange(n) %>%
  mutate(mode = factor(mode, levels = mode)) %>%
  ggplot(aes(x = "", y = n, fill = mode)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), size = 5,
            position = position_stack(vjust = 0.5)) +
  labs(title = "RMapDB by Modes") +
  theme_void(base_size = 16) +
  theme(plot.title = element_text(hjust = .5, size = 18))

# Moeity
rmap_samples %>%
  group_by(moeity) %>%
  tally() %>%
  mutate(pct = round(n/sum(n)*100)) %>%
  mutate(label = case_when(pct > 1 ~ as.character(n),
                           TRUE ~ "")) %>%
  arrange(n) %>%
  mutate(moeity = factor(moeity, levels = moeity)) %>%
  ggplot(aes(x = "", y = n, fill = moeity)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), size = 5,
            position = position_stack(vjust = 0.5)) +
  labs(title = "RMapDB by Moeity") +
  theme_void(base_size = 16) +
  theme(plot.title = element_text(hjust = .5, size = 18))

# Genome
rmap_samples %>%
  group_by(genome) %>%
  tally() %>%
  mutate(pct = round(n/sum(n)*100)) %>%
  mutate(label = case_when(pct > 1 ~ as.character(n),
                           TRUE ~ "")) %>%
  arrange(n) %>%
  mutate(genome = factor(genome, levels = genome)) %>%
  ggplot(aes(x = "", y = n, fill = genome)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), size = 5,
            position = position_stack(vjust = 0.5)) +
  labs(title = "RMapDB by Genomes") +
  theme_void(base_size = 16) +
  theme(plot.title = element_text(hjust = .5, size = 18))

# Lib Type
rmap_samples %>%
  group_by(paired_end, strand_specific) %>%
  tally() %>%
  ungroup() %>%
  mutate(pct = round(n/sum(n)*100)) %>%
  mutate(label = case_when(pct > 1 ~ as.character(n),
                           TRUE ~ "")) %>%
  mutate(strand = ifelse(strand_specific, "stranded", "unstranded")) %>%
  mutate(end = ifelse(paired_end, "paired-end", "single-end")) %>%
  mutate(LibType = paste(end, ":", strand)) %>%
  arrange(n) %>%
  mutate(LibType = factor(LibType, levels = LibType)) %>%
  ggplot(aes(x = "", y = n, fill = LibType)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), size = 5,
            position = position_stack(vjust = 0.5)) +
  labs(title = "RMapDB by LibType") +
  theme_void(base_size = 16) +
  theme(plot.title = element_text(hjust = .5, size = 18))


# Group
rmap_samples %>%
  group_by(Group) %>%
  tally() %>%
  mutate(pct = round(n/sum(n)*100)) %>%
  mutate(label = case_when(pct > 1 ~ as.character(n),
                           TRUE ~ "")) %>%
  arrange(n) %>%
  mutate(Group = factor(Group, levels = Group)) %>%
  ggplot(aes(x = "", y = n, fill = Group)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), size = 5,
            position = position_stack(vjust = 0.5)) +
  labs(title = "RMapDB by Group") +
  theme_void(base_size = 16) +
  theme(plot.title = element_text(hjust = .5, size = 18))


# Control Type
rmap_samples %>%
  group_by(ControlType) %>%
  tally() %>%
  mutate(pct = round(n/sum(n)*100)) %>%
  mutate(label = case_when(pct > 1 ~ as.character(n),
                           TRUE ~ "")) %>%
  arrange(n) %>%
  mutate(ControlType = factor(ControlType, levels = ControlType)) %>%
  ggplot(aes(x = "", y = n, fill = ControlType)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), size = 5,
            position = position_stack(vjust = 0.5)) +
  labs(title = "RMapDB by ControlType") +
  theme_void(base_size = 16) +
  theme(plot.title = element_text(hjust = .5, size = 18))



## bin the unstranded genome into length-10000 tiles ##
# These are the large regions which might ever have an R - loop :: then we can narrow later...
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)


#bedtools makewindows -g RSeq_CLI/analysis/hg38.chrom.sizes -w 100 -i srcwinnum > "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/analysis/100bp_tiles_hg38.bed"
#bedtools intersect -a "RSeq_CLI/analysis/100bp_tiles_hg38.bed" -b RSeq_CLI/tests/RSeq_out16/*/peaks_macs_unstranded/*hg38*.broadPeak -c > RSeq_CLI/analysis/100bp_tiles_hg38_macs_intersected.bed
#bedtools intersect -a "RSeq_CLI/analysis/100bp_tiles_hg38.bed" -b RSeq_CLI/tests/RSeq_out16/*/peaks_epic_unstranded/*hg38.bed -c > RSeq_CLI/analysis/100bp_tiles_hg38_epic_intersected.bed
#awk '{ print $1"\t"$2"\t"$3"\t"$5 }' RSeq_CLI/analysis/100bp_tiles_hg38_macs_intersected.bed > RSeq_CLI/analysis/100bp_tiles_hg38_macs_intersected.bedGraph
#awk '{ print $1"\t"$2"\t"$3"\t"$5 }' RSeq_CLI/analysis/100bp_tiles_hg38_epic_intersected.bed > RSeq_CLI/analysis/100bp_tiles_hg38_epic_intersected.bedGraph


fs <- list.files(path = "RSeq_CLI/tests/RSeq_out16/", recursive = TRUE)
fs <- fs[grep(fs, pattern = ".+/peaks_epic_unstranded/.+hg38.bed")]
kb_tiles <- import("RSeq_CLI/analysis/1kb_tiles_hg38_epic_intersected_conserved_only.bed")
kb_tiles <- kb_tiles[kb_tiles$score > 0,]
kb_tiles$conservation_pct <- (kb_tiles$score/length(fs)) * 100
hist(kb_tiles$conservation_pct, breaks = 100, main = "Conservation of R-loop regions", xlab = "Percent conservation")

# Overlaped with gold standard genes
#bedtools intersect -a "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/analysis/1kb_tiles_hg38_epic_intersected.bed" -b "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/analysis/SMRFseq_Regions.hg38.bed" -wa > "/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/RSeq_CLI/analysis/1kb_tiles_hg38_epic_intersected.gold.bed"
kb_tiles_gold <- import("RSeq_CLI/analysis/1kb_tiles_hg38_epic_intersected.gold.bed")
kb_tiles$type <- "Not gold"
kb_tiles$type[kb_tiles$name %in% kb_tiles_gold$name] <- "gold"

library(tidyverse)
library(ggpubr)
as_tibble(mcols(kb_tiles)) %>%
  ggplot(mapping = aes(x = type, y = conservation_pct)) +
  geom_boxplot() +
  ylab("Percent conservation") +
  xlab(NULL) +
  labs(title = "Gold standard R-loop conservation")

# TODO: GC skew
# TODO: RLFS
# TODO: overlap with repeat regions



# Correlation module
# 1. Gold standard genes with 1kb windows
# bedtools intersect -a RSeq_CLI/analysis/1kb_tiles_hg38.bed -b RSeq_CLI/analysis/correlation_genes_100kb.hg38.bed > /RSeq_CLI/analysis/correlation_genes_100kb.hg38.1kbwindow.bed
# 2. Get all of the coverage across the gold-standard genes with 1kb windows
# multiBigwigSummary BED-file --BED RSeq_CLI/analysis/correlation_genes_100kb.hg38.1kbwindow.bed -o res.npz -b RSeq_CLI/tests/RSeq_out16/*/coverage_unstranded/*hg38.bw --outRawCounts RSeq_CLI/analysis/gold_standard_bw_coverage.tab -p 80
# 3. Do the correlation analysis in R
covSum <- read_tsv("RSeq_CLI/analysis/gold_standard_bw_coverage.tab")

mat <- as.matrix(covSum %>%
                   select(-1, -2, -3))
colnames(mat) <- gsub(colnames(mat), pattern = "'([ES]RX[0-9]+_.+).hg38.bw'", replacement = "\\1")
rownames(mat) <- paste0(covSum[,1, drop = TRUE], "_", covSum[,2, drop = TRUE], "_", covSum[,3, drop = TRUE])
library(pheatmap)
set.seed(42)
colData <- rmap_samples %>%
  mutate(newname = paste0(mode, "_", Cell, "_", Group, "_", Condition, "_", gsub(sample_name, pattern = "([ES]RX[0-9]+)_.+", replacement = "\\1"))) %>%
  filter(sample_name %in% colnames(mat)) %>%
  arrange(match(sample_name, colnames(mat))) %>%
  column_to_rownames(var = "newname") %>%
  select(mode)
colnames(mat) <- rownames(colData)
all(rownames(colData) == colnames(mat))
save(mat, colData, file = "RSeq_CLI/helpers/data/gold_standard_bw_coverage.rda")

samples <- sample(colnames(mat), 40)
matsmall <- mat[, samples]
colSmall <- colData[samples,, drop = FALSE]
corSmall <- cor(matsmall)

pheatmap(corSmall, show_rownames = TRUE, show_colnames = TRUE, height = 18, width = 18, 
         annotation_col = colSmall,
         filename = "RSeq_CLI/analysis/heatmap_small.png")

pres <- prcomp(mat)
ggscatter(as.data.frame(pres$rotation), x = "PC1", y = "PC2")



mat <- read_tsv("RSeq_CLI/analysis/all_epic_tiles_hg38_1kb_bins_summary.tab")




# # Okay... so this is going to be a rough one. If you're reading this, good luck.
# # The problem: we needed to change some axis labels to different colors...
# # SO says it's possible but difficult: https://stackoverflow.com/questions/57486547/coloring-the-axis-tick-text-by-multiple-colors/58643916#58643916
# # That solution didn't actually work here because we already have an x1 and x2 axis set.
# # We needed to make an x3 axis. However, using add_trace(xaxis="x3", showscale=F) gives an error.
# # Oh no!
# # So... the solution is to basically:
# # 1. copy the "layout" for xaxis into a new list entry in layout "xaxis3"
# # 2. set the indices to change color for ticktext and tickvals
# # 3. Change the color using the tickfont param for the next xaxis3
# # 4. Set the heatmap data layer (layer 5) to map to "x3" which is xaxis3.
# # The reason step #4 is necessary is because it's not enough to simply make xaxis3 layout -- there has to be data associated with it or it wont be drawn. We can associate the heatmap layer data to xaxis3 in order to fulfill this requirement. Email Henry at millerh1@uthscsa.edu if you want any more help with this...
# hp <- heatmaply_cor(x = corSmall, k_row = NA, k_col =  NA)
# hp$height <- 1200
# hp$width <- 1500
# hp$x$layout$xaxis3 <- hp$x$layout$xaxis  # step #1
# hp$x$layout$xaxis3$ticktext <- hp$x$layout$xaxis3$ticktext[1:5]  # Step #2
# hp$x$layout$xaxis3$tickvals <- hp$x$layout$xaxis3$tickvals[1:5]
# hp$x$layout$xaxis3$tickfont$color <- "green"  # Step #3
# hp$x$data[[10]] <- hp$x$data[[5]]
# hp$x$data[[10]]$xaxis <- "x3"  # Step #4
# # Rinse and repeat for the y axis (it's now y2 instead of y)
# hp$x$layout$yaxis3 <- hp$x$layout$yaxis2  # step #1
# hp$x$layout$yaxis3$overlapping <- "y2"
# hp$x$layout$yaxis3$ticktext <- hp$x$layout$yaxis3$ticktext[1:5]  # Step #2
# hp$x$layout$yaxis3$tickvals <- hp$x$layout$yaxis3$tickvals[1:5]
# hp$x$layout$yaxis3$tickfont$color <- "green"  # Step #3
# hp$x$data[[5]]$yaxis <- "y3"  # Step #4
# 
# hp %>% add_heatmap(z = corSmall, x = 1:length(rownames(corSmall)))


















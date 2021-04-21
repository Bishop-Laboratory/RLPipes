library(tidyverse)
library(pheatmap)

# Correlation module
analysis_folder <- "analysis/ext_data/"
rmapdb_folder <- "data/"
# 1. Make 1kb tiles of hg38
cmd <- paste0("wget -O ", analysis_folder, "/hg38.chrom.sizes http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes && ", 
              "bedtools makewindows -g ", analysis_folder, "/hg38.chrom.sizes -w 1000 -i srcwinnum > ", 
              analysis_folder, "1kb_tiles_hg38.bed")
system(cmd)
# 2. Gold standard genes with 1kb windows
cmd <- paste0("bedtools intersect -a ", analysis_folder, "1kb_tiles_hg38.bed -b ", analysis_folder,
              "correlation_genes_100kb.hg38.bed > ", analysis_folder,
              "correlation_genes_100kb.hg38.1kbwindow.bed")
system(cmd)
# 2. Get all of the coverage across the gold-standard genes with 1kb windows
cmd <- paste0("multiBigwigSummary BED-file --BED ", analysis_folder, 
              "/correlation_genes_100kb.hg38.1kbwindow.bed -o ", analysis_folder, 
              "/corr_gene_bw_res.npz -b ", rmapdb_folder, 
              "/*/coverage_unstranded/*hg38.bw --outRawCounts ", analysis_folder, 
              "/gold_standard_bw_coverage.tab -p 80")
system(cmd)
# 3. Do the correlation analysis in R
covSum <- read_tsv(paste0(analysis_folder, "/gold_standard_bw_coverage.tab"))
mat <- as.matrix(covSum %>%
                   select(-1, -2, -3))
colnames(mat) <- gsub(colnames(mat), pattern = "'([ES]RX[0-9]+_.+).hg38.bw'", replacement = "\\1")
rownames(mat) <- paste0(covSum[,1, drop = TRUE], "_", covSum[,2, drop = TRUE], "_", covSum[,3, drop = TRUE])
colData <- rmap_samples %>%
  filter(sample_name %in% colnames(mat)) %>%
  arrange(match(sample_name, colnames(mat))) %>%
  column_to_rownames(var = "clean_name") %>%
  select(mode)
colnames(mat) <- rownames(colData)
all(rownames(colData) == colnames(mat))
save(mat, colData, file = "../../helpers/data/gold_standard_bw_coverage.rda")

corrmat <- cor(mat)
pheatmap(corrmat, annotation_col = colData, 
         filename = "results/full_corr_heatmap.png", fontsize = 6,
         height = 18, width = 20)





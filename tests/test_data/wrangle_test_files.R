# Script for wrangling test files to correct naming conventions
library(tidyverse)
library(parallel)

MC_CORES <- 44
S3_BAM_URI <- "s3://rseq-testing/bam-files/"
bamsAvail <- system(paste0("aws s3 ls ", S3_BAM_URI), intern = TRUE)

oldnew <- tibble(
  oldfls = gsub(bamsAvail, pattern = ".+ ([ES]{1}RX[0-9]+_.+\\.[hgmm]{2}[0-9]+\\.bam)", replacement = "\\1"),
  newfls = gsub(bamsAvail, pattern = ".+ ([ES]{1}RX[0-9]+)_.+\\.([hgmm]{2}[0-9]+\\.bam)", replacement = "\\1_\\2")
) %>%
  mutate(across(everything(), function(x) {paste0(S3_BAM_URI, x)})) 

mclapply(seq(rownames(oldnew)), function(i) {
  old <- oldnew$oldfls[i]
  new <- oldnew$newfls[i]
  system(paste0("aws s3 mv ", old, " ", new))
}, mc.cores = MC_CORES)

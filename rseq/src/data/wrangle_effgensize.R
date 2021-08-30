library(tidyverse)
av_gen <- read_tsv("available_genomes.tsv.xz")
av_gen %>%
  select(UCSC_orgID, contains("eff_gen")) %>%
  pivot_longer(cols = contains('eff_gen')) %>%
  mutate(read_length = as.numeric(gsub(name, pattern = ".+_([0-9]+)bp$", replacement = "\\1"))) %>%
  fill(value, .direction = "down") %>%
  select(-name, eff_genome_size=value) %>%
  relocate(read_length, .before=eff_genome_size) %>%
  write_tsv("eff_gen_size.tsv")
system("xz -f eff_gen_size.tsv")

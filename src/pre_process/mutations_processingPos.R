library(tidyverse)

source("./src/utils/read_hgvsp.R")


# mutation data processing


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "protein")


mut <- read_tsv(file = "./output/files/mutations.txt.gz") %>%
  inner_join(samples[, c("sample", "batch")], by = "sample") %>%
  select(sample, batch, everything()) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, sample) %>%
  nest() %>%
  ungroup() %>%
  mutate(pos = map(.x = data, .f = read_hgvsp)) %>%
  unnest() %>%
  mutate(prot_pos = as.numeric(prot_pos)) %>%
  mutate(same = replace_na(same, 1)) %>%
  filter(same == 1) %>%
  select(-n, -batch, -same, -gene_id, -transcript_id, -protein_id)


gz1 <- gzfile("./output/files/mutations_protpos.txt.gz", "w")
write.table(mut, gz1, sep="\t", quote=F, row.names=F)
close(gz1)


library(tidyverse)


# construction of a mutation matrix
mut <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")


# remove splice site mutations
mut <- mut %>%
  filter(variant_class != "Splice_Site")


mat <- mut %>%
  group_by(sample, gene_symbol) %>%
  summarise() %>%
  ungroup() %>%
  mutate(state = 1) %>%
  group_by(gene_symbol) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = "sample", values_from = "state", values_fill = list(state = 0)) %>%
  select(n, gene=gene_symbol, everything())

write_tsv(mat, "./output/files/mutations_matrix.txt.gz")

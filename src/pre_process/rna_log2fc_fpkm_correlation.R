library(tidyverse)

source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("protein", "mRNA"))


fpkm <- read_tsv(file = "./output/files/transcriptomics_log2fpkm.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "fpkm")

log2fc <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fc")

rna <- inner_join(fpkm, log2fc, by = c("gene", "sample"))


corr <- rna %>%
  group_by(gene) %>%
  summarise(
    p = cor(fpkm, log2fc, method = "pearson", use = "p"),
    s = cor(fpkm, log2fc, method = "spearman", use = "p"),
    n = sum(!(is.na(fpkm) | is.na(log2fc)))) %>%
  ungroup()

histogram <- corr %>%
  pivot_longer(-c(gene, n), names_to = "correlation", values_to = "value") %>%
  mutate(correlation = if_else(correlation == "p", "pearson", "spearman")) %>%
  ggplot(mapping = aes(x = value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~ correlation)
histogram


corr <- rna %>%
  inner_join(samples_annotation[, c("sample", "batch", "tissue")], by = "sample") %>%
  group_by(batch, tissue, gene) %>%
  summarise(
    p = cor(fpkm, log2fc, method = "pearson", use = "p"),
    s = cor(fpkm, log2fc, method = "spearman", use = "p")) %>%
  ungroup() %>%
  unite(col = "batch", "batch", "tissue", sep = "-")

histogram <- corr %>%
  pivot_longer(-c(gene, batch), names_to = "correlation", values_to = "value") %>%
  mutate(correlation = if_else(correlation == "p", "pearson", "spearman")) %>%
  ggplot(mapping = aes(x = value)) +
  geom_histogram(bins = 10) +
  facet_grid(batch ~ correlation)
histogram



protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "protein_log2fc")

protein_rna <- rna %>%
  rename(rna_fpkm = fpkm, rna_log2fc = log2fc) %>%
  inner_join(protein, by = c("gene", "sample"))

corr_protein_rna <- protein_rna %>%
  group_by(gene) %>%
  summarise(
    p_fpkm = cor(rna_fpkm, protein_log2fc, method = "pearson", use = "p"),
    s_fpkm = cor(rna_fpkm, protein_log2fc, method = "spearman", use = "p"),
    p_log2fc = cor(rna_log2fc, protein_log2fc, method = "pearson", use = "p"),
    s_log2fc = cor(rna_log2fc, protein_log2fc, method = "spearman", use = "p")) %>%
  ungroup()


histogram <- corr_protein_rna %>%
  pivot_longer(-gene, names_to = c("correlation", "type"), names_sep = "_", values_to = "value") %>%
  mutate(correlation = if_else(correlation == "p", "pearson", "spearman")) %>%
  ggplot(mapping = aes(x = value)) +
  geom_histogram(bins = 10) +
  #geom_density(fill = "grey", color = "black") +
  facet_grid(rows = vars(type), cols = vars(correlation), scales = "fixed")
histogram




corr_protein_rna <- protein_rna %>%
  inner_join(samples_annotation[, c("sample", "batch", "tissue")], by = "sample") %>%
  group_by(batch, tissue, gene) %>%
  summarise(
    p_fpkm = cor(rna_fpkm, protein_log2fc, method = "pearson", use = "p"),
    s_fpkm = cor(rna_fpkm, protein_log2fc, method = "spearman", use = "p"),
    p_log2fc = cor(rna_log2fc, protein_log2fc, method = "pearson", use = "p"),
    s_log2fc = cor(rna_log2fc, protein_log2fc, method = "spearman", use = "p")) %>%
  ungroup() %>%
  unite(col = "batch", "batch", "tissue", sep = "-")


histogram <- corr_protein_rna %>%
  select(-contains("fpkm")) %>%
  pivot_longer(-c(gene, batch), names_to = "correlation", values_to = "value") %>%
  mutate(correlation = if_else(correlation == "p_log2fc", "pearson", "spearman")) %>%
  ggplot(mapping = aes(x = value)) +
  geom_histogram(bins = 10) +
  #geom_density(fill = "grey", color = "black") +
  facet_grid(rows = vars(batch), cols = vars(correlation), scales = "free")
histogram



# load R packages
library(tidyverse)


# source R functions
source("./src/utils/psite_prot_seq_match.R")


# load protein sequences for the canonical uniprot transcripts
gene_prot_seq <- read_tsv(file = "./output/files/uniprot_canTrpt_protSeqGRCh37.txt")


# load compiled data
load(file = "./data/phospho/esetNR.Rdata")

# condition annotation
cond_anno <- pData(esetNR)
cond_anno <- as_tibble(cond_anno)

# phosphopeptide annotation
phospho_anno <- fData(esetNR)
phospho_anno <- as_tibble(phospho_anno)

# phosphorylation data (fold-changes)
phospho_vals <- exprs(esetNR)
phospho_vals <- as_tibble(phospho_vals)


# remove CPTAC samples from dataset
cond_anno <- cond_anno %>%
  filter(!str_detect(biological_sample, "CPTAC")) %>%
  mutate(sample = str_c(condition_id, experiment_id, sep="_")) %>%
  select(sample, everything())

phospho_vals <- phospho_vals %>%
  select_if(.predicate = colnames(.) %in% cond_anno$sample)


# phospho-sites
bk_pho <- phospho_anno %>%
  select(gene_name, position=positions, residue=residues) %>%
  bind_cols(phospho_vals)

duplicated <- bk_pho %>%
  select(1:3) %>%
  filter(duplicated(paste(gene_name, position, residue, sep = "_")))

bk_pho <- bk_pho %>%
  anti_join(duplicated, by=c("gene_name", "position", "residue")) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, gene=gene_name, everything()) %>%
  group_by(n, gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite2)) %>%
  filter(same == 1) %>%
  unnest(cols = data) %>%
  select(-n, -seq, -same)


# remove phosphosites without expression in all samples
bk_pho <- bk_pho %>%
  pivot_longer(-c(gene, position, residue), names_to = "sample", values_to = "log2fc") %>%
  group_by(gene, position, residue) %>%
  filter(!(sum(is.na(log2fc)) == n())) %>%
  ungroup() %>%
  pivot_wider(names_from = "sample", values_from = "log2fc")


write_tsv(cond_anno, "./output/files/esetNR_cond_anno.txt")
write_tsv(bk_pho, "./output/files/esetNR_phospho.txt.gz")

# load R packages
library(tidyverse)


# source R functions
source("./src/utils/psite_prot_seq_match.R")


# load protein sequences for the canonical uniprot transcripts
gene_prot_seq <- read_tsv(file = "./output/files/uniprot_canTrpt_protSeqGRCh37.txt")


# load benchmark data

# fold-changes
bk_fc <- read_tsv(file = "./data/kinase_substrate/benchmark_dataset/foldchanges.txt")

# phospho-sites
bk_pho <- read_tsv(file = "./data/kinase_substrate/benchmark_dataset/phosphosites.txt") %>%
  select(-types, -locscores, -ensg, -peptide, -ensp) %>%
  select(gene_name, position=positions, residue=residues) %>%
  bind_cols(bk_fc)

duplicated <- bk_pho %>%
  select(1:3) %>%
  filter(duplicated(paste(gene_name,position,residue, sep = "_")))

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


gz1 <- gzfile("./output/files/bk_pho.txt.gz", "w")
write.table(bk_pho, gz1, sep="\t", quote=F, row.names=F)
close(gz1)

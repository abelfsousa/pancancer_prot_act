library(tidyverse)


# assembling of cnv data from "Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma"
# log2 copy-number ratios from cptac portal


# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(str_detect(batch, "ccrcc"))


# protein coding genes
protein_coding <- read_tsv("./output/files/gencode.v19.annotation.txt") %>%
  select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  select(-gene_type) %>%
  group_by(gene_name) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


ccrcc <- read_tsv(file = "./data/dna/cnv/cptac/ccrcc/all_thresholded.by_genes.txt") %>%
  rename(gene=`Gene Symbol`) %>%
  select(-`Locus ID`, -`Cytoband`) %>%
  rename_all(.funs = ~ str_replace(.x, "-", ".")) %>%
  filter(gene %in% protein_coding$gene_name)


# samples
samples <- tibble(sample = colnames(ccrcc)[-c(1)], batch = "discovery-ccrcc", cancer = "ccrcc")
write.table(samples, "./output/files/cnv_samples_ccrcc.txt", sep="\t", quote=F, row.names=F)


ccrcc <- ccrcc %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

gz1 <- gzfile("./output/files/cnv_ccrcc.txt.gz", "w")
write.table(ccrcc, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
